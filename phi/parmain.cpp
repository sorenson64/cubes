#include<iostream>
#include<fstream>
#include<bitset>
#include<cmath>
#include<gmp.h>

#include"mpi.h"

using namespace std;

const int64_t minintervalsize=128;
//const int64_t numintervals=1024; 
const int64_t numintervals=1024; 
const int64_t maskwidth=minintervalsize*numintervals;

#include "int.h"
#include "smallprimes.h"
#include "rlist.h"
#include "masks.h"
#include "primetest.h"

// for data reporting

uint64_t countprimesfound;
//uint64_t countprimetests; // this is declared in primetest.h
uint64_t countfails;
uint64_t countR[5];
uint64_t countfilter, countfix, countfix2;
uint64_t countT;
long double totalT;
uint64_t countk, totalk;
uint64_t nbegin, nend;

void countreset()
{
  countprimesfound=0;
  countprimetests=0;
  countfails=0;
  countR[0]= countR[1]= countR[2]= countR[3]= countR[4]=0;
  countfilter=countfix=countfix2=0;
  countT=0;
  totalT=0;
  countk=totalk=0;
}

void countprint(ostream & os, int id)
{
  os << nbegin <<" "<< nend 
     << " \t"<<countprimesfound<<" "<<countprimetests<<" "<<countfails
     << " \t"<<countR[0]<<" "<<countR[1]<<" "<<countR[2]
     <<" "<<countR[3]<<" "<<countR[4]
     << " \t"<<countfilter<<" "<<countfix<<" "<<countfix2
     << " \t"<<(((long double)totalT)/countT)
     << " \t"<<(((long double)totalk)/countk)
     << " \t"<<id<< endl;
}

int64_t findR(double n, int col)
{
  long double bound=powl(n,1.0L/3.0L);
  int i=0;
  while(i<rlistlen && rlist[i++][col]<=bound) ;
  return rlist[i][col];
}

int64_t findR(double n) { return findR(n,0); }

//int64_t findstart(int p, int128 M, int128 qleft)
int64_t findstart(int p, int128 M, mpz_t qleft)
{
  int Mmodp=M%p;
  //int qmodp=qleft%p;
  int qmodp=mpz_fdiv_ui(qleft,p);
  int start=p-(inv(Mmodp,p)+qmodp)%p;
  if(start==p) start=0;
  return start;
}

bool filter(int128 n)
{
  countfilter++;
  // trial division by small primes vi GCD
  // primes up thru 7 handled in modulus mR
  const int64_t m=11*13*17*19*23;
  int64_t rem=n%m;
  return (gcd(m,rem)==1);
}

mpz_t N1,N2;

bool fixit2(mpz_t nstart, mpz_t nstop)
{
  countfix2++;
  //bigint2mpz(N1,nstart);
  //bigint2mpz(N2,nstop);
  mpz_set(N1,nstart);
  mpz_set(N2,nstop);

  while(mpz_cmp(N1,N2)<0) // while N1<N2
  {
    mpz_nextprime(N1,N1);
    if(mpz_cmp(N1,N2)>0)
    {
      return false;
    }
    int val=mpz_probab_prime_p(N1,64);
    if(val==2) 
    {
      return true;
    }
    if(val==1) // need to prove primality later
    {
      ofstream f("primes2check.txt",ios::app);
      f << N1 << endl;
      f.close();
      return true;
    }
  }
  return false;
}

mpz_t qstartf, qstopf, qdeltaf, nf;

bool fixit(mpz_t nleft, mpz_t nright)
{
  countfix++;
  int64_t m=2*3*5*7; // good enough here
  for(int r=1; r<5; r++)
  {
    int64_t R=findR(mpz_get_d(nright),r);

    //int128 qstart=nleft/(m*R)+1;
    mpz_set(qstartf,nleft);
    mpz_fdiv_q_ui(qstartf,qstartf,m);
    mpz_fdiv_q_ui(qstartf,qstartf,R);
    mpz_add_ui(qstartf,qstartf,1);

    //int128 qstop=nright/(m*R);
    mpz_set(qstopf,nright);
    mpz_fdiv_q_ui(qstopf,qstopf,m);
    mpz_fdiv_q_ui(qstopf,qstopf,R);

    mpz_sub(qdeltaf,qstopf,qstartf);

    int64_t qdelta=mpz_get_ui(qdeltaf);
    //int128 n=qstart*m*R+1;
    mpz_mul_ui(nf,qstartf,m);
    mpz_mul_ui(nf,nf,R);
    mpz_add_ui(nf,nf,1);

    for(int i=0; i<qdelta; i++)
    {
      //if(filter(n))
      {
        if(primetest(nf,R)) { countR[r]++; return true; }
      }
      //n=n+m*R;
      mpz_add_ui(nf,nf,m*R);
    }
  }
  return fixit2(nleft,nright);
}

bool fixit_nogmp(int128 nleft, int128 nright)
{
  countfix++;
  int64_t m=2*3*5*7*11*13; // good enough here
  for(int r=1; r<5; r++)
  {
    int64_t R=findR((double)nright,r);
    int128 qstart=nleft/(m*R)+1;
    int128 qstop=nright/(m*R);
    int64_t qdelta=qstop-qstart;
    int128 n=qstart*m*R+1;
    for(int i=0; i<qdelta; i++)
    {
      if(filter(n))
      {
        if(primetest(n,R)) { countR[r]++; return true; }
      }
      n=n+m*R;
    }
  }
  //return fixit2(nleft,nright);
  return 1;
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);

  //cout << atoi(argv[1]) << endl;

  // get our name
  int id;
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  // how many of us are there?
  int np;
  MPI_Comm_size(MPI_COMM_WORLD,&np);

/*  Phi version - node numbers not needed
  int thisnode=atoi(argv[1]);

  np=np*4; // should be 560 but node-4 is lame
  id=id+thisnode*112; // adjust id number too
*/
 
  char filename[]="data/xxx.out";
  filename[5]=(char)(id/100+'0');
  filename[6]=(char)((id/10)%10+'0');
  filename[7]=(char)(id%10+'0');

  ofstream fout(filename,ios::app);

  //makemasks();  // no masks for cubes
  primetestinits();
  countreset();
  Smallprimes p(maskwidth);
  bitset<maskwidth> x; // bit "vector" for sieving
  //printmasks();
  //
  mpz_init(N1); mpz_init(N2); 
  mpz_t ncubed; mpz_init(ncubed);
  mpz_t nplustcubed; mpz_init(nplustcubed);
  mpz_t Mgmp; mpz_init(Mgmp);
  mpz_t QL; mpz_init(QL);
  mpz_t QR; mpz_init(QR);
  mpz_t nstart, nstop, qstart, qstop, jstart, jstop, N;
  mpz_init(nstart); mpz_init(nstop);
  mpz_init(qstart); mpz_init(qstop);
  mpz_init(jstart); mpz_init(jstop);
  mpz_init(N);
  mpz_init(qstartf); mpz_init(qstopf); mpz_init(qdeltaf); mpz_init(nf);

  int64_t loopcount=0;

  //while(true) // BIG LOOP on n
  //while(n<2e9+1e7) // BIG LOOP on n
  //while(n<1e10+1e7) // BIG LOOP on n

const uint64_t nchunk=1e7;
//for(uint64_t myn=2.244e13+id*nchunk; myn<3e13; myn+=np*nchunk)
for(uint64_t myn=2.54e13+id*nchunk; myn<3e13; myn+=np*nchunk)
{
  uint64_t n=myn; // start point
  nbegin=n;

  while(n<myn+nchunk)
{
  loopcount++;

  //int128 ncubed=((int128)n)*((int128)n)*((int128)n);
  mpz_ui_pow_ui(ncubed,n,3);

  //const int t=numintervals;
  int64_t T=numintervals;

  //int128 nplustcubed=((int128)n+T)*((int128)n+T)*((int128)n+T);
  mpz_ui_pow_ui(nplustcubed,n+T,3);

  //int128 deltan=nplustcubed-ncubed;
  int128 deltan=3*((int128)n)*n*T+3*((int128)n)*T*T+T*T*T;

  //int64_t R=findR(((int128)(n+t))*(n+t)); // need R>(n+t)^{2/3}
  int64_t R=findR(mpz_get_d(nplustcubed)); // need R>(n+t)

  // compute the modulus M=R*2*small primes
  int128 M=2*R;
  int k=1;

  while(deltan/(T*M*primes[k])>minintervalsize) { M *= primes[k++]; }

  //cout << "M=" << M << endl;

  //while(n/(2*M)>minintervalsize) { M *= 2; }
  int64_t adjust=(deltan/(T*M))/minintervalsize;
  if(adjust>1) M *= adjust;

  //cout << "M=" << M << endl;

  bigint2mpz(Mgmp,M);

  totalk+=k; countk++;

  // compute T
  int64_t intervalsize=deltan/(T*M);

  T=maskwidth/(intervalsize+1);
  //nplustcubed=((int128)n+T)*((int128)n+T)*((int128)n+T);
  //deltan=nplustcubed-ncubed;

  mpz_ui_pow_ui(nplustcubed,n+T,3);
  deltan=3*((int128)n)*n*T+3*((int128)n)*T*T+T*T*T;

  totalT+=T; countT++;

  //T=numintervals;

  //cout << "  T="<<T<<" numintervals="<<numintervals<<endl;

  // numbers are of the form M*q+1, from q=left to right
  // so qleft is bit position 0 of the bitset x
  //int128 qleft=(ncubed/M)+1; // n cubed is not prime
  mpz_fdiv_q(QL,ncubed,Mgmp);
  mpz_add_ui(QL,QL,1);

  //int128 qright=(nplustcubed/M);
  mpz_fdiv_q(QR,nplustcubed,Mgmp);

  if((n+T)%M==0) mpz_sub_ui(QR,QR,1); //qright--;

  int64_t qdelta; //=qright-qleft+1;
  mpz_sub(save,QR,QL);
  qdelta=mpz_get_ui(save);
  qdelta++;

  /* * /
  cout << "n=" << n << " R="<<R<<" k="<<k<<" M="<<M
    <<" intv size="<< deltan/M/T
	  <<endl;
  cout << "target n="<<n+T
	  << " each interval averages "<< (deltan)/(T*M) 
	  << " total size=" << (deltan)/M 
	  << endl;
  cout << "qleft=" << (QL) // qleft
	  << " qright=" <<  (QR) //qright 
	  << " qdelta=" << qdelta  // qright-qleft+1
	  << endl;
  cout << "maskwidth="<<maskwidth<<endl;
  / * */

  x.reset(); // all bits are zero
  // do masks first -- NO no masking for cubes - too many primes in M
/* * /
  for(int i=1; i<masklen; i++) // not using k here on purpose
    if(M%primes[i]!=0) // don't do primes in the modulus!
      {
	int p=primes[i];
        int start=findstart(p,M,qleft);
	// check that (M*(qleft+start)+1 is divisible by p
	//if( ((M%p)*(qleft+start)%p +1 )%p != 0) 
	//{ cout << "Error! "<< p << endl; break; }
	y=mask[i];
	y <<= start;
	x |= y;
	/ *
	cout << primes[i] << endl;
	cout << "---" << endl;
	cout << mask[i] << endl;
	cout << "---" << endl;
	cout << (mask[i]<<start) << endl;
	cout << "---" << endl;
	cout << x << endl;
	cout << endl;
      }
	* /
/ * */
  for(int i=1; i<p.length(); i++)
  if(M%p[i]!=0)
  {
    int start= findstart(p[i],M,QL);
    //if( ((M%p[i])*(qleft+start)%p[i] +1 )%p[i] != 0) 
	//{ cout << "Error! "<< p[i] << endl; break; }
    for(int j=start; j<qdelta; j+=p[i]) x.set(j);
  }
  //cout << endl << "---" << endl;
  //cout << x << endl;



  for(int t=0; t<T; t++)
  {
    // left half (A)  [(n+t)^2,(n+t)(n+t+1)]
    // [(n+t)^3,(n+t+1)^3]
    //int128 nstart=((int128)(n+t))*(n+t)*(n+t);
    mpz_ui_pow_ui(nstart,n+t,3);
    //int128 qstart=nstart/M+1;
    mpz_fdiv_q(qstart,nstart,Mgmp);
    mpz_add_ui(qstart,qstart,1);
    //int128 nstop=((int128)(n+t+1))*(n+t+1)*(n+t+1);
    mpz_ui_pow_ui(nstop,n+t+1,3);
    //int128 qstop=nstop/M;
    mpz_fdiv_q(qstop,nstop,Mgmp);

    mpz_sub(jstart,qstart,QL);
    mpz_sub(jstop,qstop,QL);

    bool done=false;
    //for(int j=qstart-qleft;!done && j<qstop-qleft; j++)
    for(int j=mpz_get_ui(jstart);!done && j<mpz_get_ui(jstop); j++)
      if(!x.test(j))
      {
        //int128 n=((int128)M)*(qleft+j)+1;
        mpz_add_ui(save,QL,j);
        mpz_mul(N,save,Mgmp);
        mpz_add_ui(N,N,1);
        if(primetest( N, R )) done=true;
        //if(done) { cout << nstart <<" "<<N<<" "<<nstop<<endl; }
      }


    if(done) countR[0]++;
    if(!done) countfails++;
    if(!done) done=fixit(nstart,nstop);

    if(done) countprimesfound++;

  // NOT NEEDED - only first half is done for cubes
    // right half (B)  [(n+t)(n+t+1),(n+t+1)^2]
/* * /
    nstart=((int128)(n+t))*(n+t+1);
    qstart=nstart/M+1;
    nstop=((int128)(n+t+1))*(n+t+1);
    qstop=nstop/M;

    done=false;
    for(int j=qstart-qleft;!done && j<qstop-qleft; j++)
      if(!x.test(j))
      {
        int128 n=((int128)M)*(qleft+j)+1;
        if(primetest( n, R )) done=true;
      }

    if(done) countR[0]++;
    if(!done) countfails++;
    if(!done) done=fixit(nstart,nstop);

    if(done) countprimesfound++;
/ * */
  }

  /*
  cout << "Total Failure count="<<failcount<<endl;
  cout << "Total Primes found="<<primesfound<<endl;
  cout << "Number of prime tests performed="<<numtests
	  << " success rate="<< (100.0*primesfound)/numtests<<endl;
  cout << "Processed from "<< n <<"^2 to "<< n+numintervals <<"^2"<<endl;
  */

  n+=T;

  if(false && loopcount%10000==0)
  {
  cout << "n="<<n<<endl;
  cout << "Failure count="<<countfails<<endl;
  cout << "Primes found="<<countprimesfound<<endl;
  cout << "Number of prime tests performed="<<countprimetests
	  << " success rate="<< (100.0*countprimesfound)/countprimetests<<endl;
  }
} // BIG loop on n
  nend=n;
  
/*
  cout << "\nFinal n="<<n<<endl;
  cout << "Failure count="<<failcount<<endl;
  cout << "Primes found="<<primesfound<<endl;
  cout << "Number of prime tests performed="<<numtests
	  << " success rate="<< (100.0*primesfound)/numtests<<endl;
*/
  countprint(fout,id);
  fout.flush();
  } // BIG for loop for parallelization

  fout.close();

  MPI_Finalize();

  return 0;
}
