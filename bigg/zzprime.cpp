#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <mkl.h>

using namespace std;
void read_genotype(ifstream&in, vector<double>&Z, vector<string>&ID, const vector<double>&twop){
  for(string id, line; in>>id; ID.push_back(id)){
    getline(in, line);
    stringstream ss(line);
    double       td;
    for(auto&x:twop){
      ss>>td;
      Z.push_back(td-x);
    }
  }
}


int main(int argc, char *argv[])
{				// Read vector of 2p
  ios_base::sync_with_stdio(false);
  vector<double> twop;
  double         s2pq{0.}, is2pq, alpha{1.}, beta{0.};
  {
    ifstream fin("2p.txt");
    for(double td; fin>>td; twop.push_back(td)) s2pq+=td*(1-td*.5);
  }
  is2pq = 1./s2pq;		// product is faster than divide
  cout.precision(12);

  vector<double> A;
  vector<string> I;
  ifstream fin(argv[1]);
  
  read_genotype(fin, A, I, twop);
  int m = static_cast<int>(I.size());
  int k = static_cast<int>(A.size()/I.size());
  if(argc==2){
    double G[m*m];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, k, alpha, &A[0], k, &A[0], k, beta, G, m);
    for(auto&x:G) x*=is2pq;
    cout.write(reinterpret_cast<char*>(G), sizeof(double)*m*m);
    //for(int i=0; i<m; ++i) for(int j=0; j<=i; ++j) cout<<I[i]<<' '<<I[j]<<' '<<G[i*m+j]<<'\n';
  }else if(argc==3){
    vector<double> B;
    vector<string> J;
    ifstream f2(argv[2]);
    read_genotype(f2, B, J, twop);
    int n = static_cast<int>(J.size());
    double G[m*n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, &A[0], k, &B[0], k, beta, G, n);
    for(auto&x:G) x*=is2pq;
    cout.write(reinterpret_cast<char*>(G), sizeof(double)*m*n);
    //auto p=G;
    //for(auto&i:I) for(auto&j:J) cout<<i<<' '<<j<<' '<<*p++<<'\n';
  }
  
  return 0;
}
