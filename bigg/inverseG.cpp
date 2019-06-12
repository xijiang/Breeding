// icpc -std=c++11 -mkl -DMKL_ILP64 inverseG.cpp -o inverseG
#include <iostream>
#include <cmath>
#include <vector>
#include <mkl.h>

using namespace std;

void trick_init_G(vector<double>&G, long nid){
  vector<double> t(nid, 0);
  for(long i{0}; i<nid; ++i) G.insert(G.end(), t.begin(), t.end());
}

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  double diag{0.};
  if(argc==2) diag = stof(argv[1]);

  // Read G which is in binary storage
  cin.seekg(0, cin.end);
  MKL_INT sz  = cin.tellg();
  MKL_INT len = sz/sizeof(double);
  MKL_INT dim(sqrt(len)), info;
  cin.seekg(0, cin.beg);
  
  vector<double> G;
  trick_init_G(G, dim);
  
  cin.read(reinterpret_cast<char*>(&G[0]), sz);
  for(MKL_INT i=0; i<dim; ++i) G[i*dim+i]+=diag;
  
  info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', dim, &G[0], dim);
  if(info){
    cerr<<"Decomposition error: "<<info<<endl;
    return 1;
  }
  info = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'L', dim, &G[0], dim);

  if(!info) cout.write(reinterpret_cast<char*>(&G[0]), sz);
  else{
    cerr<<"ERROR: "<<info<<" !!!"<<endl;
    return 2;
  }
  return 0;
}
