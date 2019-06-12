#include <iostream>
#include <fstream>
#include <vector>
/**
 * Note: one can declare G as static,
 *  with a predefined size,
 *  and compile with option -mcmodel=large
 */
using namespace std;

void readdiag(double *G, long k, long dim, long bsz, long nid){
  long     len(dim*dim);
  long     beg(k*bsz*nid+k*bsz);	// 1st element to be updated in G
  double   g[len];
  long     i, j;
  ifstream fin(to_string(k+1)+'-'+to_string(k+1)+".G");

  fin.read(reinterpret_cast<char*>(g), len*sizeof(double));
  for(i=0; i<dim; ++i){
    for(j=0; j<=i; ++j) G[beg+j] = g[i*dim+j];
    beg+=nid;
  }
}


void readoffd(double *G, long m, long n, long rws, long bsz, long nid){
  long     len(rws*bsz);		// rows x columns
  long     beg(m*bsz*nid+n*bsz);
  double   g[len];
  long     i, j;
  ifstream fin(to_string(m+1)+'-'+to_string(n+1)+".G");

  fin.read(reinterpret_cast<char*>(g), len*sizeof(double));
  for(i=0; i<rws; ++i){
    for(j=0; j<bsz; ++j) G[beg+j] = g[i*bsz+j];
    beg+=nid;
  }
}


void trick_init_G(vector<double>&G, long nid){
  vector<double> t(nid, 0);
  for(long i{0}; i<nid; ++i) G.insert(G.end(), t.begin(), t.end());
}


int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  long nid(stoi(argv[1]));	// Number of ID, e.g., 63721
  long bsz(stoi(argv[2]));	// block size,   e.g.,  1000
  long len(nid*nid);		// final G size
  long nbk(nid/bsz);		// number of blocks of size bsz
  long nrb(nid%bsz), i, j;	// the residue block size
  vector<double> G; // I only update the lower triangle, and output that also

  trick_init_G(G, nid);		// to dynamite allocate G

  // -- for those of dimension BLK x BLK
  for(i=0; i<nbk; ++i){
    for(j=0; j<i; ++j) readoffd(&G[0], i, j, bsz, bsz, nid);
    readdiag(&G[0], i, bsz, bsz, nid);
  }
  // -- for those of dimension lesser, e.g., at the lower edge
  if(nrb){
    for(i=0; i<nbk; ++i) readoffd(&G[0], nbk, i, nrb, bsz, nid);
    readdiag(&G[0], nbk, nrb, bsz, nid);
  }

  cout.write(reinterpret_cast<char*>(&G[0]), len*sizeof(double));

  return 0;
}
