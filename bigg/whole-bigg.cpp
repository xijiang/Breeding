#include <iostream>
#include <bitset>
#include <vector>
#include <sstream>
#include <map>
#include <numeric>		// inner_product
#include <thread>

using namespace std;

// better set it to a number close to you dataset
// Here it is for 1M SNP
const int MLC{2'000'000};
using MBS=bitset<MLC>;
using VBT=vector<MBS>;
using VID=vector<string>;
using VDB=vector<double>;
using MSD=map<string, double>;
using MVD=map<int, VDB>;

MSD msd{{"00", 0.}, {"01", 1.}, {"10", 2.}};

void read_genotype(istream&, VBT&, VID&, VDB&);
void bitset_2_double(const MBS&, double*, int, const VDB&);
void worker(VBT&, VDB&, double, int, int, VDB&);

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  int NT(thread::hardware_concurrency()-2);
  if(NT<1) NT=1;
  
  VBT GT;
  VID ID;
  VDB twop;
  double s2pq{0};
  
  read_genotype(cin, GT, ID, twop);
  for(auto&x:twop) s2pq+=x*(1-x*.5);

  int nid=static_cast<int>(ID.size());
  int nlc=static_cast<int>(twop.size());
  MVD G;
  
  for(auto i=0; i<nid; i+=NT){
    thread th[NT];
    VDB    gi[NT];
    for(auto j=0; j<NT; ++j) if(i+j<nid) th[j] = thread(worker, ref(GT), ref(twop), s2pq, i+j, nlc, ref(gi[j]));
    for(auto j=0; j<NT; ++j) if(i+j<nid){
	th[j].join();
	G[i+j] = gi[j];
      }
  }

  cout<<nid<<'\n';
  for(auto&id:ID) cout<<id<<'\n';
  for(auto&[i, p]:G){
    for(auto&q:p) cout.write(reinterpret_cast<char*>(&q), sizeof(double));
    cout<<'\n';
  }
  
  return 0;
}


void read_genotype(istream&in, VBT&GT, VID&ID, VDB&twop){
  bool first_line{true};
  
  for(string id, line; in>>id; ID.push_back(id)){
    getline(in, line);
    int          ilc{0};
    stringstream ss(line);
    string       gt;
    
    for(char ig; ss>>ig; ++ilc)
      switch(ig){
      case '0':
	gt+="00";
	if(first_line) twop.push_back(0.);
	break;
      case '1': gt+="01";
	if(first_line) twop.push_back(1.);
	else ++twop[ilc];
	break;
      case '2': gt+="10";
	if(first_line) twop.push_back(2.);
	else twop[ilc]+=2.;
	break;
      default:
	throw runtime_error("Bad genotype, only 0, 1 or 2 is allowd");
      }
    
    first_line=false;
    GT.push_back(MBS(gt));
  }
  
  for(auto&x:twop) x/=ID.size();
}


void bitset_2_double(const MBS&bs, double*dd, int nlc, const VDB&twop){
  string str=bs.to_string().substr(MLC-nlc*2);
  for(auto i=0; i<nlc; ++i) dd[i] = msd[str.substr(i*2, 2)] - twop[i];
}


void worker(VBT&GT, VDB&twop, double s2pq, int i, int nlc, VDB&gi){
  double di[nlc];
  bitset_2_double(GT[i], di, nlc, twop);
  for(auto j=0; j<i; ++j){
    double dj[nlc];
    bitset_2_double(GT[j], dj, nlc, twop);
    gi.push_back(inner_product(di, di+nlc, dj, 0.)/s2pq);
  }
  gi.push_back(inner_product(di, di+nlc, di, 0.)/s2pq);
}
