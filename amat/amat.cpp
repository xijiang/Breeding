#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>

using namespace std;
using PM =tuple<int, int>;
using PED=vector<PM>;
using MID=map<PM, double>;	// store intermediate resultes

bool nARG(int argc, string cmd){
  if(argc!=3){
    cerr<<"Usage: cat a-pedigree | " << cmd << " option ID-list\n";
    cerr<<"\n where options:\n";
    cerr<<"     A/a: calculate A matrix, output binary format\n";
    cerr<<"     F/f: inbreeding values of the ID\n";
    cerr<<"     T/t: diagonal D inverse and T inverse for calculation A inverse\n";
    //cerr<<"     V/v: inverse of the A matrix\n";
    return false;
  }
  return true;
}


bool read_ped(istream&in, PED&ped){
  clog<<"LOG: Reading the pedigree ...\n";
  ped.push_back({0,0});		// the magic dummy
				// hence, row number is ID number id starts from 1.
  for(int pa, ma; in>>pa>>ma; ped.push_back({pa, ma}))
    if(pa>=static_cast<int>(ped.size())
       || ma>=static_cast<int>(ped.size())
       || pa<0 || ma<0){
      cerr<<"ERROR: Invalid pa / ma ID\n";
      return false;
    }
  return true;
}


bool read_list(string file, vector<int>&ilist, size_t nid){
  clog<<"LOG: Reading the ID list ...\n";
  
  int oid{0};
  ifstream fin(file);
  if(!fin.is_open()){
    cerr << "File not exists\n";
    return false;
  }
  for(int id; fin>>id; ilist.push_back(id)){
    if(id>static_cast<int>(nid) || id<0){
      cerr<<"ERROR: ID not in the pedigree\n";
      return false;		// id must be in 1:nid and sorted
    }
    if(id<=oid) clog<<"Warning: ID list not sorted\n";
    oid = id;
  }
  return true;
}


double Amat(int i, int j, const PED &ped, MID&mid){
  // clog<<i<<' '<<j<<'\n'; //uncomment this to check the stack
  if(i==0 || j==0) return 0;	// so that relationship with un unknown is not stored in mid
  
  if(i>j) swap(i, j);
  
  // Look up if {i,j} was calculated or not
  if(mid.find({i, j})!=mid.end()) return mid[{i, j}];

  const auto &[pa, ma] = ped[j];
  if(i==j)
    mid[{j, j}] = 1 + Amat(pa, ma, ped, mid) / 2.;
  else
    mid[{i, j}] = (Amat(i, pa, ped, mid) + Amat(i, ma, ped, mid)) / 2.;

  return mid[{i, j}];
}

void DnT(const PED&ped, MID&mid, vector<int>&ilist){
  ofstream dfile{"D.vec"}, tfile{"T.mat"};
  // Inbred coefficients of parents of ilist
  map<int, double> pma;
  pma[0] = -1;
  for(const auto&id:ilist){
    const auto&[pa, ma] = ped[id];
    if(pma.find(pa)==pma.end()) pma[pa] = Amat(pa, pa, ped, mid) -1;
    if(pma.find(ma)==pma.end()) pma[ma] = Amat(ma, ma, ped, mid) -1;
  }

  dfile.precision(12);

  // directly write D and Ti to file
  for(const auto&id:ilist){
    const auto&[pa, ma] = ped[id];

    // D inverse
    dfile << 1./(.5 - .25*(pma[pa] + pma[ma])) << '\n';

    // T inverse
    auto putT = [](int a, int b, int c, ostream&tfile){
		  if(a) tfile << c << ' ' << a << ' ' << -.5 << '\n';
		  if(b) tfile << c << ' ' << b << ' ' << -.5 << '\n';
		  tfile       << c << ' ' << c << ' ' <<   1 << '\n';
		};
    if(pa<ma) putT(pa, ma, id, tfile);
    else      putT(ma, pa, id, tfile);
  }
}

int main(int argc, char *argv[])
{
  if(!nARG(argc, argv[0])) return 1;
  
  PED ped;			// the pedigree to look up
  if(!read_ped(cin, ped)) return 2;

  vector<int> ilist;
  map<PM, double> mid;		// store the mid results of Amat
  
  switch(argv[1][0]){
  case 'A':
  case 'a':
    clog<<"Calculate A matrix of listed ID\n";
    if(!read_list(argv[2], ilist, ped.size()-1)) return 3;
    for(auto&i:ilist)
      for(auto&j:ilist){
	if(j>i) break;
	cout<<i<<'\t'<<j<<'\t'<<Amat(i, j, ped, mid)<<'\n';
      }
    break;
    
  case 'F':
  case 'f':
    clog<<"Calculate inbreeding values of listed ID\n";
    if(!read_list(argv[2], ilist, ped.size()-1)) return 3;
    for(auto&id:ilist) cout<<id<<'\t'<<Amat(id, id, ped, mid) - 1<<'\n';
    break;

  case 'T':
  case 't':
    clog<<"Calculate inverse D and T matrix for inverse A construction\n";
    clog<<"Results are in D.vec and Ti.mat.\n";
    if(!read_list(argv[2], ilist, ped.size()-1)) return 3;
    DnT(ped, mid, ilist);
    break;
    
  default:
    cerr<<"ERROR: Invalid option \'"<<argv[1]<<"\'\n";
    return 4;
  }

  clog<<"LOG: Number of items of intermediate results: "<<mid.size()<<'\n';
  return 0;
}
