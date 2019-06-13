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
    
  for(int id; fin>>id; ilist.push_back(id)){
    if(id>static_cast<int>(nid) || id<0 || id<=oid){
      cerr<<"ERROR: ID not in the pedigree\n";
      return false;		// id must be in 1:nid and sorted
    }
    if(id<=oid){
      cerr<<"ERROR: ID list not sorted\n";
      return false;
    }
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


int main(int argc, char *argv[])
{
  if(!nARG(argc, argv[0])) return 1;
  
  PED ped;			// the pedigree to look up
  if(!read_ped(cin, ped)) return 2;

  vector<int> ilist;
  if(!read_list(argv[2], ilist, ped.size()-1)) return 3;
  map<PM, double> mid;		// store the mid results of Amat
  
  switch(argv[1][0]){
  case 'A':
  case 'a':
    // calculate A of listed ID
    for(auto&i:ilist)
      for(auto&j:ilist){
	if(j>i) break;
	cout<<i<<'\t'<<j<<'\t'<<Amat(i, j, ped, mid)<<'\n';
      }
    break;
    
  case 'F':
  case 'f':
    // calculate inbreeding values of listed ID
    for(auto&id:ilist) cout<<id<<'\t'<<Amat(id, id, ped, mid) - 1<<'\n';
    break;

  default:
    cerr<<"ERROR: Invalid option \'"<<argv[1]<<"\'\n";
    return 4;
  }

  clog<<"LOG: Number of items of intermediate results: "<<mid.size()<<'\n';
  return 0;
}
