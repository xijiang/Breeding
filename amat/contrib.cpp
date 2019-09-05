/**
 * This program calculate the genetic contribution of ID in list-a to ID in list-b.
 * Both lists are from some.ped of format {pa, ma}_nid.
 * Usage: cat some.ped | this-program list-a list-b >results
 */
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
  if(argc!=3){
    cerr<<"Usage: \n";
    cerr<<"   cat (2-col)pedigree | "<<argv[0]<<" list-1(ans) list-2(off)\n";
    return 1;
  }
  
  PED ped;			// the pedigree to look up
  if(!read_ped(cin, ped)) return 2;

  vector<int> ans, off;
  if(!read_list(argv[1], ans, ped.size()-1)) return 3;
  if(!read_list(argv[2], off, ped.size()-1)) return 4;
  sort(ans.begin(), ans.end());
  sort(off.begin(), off.end());
  if(ans.back() > off.front()){
    cerr<<"Offspring can't contribute to ancestors\n";
    return 5;
  }

  map<PM, double> mid;		// store the mid results of Amat
  
  cout.precision(12);
  for(const auto&i:ans)
    for(const auto&j:off)
      cout<<i<<'\t'<<j<<'\t'<<Amat(i, j, ped, mid)<<'\n';

  return 0;
}
