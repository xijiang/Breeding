#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  ifstream fid("id.txt");
  vector<string> ID;
  for(string id; fid>>id; ID.push_back(id));
  int nid(static_cast<int>(ID.size()));

  double gi[nid];
  cout.precision(12);
  
  for(auto i=0; i<nid; ++i){
    cin.read(reinterpret_cast<char*>(gi), nid*sizeof(double));
    for(auto j=0; j<=i; ++j) cout<<ID[i]<<' '<<ID[j]<<' '<<gi[j]<<'\n';
  }
  return 0;
}
