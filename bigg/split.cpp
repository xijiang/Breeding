#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  const int MLC{1'000'000}
  int       NL{1'000};
  if(argc==2) NL = stoi(argv[1]);
  double    frq[MLC]{0};
  ofstream  foo, fid("id.txt");
  int       nid{0}, fn{1}, nlc{0};

  for(string line; getline(cin, line); ++nid){
    if(nid % NL == 0){
      foo.close();
      foo.open(to_string(fn++)+".gt");
    }
    stringstream ss(line);
    string id;
    foo<<line<<'\n';
    ss>>id;
    fid<<id<<'\n';
    
    nlc = 0;
    for(double td; ss>>td; frq[nlc++]+=td);
  }
  foo.close();
  foo.open("2p.txt");
  for(auto i{0}; i<nlc; ++i) foo<<frq[i]/nid<<'\n';
  return 0;
}
