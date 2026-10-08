// Diagnose neutral escape accounting using existing, unmerged grid outputs.
// Recommended:
// root -l -b -q 'DiagnoseNeutralEscape.C("clean_input_files.txt",-1)'
// Matching is reset for each file because EscapeTree lacks Run/SubRun.

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {
TTree* findTree(TDirectory* d, const char* name) {
  if (!d) return nullptr;
  if (auto* t=dynamic_cast<TTree*>(d->Get(name))) return t;
  TIter next(d->GetListOfKeys());
  while (auto* k=dynamic_cast<TKey*>(next())) {
    if (std::string(k->GetClassName()).find("TDirectory")==std::string::npos) continue;
    if (auto* t=findTree(dynamic_cast<TDirectory*>(k->ReadObj()),name)) return t;
  }
  return nullptr;
}
std::string quote(const std::string& s) {
  std::string q="\""; for(char c:s){if(c=='"')q+='"';q+=c;} return q+'"';
}
bool endsWith(const std::string& s,const std::string& x) {
  return s.size()>=x.size() && s.compare(s.size()-x.size(),x.size(),x)==0;
}
struct Evt { Long64_t entry; int run,subrun,event; double ledger; };
struct Esc { Long64_t entry; int track,pdg; double birth,exit; bool charged; };
struct Counts { long long files=0,failed=0,events=0,withN=0,negative=0,ambiguous=0,mismatch=0; };

void processFile(const std::string& fn,long long limit,Counts& n,
                 std::ofstream& out,std::ofstream& details) {
  TFile f(fn.c_str(),"READ");
  if(f.IsZombie()){std::cerr<<"ERROR opening "<<fn<<'\n';++n.failed;return;}
  TTree* et=findTree(&f,"EventTree"); TTree* xt=findTree(&f,"EscapeTree");
  if(!et||!xt){std::cerr<<"ERROR missing tree in "<<fn<<'\n';++n.failed;return;}

  int ev=0,run=0,sub=0; double ledger=0;
  et->SetBranchAddress("Event",&ev); et->SetBranchAddress("Run",&run);
  et->SetBranchAddress("SubRun",&sub); et->SetBranchAddress("E_escape_neutral",&ledger);
  std::vector<Evt> events; std::map<int,int> multiplicity;
  for(Long64_t i=0;i<et->GetEntries();++i){
    if(limit>=0 && n.events+(long long)events.size()>=limit) break;
    et->GetEntry(i); events.push_back({i,run,sub,ev,ledger}); ++multiplicity[ev];
  }

  int xev=0,tr=0,pdg=0; double birth=0,exit=0; bool charged=false;
  xt->SetBranchAddress("Event",&xev); xt->SetBranchAddress("TrackID",&tr);
  xt->SetBranchAddress("PDG",&pdg); xt->SetBranchAddress("BirthKE",&birth);
  xt->SetBranchAddress("ExitKE",&exit); xt->SetBranchAddress("Charged",&charged);
  std::map<int,std::vector<Esc>> byEvent;
  for(Long64_t i=0;i<xt->GetEntries();++i){xt->GetEntry(i);byEvent[xev].push_back({i,tr,pdg,birth,exit,charged});}

  for(const auto& e:events){
    bool amb=multiplicity[e.event]>1; if(amb)++n.ambiguous;
    double maxN=-std::numeric_limits<double>::infinity(),sum=0; int maxTrack=-1,nNeutral=0,nNeutron=0;
    auto it=byEvent.find(e.event);
    if(it!=byEvent.end()) for(const auto& x:it->second){
      if(x.charged)continue; ++nNeutral; sum+=x.exit;
      if(x.pdg==2112){++nNeutron;if(x.exit>maxN){maxN=x.exit;maxTrack=x.track;}}
    }
    bool hasN=maxTrack>=0; double diff=hasN?e.ledger-maxN:std::numeric_limits<double>::quiet_NaN();
    double ledgerDiff=e.ledger-sum; bool neg=hasN&&diff < -1e-9; bool bad=std::fabs(ledgerDiff)>1e-8;
    ++n.events; if(hasN)++n.withN; if(neg)++n.negative; if(bad)++n.mismatch;
    out<<quote(fn)<<','<<e.entry<<','<<e.run<<','<<e.subrun<<','<<e.event<<','<<multiplicity[e.event]
       <<','<<(amb?1:0)<<','<<std::setprecision(17)<<e.ledger<<',';
    if(hasN)out<<maxN<<','<<maxTrack<<','<<diff;else out<<"nan,-1,nan";
    out<<','<<sum<<','<<ledgerDiff<<','<<nNeutral<<','<<nNeutron<<','<<(neg?1:0)<<','<<(bad?1:0)<<'\n';
    if(neg&&it!=byEvent.end())for(const auto& x:it->second)if(!x.charged)
      details<<quote(fn)<<','<<e.entry<<','<<e.run<<','<<e.subrun<<','<<e.event<<','<<x.entry
             <<','<<x.track<<','<<x.pdg<<','<<std::setprecision(17)<<x.birth<<','<<x.exit<<'\n';
  }
  ++n.files;
}
}

void DiagnoseNeutralEscape(const char* inputName,int maxTotalEvents=-1) {
  std::string input=inputName?inputName:""; std::vector<std::string> files;
  if(endsWith(input,".txt")){
    std::ifstream list(input); if(!list){std::cerr<<"ERROR opening list "<<input<<'\n';return;}
    std::string s; while(std::getline(list,s)){
      auto a=s.find_first_not_of(" \t\r\n"); if(a==std::string::npos||s[a]=='#')continue;
      auto b=s.find_last_not_of(" \t\r\n"); files.push_back(s.substr(a,b-a+1));
    }
  } else files.push_back(input);

  std::ofstream out("neutral_escape_diagnostic_all_files.csv");
  out<<"source_file,event_tree_entry,run,subrun,event,event_id_multiplicity_within_file,ambiguous_within_file,"
       "E_escape_neutral,max_neutron_ExitKE_field,max_neutron_track,E_escape_neutral_minus_max_neutron,"
       "sum_neutral_ExitKE_fields,E_escape_neutral_minus_sum_fields,n_neutral_escape_rows,n_neutron_escape_rows,is_negative,ledger_mismatch\n";
  std::ofstream details("neutral_escape_negative_contributors.csv");
  details<<"source_file,event_tree_entry,run,subrun,event,escape_tree_entry,track,pdg,birthKE,ExitKE_field\n";
  Counts n; std::cout<<"Processing "<<files.size()<<" unmerged ROOT files.\n";
  for(size_t i=0;i<files.size();++i){
    if(maxTotalEvents>=0&&n.events>=maxTotalEvents)break;
    processFile(files[i],maxTotalEvents,n,out,details);
    if((i+1)%25==0||i+1==files.size())std::cout<<"  attempted "<<i+1<<'/'<<files.size()<<", events "<<n.events<<'\n';
  }
  std::cout<<"\nSummary\n  files read: "<<n.files<<"\n  files failed: "<<n.failed
           <<"\n  events: "<<n.events<<"\n  events with escaping neutron: "<<n.withN
           <<"\n  ambiguous IDs within a file: "<<n.ambiguous
           <<"\n  negative E_escape_neutral - max neutron ExitKE: "<<n.negative
           <<"\n  ledger mismatches: "<<n.mismatch
           <<"\nCreated neutral_escape_diagnostic_all_files.csv and neutral_escape_negative_contributors.csv\n";
}
