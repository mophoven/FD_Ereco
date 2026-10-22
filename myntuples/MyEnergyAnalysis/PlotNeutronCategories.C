// Summarize the per-neutron classification produced by MyEnergyAnalysis.
// Usage:
//   root -l -b -q 'PlotNeutronCategories.C("neutron_category_test.root")'

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TStyle.h>
#include <TTree.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
TTree* findNeutronTree(TDirectory* dir) {
  if(!dir) return nullptr;
  if(auto* tree=dynamic_cast<TTree*>(dir->Get("NeutronTree"))) return tree;
  TIter next(dir->GetListOfKeys());
  while(auto* key=dynamic_cast<TKey*>(next())) {
    if(std::string(key->GetClassName()).find("TDirectory")==std::string::npos) continue;
    if(auto* tree=findNeutronTree(dynamic_cast<TDirectory*>(key->ReadObj()))) return tree;
  }
  return nullptr;
}
void style(TH1D& h) {
  h.SetStats(false); h.SetFillColor(kAzure+1); h.SetLineColor(kAzure+3);
  h.SetLineWidth(2); h.GetYaxis()->SetTitle("Number of neutrons");
  h.GetYaxis()->SetTitleOffset(1.25); h.GetXaxis()->SetLabelSize(0.045);
}
}

void PlotNeutronCategories(const char* filename="neutron_category_test.root") {
  TFile file(filename,"READ");
  if(file.IsZombie()){std::cerr<<"ERROR: cannot open "<<filename<<'\n';return;}
  TTree* tree=findNeutronTree(&file);
  if(!tree){std::cerr<<"ERROR: NeutronTree not found\n";return;}

  bool invisible=false,hadronic=false,escaped=false;
  int outcome=-1,nInteractions=0,nElastic=0,nInelastic=0,nCapture=0,nOther=0;
  tree->SetBranchAddress("Invisible",&invisible);
  tree->SetBranchAddress("HadronicInteraction",&hadronic);
  tree->SetBranchAddress("Escaped",&escaped);
  tree->SetBranchAddress("Outcome",&outcome);
  tree->SetBranchAddress("NInteractions",&nInteractions);
  tree->SetBranchAddress("NElastic",&nElastic);
  tree->SetBranchAddress("NInelastic",&nInelastic);
  tree->SetBranchAddress("NCapture",&nCapture);
  tree->SetBranchAddress("NOther",&nOther);

  long long outcomes[4]={0,0,0,0};
  long long nInvisible=0,nAnyInteraction=0,nHadronic=0;
  long long hasElastic=0,hasInelastic=0,hasCapture=0,hasOther=0;
  for(Long64_t i=0;i<tree->GetEntries();++i){
    tree->GetEntry(i);
    if(outcome>=0&&outcome<4)++outcomes[outcome];
    if(invisible)++nInvisible;
    if(nInteractions>0)++nAnyInteraction;
    if(hadronic)++nHadronic;
    if(nElastic>0)++hasElastic;
    if(nInelastic>0)++hasInelastic;
    if(nCapture>0)++hasCapture;
    if(nOther>0)++hasOther;
  }
  const long long total=tree->GetEntries();
  std::cout<<"\nNeutron summary (neutrons that entered the active volume)\n"
           <<"  total: "<<total<<'\n'
           <<"  leaves without interaction: "<<outcomes[0]<<'\n'
           <<"  interacts then leaves: "<<outcomes[1]<<'\n'
           <<"  interacts and stays: "<<outcomes[2]<<'\n'
           <<"  no recorded interaction and stays: "<<outcomes[3]<<'\n'
           <<"  invisible: "<<nInvisible<<'\n'
           <<"  any interaction: "<<nAnyInteraction<<'\n'
           <<"  any hadronic interaction: "<<nHadronic<<'\n';

  std::ofstream csv("neutron_category_summary.csv");
  csv<<"metric,count,fraction\n";
  auto row=[&](const char* name,long long count){csv<<name<<','<<count<<','<<std::setprecision(8)<<(total?double(count)/total:0.0)<<'\n';};
  row("total",total); row("leaves_without_interaction",outcomes[0]);
  row("interacts_then_leaves",outcomes[1]); row("interacts_and_stays",outcomes[2]);
  row("no_interaction_and_stays",outcomes[3]); row("invisible",nInvisible);
  row("any_interaction",nAnyInteraction); row("any_hadronic_interaction",nHadronic);
  row("hadElastic",hasElastic); row("neutronInelastic",hasInelastic);
  row("nCapture",hasCapture); row("other",hasOther);

  gStyle->SetOptStat(0);
  TH1D hOutcome("hOutcome","Exclusive neutron outcomes",4,0,4); style(hOutcome);
  const char* outcomeLabels[4]={"Leaves directly","Interacts + leaves","Interacts + stays","No interaction + stays"};
  for(int i=0;i<4;++i){hOutcome.SetBinContent(i+1,outcomes[i]);hOutcome.GetXaxis()->SetBinLabel(i+1,outcomeLabels[i]);}
  TCanvas c1("c1","Exclusive neutron outcomes",1100,760); c1.SetLeftMargin(0.12); c1.SetBottomMargin(0.17);
  hOutcome.Draw("HIST TEXT0"); c1.SaveAs("neutron_outcomes_exclusive.png");

  TH1D hRequested("hRequested","Requested neutron counts (categories can overlap)",4,0,4); style(hRequested);
  const char* labels[4]={"Leaves directly","Invisible","Any interaction","Hadronic interaction"};
  const long long values[4]={outcomes[0],nInvisible,nAnyInteraction,nHadronic};
  for(int i=0;i<4;++i){hRequested.SetBinContent(i+1,values[i]);hRequested.GetXaxis()->SetBinLabel(i+1,labels[i]);}
  TCanvas c2("c2","Requested neutron counts",1100,760); c2.SetLeftMargin(0.12); c2.SetBottomMargin(0.17);
  hRequested.Draw("HIST TEXT0"); c2.SaveAs("neutron_requested_counts.png");

  std::cout<<"Created neutron_category_summary.csv, neutron_outcomes_exclusive.png, and neutron_requested_counts.png\n";
}
