// Summarize the per-neutron classification produced by MyEnergyAnalysis.
// Usage:
//   root -l -b -q 'PlotNeutronCategories.C("neutron_category_test.root")'

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TKey.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
  h.SetStats(false);
  h.SetLineWidth(2);
  h.GetYaxis()->SetTitle("Number of neutrons");
  h.GetYaxis()->SetTitleOffset(1.25);
  h.GetXaxis()->SetLabelSize(0.04);
}

std::string processCombinationLabel(int mask) {
  if(mask==0) return "No recorded interaction";
  std::vector<std::string> names;
  if(mask&1) names.emplace_back("Elastic");
  if(mask&2) names.emplace_back("Inelastic");
  if(mask&4) names.emplace_back("Capture");
  if(mask&8) names.emplace_back("Other");
  std::ostringstream out;
  for(std::size_t i=0;i<names.size();++i) {
    if(i) out<<" + ";
    out<<names[i];
  }
  return out.str();
}

void csvRow(std::ofstream& out,
            const std::string& section,
            const std::string& parent,
            const std::string& category,
            long long count,
            long long total,
            long long parentCount,
            const std::string& note) {
  out<<'"'<<section<<"\",\""<<parent<<"\",\""<<category<<"\","<<count<<','
     <<std::setprecision(8)<<(total?double(count)/total:0.0)<<','
     <<(parentCount?double(count)/parentCount:0.0)<<",\""<<note<<"\"\n";
}

struct NeutronRecord {
  double birthKE=0.0;
  int outcome=-1;
};
}

void PlotNeutronCategories(const char* filename="neutron_category_test.root") {
  TFile file(filename,"READ");
  if(file.IsZombie()){std::cerr<<"ERROR: cannot open "<<filename<<'\n';return;}
  TTree* tree=findNeutronTree(&file);
  if(!tree){std::cerr<<"ERROR: NeutronTree not found\n";return;}

  const char* requiredBranches[]={"BirthKE","Invisible","HadronicInteraction","Escaped","Outcome",
                                  "NInteractions","NElastic","NInelastic","NCapture","NOther"};
  for(const char* branch:requiredBranches) {
    if(!tree->GetBranch(branch)) {
      std::cerr<<"ERROR: NeutronTree is missing required branch "<<branch
               <<". This ROOT file predates the neutron-category output.\n";
      return;
    }
  }

  bool invisible=false,hadronic=false,escaped=false;
  int outcome=-1,nInteractions=0,nElastic=0,nInelastic=0,nCapture=0,nOther=0;
  double birthKE=0.0;
  tree->SetBranchAddress("BirthKE",&birthKE);
  tree->SetBranchAddress("Invisible",&invisible);
  tree->SetBranchAddress("HadronicInteraction",&hadronic);
  tree->SetBranchAddress("Escaped",&escaped);
  tree->SetBranchAddress("Outcome",&outcome);
  tree->SetBranchAddress("NInteractions",&nInteractions);
  tree->SetBranchAddress("NElastic",&nElastic);
  tree->SetBranchAddress("NInelastic",&nInelastic);
  tree->SetBranchAddress("NCapture",&nCapture);
  tree->SetBranchAddress("NOther",&nOther);

  std::array<long long,4> outcomes{{0,0,0,0}};
  std::array<long long,16> processCombinations{};
  long long nInvisible=0,nAnyInteraction=0,nHadronic=0;
  long long hasElastic=0,hasInelastic=0,hasCapture=0,hasOther=0;
  std::vector<NeutronRecord> records;
  records.reserve(tree->GetEntries());
  double maxBirthKE=0.0;

  for(Long64_t i=0;i<tree->GetEntries();++i) {
    tree->GetEntry(i);
    if(outcome>=0&&outcome<4) ++outcomes[outcome];
    if(invisible) ++nInvisible;
    if(nInteractions>0) ++nAnyInteraction;
    if(hadronic) ++nHadronic;
    if(nElastic>0) ++hasElastic;
    if(nInelastic>0) ++hasInelastic;
    if(nCapture>0) ++hasCapture;
    if(nOther>0) ++hasOther;

    int processMask=0;
    if(nElastic>0) processMask|=1;
    if(nInelastic>0) processMask|=2;
    if(nCapture>0) processMask|=4;
    if(nOther>0) processMask|=8;
    ++processCombinations[processMask];

    records.push_back({birthKE,outcome});
    if(std::isfinite(birthKE)) maxBirthKE=std::max(maxBirthKE,birthKE);
  }

  const long long total=tree->GetEntries();
  const long long nVisible=total-nInvisible;
  const long long nNoInteraction=total-nAnyInteraction;
  const long long nLeaves=outcomes[0]+outcomes[1];
  const long long nStays=outcomes[2]+outcomes[3];
  const long long outcomeCheck=outcomes[0]+outcomes[1]+outcomes[2]+outcomes[3];
  long long processCheck=0;
  for(long long count:processCombinations) processCheck+=count;
  if(outcomeCheck!=total || processCheck!=total)
    std::cerr<<"WARNING: category totals do not reconcile with NeutronTree entries.\n";

  std::cout<<"\nNeutron category hierarchy (neutrons entering the active volume)\n"
           <<"  all neutrons: "<<total<<'\n'
           <<"    leaves active volume: "<<nLeaves<<'\n'
           <<"      leaves without interaction: "<<outcomes[0]<<'\n'
           <<"      interacts then leaves: "<<outcomes[1]<<'\n'
           <<"    stays in active volume: "<<nStays<<'\n'
           <<"      interacts and stays: "<<outcomes[2]<<'\n'
           <<"      no interaction and stays: "<<outcomes[3]<<'\n'
           <<"    visible: "<<nVisible<<"; invisible: "<<nInvisible<<'\n'
           <<"    any interaction: "<<nAnyInteraction
           <<"; no interaction: "<<nNoInteraction<<'\n'
           <<"    any hadronic interaction: "<<nHadronic<<"\n\n"
           <<"  Non-overlapping process combinations among interacting neutrons:\n";
  for(int mask=1;mask<16;++mask)
    if(processCombinations[mask]>0)
      std::cout<<"    "<<processCombinationLabel(mask)<<": "<<processCombinations[mask]<<'\n';

  std::ofstream summary("neutron_category_summary.csv");
  summary<<"metric,count,fraction\n";
  auto row=[&](const char* name,long long count) {
    summary<<name<<','<<count<<','<<std::setprecision(8)
           <<(total?double(count)/total:0.0)<<'\n';
  };
  row("total",total);
  row("leaves_without_interaction",outcomes[0]);
  row("interacts_then_leaves",outcomes[1]);
  row("interacts_and_stays",outcomes[2]);
  row("no_interaction_and_stays",outcomes[3]);
  row("visible",nVisible);
  row("invisible",nInvisible);
  row("no_interaction",nNoInteraction);
  row("any_interaction",nAnyInteraction);
  row("any_hadronic_interaction",nHadronic);
  row("hadElastic",hasElastic);
  row("neutronInelastic",hasInelastic);
  row("nCapture",hasCapture);
  row("other",hasOther);

  // Each of the first four sections is exclusive within itself. Process flags overlap.
  std::ofstream relationships("neutron_category_relationships.csv");
  relationships<<"section,parent,category,count,fraction_of_all,fraction_of_parent,note\n";
  csvRow(relationships,"Outcome","All neutrons","Leaves without interaction",outcomes[0],total,total,"Exclusive outcome");
  csvRow(relationships,"Outcome","All neutrons","Interacts then leaves",outcomes[1],total,total,"Exclusive outcome");
  csvRow(relationships,"Outcome","All neutrons","Interacts and stays",outcomes[2],total,total,"Exclusive outcome");
  csvRow(relationships,"Outcome","All neutrons","No interaction and stays",outcomes[3],total,total,"Exclusive outcome");
  csvRow(relationships,"Visibility","All neutrons","Visible",nVisible,total,total,"FamilyEdep is at or above the configured threshold");
  csvRow(relationships,"Visibility","All neutrons","Invisible",nInvisible,total,total,"FamilyEdep is below the configured threshold");
  csvRow(relationships,"Interaction status","All neutrons","No recorded interaction",nNoInteraction,total,total,"Exclusive interaction status");
  csvRow(relationships,"Interaction status","All neutrons","Any recorded interaction",nAnyInteraction,total,total,"Exclusive interaction status");
  for(int mask=1;mask<16;++mask)
    if(processCombinations[mask]>0)
      csvRow(relationships,"Process combination","Any recorded interaction",processCombinationLabel(mask),
             processCombinations[mask],total,nAnyInteraction,"Mutually exclusive process combination");
  csvRow(relationships,"Process flag","All neutrons","Has elastic scattering",hasElastic,total,total,"Overlaps other process flags");
  csvRow(relationships,"Process flag","All neutrons","Has inelastic interaction",hasInelastic,total,total,"Overlaps other process flags");
  csvRow(relationships,"Process flag","All neutrons","Has neutron capture",hasCapture,total,total,"Overlaps other process flags");
  csvRow(relationships,"Process flag","All neutrons","Has other interaction",hasOther,total,total,"Overlaps other process flags");

  gStyle->SetOptStat(0);
  TH1D hOutcome("hOutcome","Exclusive neutron outcomes",4,0,4);
  style(hOutcome);
  hOutcome.SetFillColor(kAzure+1);
  hOutcome.SetLineColor(kAzure+3);
  const char* outcomeLabels[4]={"Leaves directly","Interacts + leaves","Interacts + stays","No interaction + stays"};
  for(int i=0;i<4;++i) {
    hOutcome.SetBinContent(i+1,outcomes[i]);
    hOutcome.GetXaxis()->SetBinLabel(i+1,outcomeLabels[i]);
  }
  TCanvas c1("c1","Exclusive neutron outcomes",1100,760);
  c1.SetLeftMargin(0.12);
  c1.SetBottomMargin(0.17);
  hOutcome.Draw("HIST TEXT0");
  c1.SaveAs("neutron_outcomes_exclusive.png");

  std::vector<int> observedMasks;
  for(int mask=1;mask<16;++mask)
    if(processCombinations[mask]>0) observedMasks.push_back(mask);
  const int nCombinationBins=std::max(1,static_cast<int>(observedMasks.size()));
  TH1D hCombinations("hCombinations","Non-overlapping interaction-process combinations",
                     nCombinationBins,0,nCombinationBins);
  style(hCombinations);
  hCombinations.SetFillColor(kGreen+1);
  hCombinations.SetLineColor(kGreen+3);
  for(std::size_t i=0;i<observedMasks.size();++i) {
    hCombinations.SetBinContent(i+1,processCombinations[observedMasks[i]]);
    hCombinations.GetXaxis()->SetBinLabel(i+1,processCombinationLabel(observedMasks[i]).c_str());
  }
  if(observedMasks.empty()) hCombinations.GetXaxis()->SetBinLabel(1,"No interacting neutrons");
  TCanvas c2("c2","Interaction-process combinations",1250,780);
  c2.SetLeftMargin(0.12);
  c2.SetBottomMargin(0.25);
  hCombinations.GetXaxis()->LabelsOption("v");
  hCombinations.Draw("HIST TEXT0");
  c2.SaveAs("neutron_process_combinations.png");

  // The stack uses mutually exclusive outcomes, so it sums to the full sample.
  const double xMax=(maxBirthKE>0.0)?std::max(0.01,std::ceil(maxBirthKE*10.0)/10.0):1.0;
  const int nBins=60;
  const int colors[4]={kAzure+1,kOrange+1,kGreen+1,kViolet+1};
  TH1D hBirth0("hBirth0","",nBins,0,xMax);
  TH1D hBirth1("hBirth1","",nBins,0,xMax);
  TH1D hBirth2("hBirth2","",nBins,0,xMax);
  TH1D hBirth3("hBirth3","",nBins,0,xMax);
  TH1D* birthHists[4]={&hBirth0,&hBirth1,&hBirth2,&hBirth3};
  for(int i=0;i<4;++i) {
    birthHists[i]->SetFillColor(colors[i]);
    birthHists[i]->SetLineColor(colors[i]);
  }
  for(const auto& record:records)
    if(record.outcome>=0&&record.outcome<4&&std::isfinite(record.birthKE))
      birthHists[record.outcome]->Fill(record.birthKE);

  THStack energyStack("energyStack","Initial neutron kinetic energy by exclusive outcome;Initial neutron kinetic energy [GeV];Number of neutrons");
  for(auto* hist:birthHists) energyStack.Add(hist);
  TCanvas c3("c3","Initial neutron kinetic energy",1150,780);
  c3.SetLeftMargin(0.12);
  c3.SetRightMargin(0.04);
  energyStack.Draw("HIST");
  energyStack.GetYaxis()->SetTitleOffset(1.25);
  TLegend legend(0.58,0.68,0.94,0.91);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  for(int i=0;i<4;++i) legend.AddEntry(birthHists[i],outcomeLabels[i],"f");
  legend.Draw();
  c3.SaveAs("neutron_initial_energy_by_outcome.png");

  std::cout<<"\nCreated:\n"
           <<"  neutron_category_summary.csv\n"
           <<"  neutron_category_relationships.csv\n"
           <<"  neutron_outcomes_exclusive.png\n"
           <<"  neutron_process_combinations.png\n"
           <<"  neutron_initial_energy_by_outcome.png\n";
}

