// Meeting-ready neutron plots. Run with:
// root -l -b -q 'Neutron20EventAnalysis.C("MiloEnergy_merged_500.root",20)'
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
TTree *findTree(TDirectory *d, const char *name) {
  if (!d) return nullptr;
  if (auto *t = dynamic_cast<TTree *>(d->Get(name))) return t;
  TIter next(d->GetListOfKeys());
  while (auto *k = dynamic_cast<TKey *>(next())) {
    if (std::string(k->GetClassName()).find("TDirectory") == std::string::npos) continue;
    if (auto *t = findTree(dynamic_cast<TDirectory *>(k->ReadObj()), name)) return t;
  }
  return nullptr;
}
struct Step { double time, ke; };
struct Neutron { int event, track; double birth, exit; };
std::vector<double> history(const Neutron &n, const std::map<int,std::map<int,std::vector<Step>>> &steps) {
  std::vector<double> h{n.birth};
  auto e=steps.find(n.event);
  if(e!=steps.end()) { auto t=e->second.find(n.track); if(t!=e->second.end()) for(const auto&s:t->second) h.push_back(s.ke); }
  h.push_back(n.exit); return h;
}
TGraph *makeGraph(const std::vector<double>&h, Color_t c) {
  auto*g=new TGraph(h.size()); for(size_t i=0;i<h.size();++i) g->SetPoint(i,i,h[i]);
  g->SetLineColor(c); g->SetMarkerColor(c); g->SetLineWidth(3); g->SetMarkerStyle(20); g->SetMarkerSize(1.1); return g;
}
}

void Neutron20EventAnalysis(const char *fileName, int maxEvents=20) {
  TFile f(fileName,"READ"); if(f.IsZombie()){std::cerr<<"Cannot open "<<fileName<<'\n';return;}
  auto*et=findTree(&f,"EventTree"),*vt=findTree(&f,"VertexTree"),*xt=findTree(&f,"EscapeTree");
  if(!et||!vt||!xt){std::cerr<<"Need EventTree, VertexTree and EscapeTree.\n";return;}
  int ev=0; et->SetBranchAddress("Event",&ev); std::vector<int>ids; std::set<int>chosen;
  for(Long64_t i=0;i<et->GetEntries()&&(int)ids.size()<maxEvents;++i){et->GetEntry(i);if(chosen.insert(ev).second)ids.push_back(ev);} et->ResetBranchAddresses();

  int ve=0,pdg=0,trk=0; double ke=0,time=0;
  vt->SetBranchAddress("Event",&ve);vt->SetBranchAddress("In_PDG",&pdg);vt->SetBranchAddress("In_TrackID",&trk);vt->SetBranchAddress("In_KE",&ke);vt->SetBranchAddress("Vtx_t",&time);
  std::map<int,std::map<int,std::vector<Step>>>steps;
  for(Long64_t i=0;i<vt->GetEntries();++i){vt->GetEntry(i);if(chosen.count(ve)&&pdg==2112)steps[ve][trk].push_back({time,ke});}
  for(auto&e:steps)for(auto&t:e.second)std::sort(t.second.begin(),t.second.end(),[](const Step&a,const Step&b){return a.time<b.time;});

  int xe=0,xtrk=0,xpdg=0;double birth=0,exit=0;
  xt->SetBranchAddress("Event",&xe);xt->SetBranchAddress("TrackID",&xtrk);xt->SetBranchAddress("PDG",&xpdg);xt->SetBranchAddress("BirthKE",&birth);xt->SetBranchAddress("ExitKE",&exit);
  std::map<int,Neutron>bestByEvent;
  for(Long64_t i=0;i<xt->GetEntries();++i){xt->GetEntry(i);if(!chosen.count(xe)||xpdg!=2112)continue;auto old=bestByEvent.find(xe);if(old==bestByEvent.end()||birth>old->second.birth)bestByEvent[xe]={xe,xtrk,birth,exit};}
  std::vector<Neutron>ns;for(int id:ids)if(bestByEvent.count(id))ns.push_back(bestByEvent[id]);
  if(ns.empty()){std::cout<<"No escaping neutrons in selected events.\n";return;}
  gStyle->SetOptStat(0);gStyle->SetTitleSize(.045,"XY");gStyle->SetLabelSize(.038,"XY");
  gStyle->SetTitleOffset(1.35,"Y");gStyle->SetTitleOffset(1.15,"X");

  int n=ns.size();double ymax=0;
  auto*hb=new TH1D("hb","Highest-energy escaping neutron per event;Event;Kinetic energy [GeV]",n,.5,n+.5);
  auto*he=new TH1D("he","",n,.5,n+.5);
  for(int i=0;i<n;++i){hb->SetBinContent(i+1,ns[i].birth);he->SetBinContent(i+1,ns[i].exit);hb->GetXaxis()->SetBinLabel(i+1,Form("%d",ns[i].event));ymax=std::max(ymax,ns[i].birth);}
  hb->SetMinimum(0);hb->SetMaximum(1.18*ymax);hb->SetMarkerStyle(20);hb->SetMarkerSize(1.3);hb->SetMarkerColor(kBlue+1);he->SetMarkerStyle(21);he->SetMarkerSize(1.3);he->SetMarkerColor(kOrange+7);
  auto*c1=new TCanvas("cSummary","Summary",1250,750);c1->SetGridy();c1->SetLeftMargin(.16);c1->SetRightMargin(.04);c1->SetBottomMargin(.15);hb->Draw("P");he->Draw("P SAME");
  auto*leg=new TLegend(.68,.76,.89,.89);leg->SetBorderSize(0);leg->AddEntry(hb,"Birth KE","p");leg->AddEntry(he,"Exit KE","p");leg->Draw();c1->SaveAs("meeting_highest_escaping_neutron_per_event.png");

  auto*hl=new TH1D("hl","Energy lost by highest-energy escaping neutron;Event;Birth KE - exit KE [MeV]",n,.5,n+.5);
  for(int i=0;i<n;++i){hl->SetBinContent(i+1,1000*(ns[i].birth-ns[i].exit));hl->GetXaxis()->SetBinLabel(i+1,Form("%d",ns[i].event));}
  hl->SetFillColor(kAzure-9);hl->SetLineColor(kBlue+2);hl->SetLineWidth(2);auto*c2=new TCanvas("cLoss","Energy loss",1250,750);c2->SetGridy();c2->SetLeftMargin(.16);c2->SetRightMargin(.04);c2->SetBottomMargin(.15);hl->Draw("HIST");c2->SaveAs("meeting_escaping_neutron_energy_loss.png");

  std::vector<Neutron>examples=ns;std::sort(examples.begin(),examples.end(),[](const Neutron&a,const Neutron&b){return a.birth>b.birth;});if(examples.size()>3)examples.resize(3);
  Color_t colors[]={kBlue+1,kRed+1,kGreen+2};
  for(size_t i=0;i<examples.size();++i){
    auto*c3=new TCanvas(Form("cSteps%zu",i),"Representative history",900,700);c3->SetGrid();c3->SetLeftMargin(.18);c3->SetRightMargin(.04);c3->SetBottomMargin(.17);
    auto h=history(examples[i],steps);auto*g=makeGraph(h,colors[i]);
    g->SetMinimum(0);g->SetMaximum(1.12*(*std::max_element(h.begin(),h.end())));
    g->SetTitle(Form("Event %d, track %d;Recorded point: birth #rightarrow interactions #rightarrow exit;Neutron KE [GeV]",examples[i].event,examples[i].track));g->Draw("ALP");
    c3->SaveAs(Form("meeting_neutron_steps_event%d_track%d.png",examples[i].event,examples[i].track));
  }

  const Neutron&best=examples.front();auto h=history(best,steps);auto range=std::minmax_element(h.begin(),h.end());double span=*range.second-*range.first;double padding=span>0?.1*span:std::max(.001,.05*std::abs(*range.first));
  auto*c4full=new TCanvas("cFull","Full scale",900,700);c4full->SetGrid();c4full->SetLeftMargin(.18);c4full->SetRightMargin(.04);c4full->SetBottomMargin(.16);auto*full=makeGraph(h,kOrange+1);
  full->SetMinimum(0);full->SetMaximum(1.1*(*range.second));full->SetTitle(Form("Full scale: event %d, track %d;Recorded point;Neutron KE [GeV]",best.event,best.track));full->Draw("ALP");
  c4full->SaveAs("meeting_highest_escaping_neutron_full_scale.png");
  auto*c4zoom=new TCanvas("cZoom","Zoomed energy change",900,700);c4zoom->SetGrid();c4zoom->SetLeftMargin(.25);c4zoom->SetRightMargin(.04);c4zoom->SetBottomMargin(.16);auto*zoom=makeGraph(h,kOrange+1);
  zoom->SetMinimum(std::max(0.,*range.first-padding));zoom->SetMaximum(*range.second+padding);zoom->SetTitle(Form("Zoomed energy change: event %d, track %d;Recorded point;Neutron KE [GeV]",best.event,best.track));zoom->Draw("ALP");
  zoom->GetYaxis()->SetLabelSize(.032);zoom->GetYaxis()->SetTitleOffset(2.25);c4zoom->Modified();c4zoom->Update();
  c4zoom->SaveAs("meeting_highest_escaping_neutron_zoom.png");

  std::cout<<"Selected "<<ids.size()<<" events; "<<ns.size()<<" had an escaping neutron.\n";
  std::cout<<"Highest: event "<<best.event<<", track "<<best.track<<", birth "<<best.birth<<" GeV, exit "<<best.exit<<" GeV, loss "<<1000*(best.birth-best.exit)<<" MeV.\n";
  std::cout<<"Created "<<4+examples.size()<<" meeting_*.png files, one plot per image.\n";
}
