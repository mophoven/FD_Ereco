// ROOT macro for the last two analysis requests in MyEnergyAnalysis_module.cc.
// Run with:
//   root -l -q 'Neutron20EventAnalysis.C("your_500_event_file.root",20)'

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
TTree *findTree(TDirectory *dir, const char *wanted)
{
  if (!dir) return nullptr;
  if (auto *tree = dynamic_cast<TTree *>(dir->Get(wanted))) return tree;
  TIter next(dir->GetListOfKeys());
  while (auto *key = dynamic_cast<TKey *>(next())) {
    if (std::string(key->GetClassName()).find("TDirectory") == std::string::npos) continue;
    if (auto *tree = findTree(dynamic_cast<TDirectory *>(key->ReadObj()), wanted)) return tree;
  }
  return nullptr;
}

struct Step {
  double time;
  double ke;
};
}

void Neutron20EventAnalysis(const char *fileName, int maxEvents = 20)
{
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Cannot open " << fileName << '\n';
    return;
  }

  TTree *eventTree = findTree(&input, "EventTree");
  TTree *vertexTree = findTree(&input, "VertexTree");
  TTree *escapeTree = findTree(&input, "EscapeTree");
  if (!eventTree || !vertexTree || !escapeTree) {
    std::cerr << "The file must contain EventTree, VertexTree, and EscapeTree.\n";
    return;
  }

  // Select event IDs from the file instead of assuming that they are 0,1,2,...
  int event = 0;
  eventTree->SetBranchAddress("Event", &event);
  std::vector<int> eventIds;
  std::set<int> selected;
  for (Long64_t i = 0; i < eventTree->GetEntries() && (int)eventIds.size() < maxEvents; ++i) {
    eventTree->GetEntry(i);
    if (selected.insert(event).second) eventIds.push_back(event);
  }
  eventTree->ResetBranchAddresses();
  if (eventIds.empty()) {
    std::cerr << "No events were found.\n";
    return;
  }

  // Collect incoming neutron KE at each recorded physics-interaction vertex.
  int vEvent = 0, inPdg = 0, trackId = 0;
  double inKE = 0.0, vertexTime = 0.0;
  vertexTree->SetBranchAddress("Event", &vEvent);
  vertexTree->SetBranchAddress("In_PDG", &inPdg);
  vertexTree->SetBranchAddress("In_TrackID", &trackId);
  vertexTree->SetBranchAddress("In_KE", &inKE);
  vertexTree->SetBranchAddress("Vtx_t", &vertexTime);

  std::map<int, std::map<int, std::vector<Step>>> steps;
  for (Long64_t i = 0; i < vertexTree->GetEntries(); ++i) {
    vertexTree->GetEntry(i);
    if (!selected.count(vEvent) || inPdg != 2112) continue;
    steps[vEvent][trackId].push_back({vertexTime, inKE});
  }
  for (auto &[ev, tracks] : steps)
    for (auto &[trk, values] : tracks)
      std::sort(values.begin(), values.end(), [](const Step &a, const Step &b) { return a.time < b.time; });

  gStyle->SetOptStat(0);
  auto *allCanvas = new TCanvas("cNeutron20Events", "Neutron energy at each interaction", 1800, 1200);
  allCanvas->Divide(4, 5, 0.002, 0.002);
  const int colors[] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kViolet, kGray+2};
  for (size_t pad = 0; pad < eventIds.size(); ++pad) {
    allCanvas->cd(pad + 1);
    auto *multi = new TMultiGraph();
    int color = 0;
    for (const auto &[trk, values] : steps[eventIds[pad]]) {
      if (values.empty()) continue;
      auto *graph = new TGraph(values.size());
      for (size_t s = 0; s < values.size(); ++s) graph->SetPoint(s, s + 1, values[s].ke);
      graph->SetLineColor(colors[color % 8]);
      graph->SetMarkerColor(colors[color % 8]);
      graph->SetMarkerStyle(20 + color % 5);
      graph->SetTitle(Form("track %d", trk));
      multi->Add(graph, "LP");
      ++color;
    }
    multi->SetTitle(Form("Event %d;recorded interaction step;incoming neutron KE [GeV]", eventIds[pad]));
    multi->Draw("A");
    if (color) gPad->BuildLegend(0.58, 0.68, 0.88, 0.88, "neutron tracks");
  }
  allCanvas->SaveAs("neutron_energy_steps_first20.png");
  allCanvas->SaveAs("neutron_energy_steps_first20.pdf");

  // Find the highest-birth-energy escaping neutron among those same events.
  int xEvent = 0, xTrack = 0, xPdg = 0;
  double birthKE = 0.0, exitKE = 0.0;
  escapeTree->SetBranchAddress("Event", &xEvent);
  escapeTree->SetBranchAddress("TrackID", &xTrack);
  escapeTree->SetBranchAddress("PDG", &xPdg);
  escapeTree->SetBranchAddress("BirthKE", &birthKE);
  escapeTree->SetBranchAddress("ExitKE", &exitKE);
  double bestBirth = -1.0, bestExit = 0.0;
  int bestEvent = -1, bestTrack = -1;
  for (Long64_t i = 0; i < escapeTree->GetEntries(); ++i) {
    escapeTree->GetEntry(i);
    if (!selected.count(xEvent) || xPdg != 2112 || birthKE <= bestBirth) continue;
    bestBirth = birthKE;
    bestExit = exitKE;
    bestEvent = xEvent;
    bestTrack = xTrack;
  }
  if (bestTrack < 0) {
    std::cout << "No escaping neutron was found in the selected events.\n";
    return;
  }

  std::vector<double> history{bestBirth};
  for (const Step &s : steps[bestEvent][bestTrack]) history.push_back(s.ke);
  history.push_back(bestExit);
  auto makeGraph = [&history]() {
    auto *g = new TGraph(history.size());
    for (size_t i = 0; i < history.size(); ++i) g->SetPoint(i, i, history[i]);
    g->SetLineColor(kYellow+1);
    g->SetMarkerColor(kYellow+1);
    g->SetLineWidth(3);
    g->SetMarkerStyle(20);
    return g;
  };

  auto *zoomCanvas = new TCanvas("cHighestEscapingNeutron", "Highest-energy escaping neutron", 1400, 650);
  zoomCanvas->Divide(2, 1);
  zoomCanvas->cd(1);
  auto *full = makeGraph();
  full->SetTitle(Form("Highest-energy escaping neutron: event %d, track %d;birth / interaction / exit point;KE [GeV]", bestEvent, bestTrack));
  full->Draw("ALP");

  zoomCanvas->cd(2);
  auto *zoom = makeGraph();
  const auto range = std::minmax_element(history.begin(), history.end());
  double span = *range.second - *range.first;
  double padding = (span > 0.0) ? 0.10 * span : std::max(0.01, 0.05 * std::abs(*range.first));
  zoom->SetMinimum(std::max(0.0, *range.first - padding));
  zoom->SetMaximum(*range.second + padding);
  zoom->SetTitle("Zoom on yellow curve;birth / interaction / exit point;KE [GeV]");
  zoom->Draw("ALP");
  zoomCanvas->SaveAs("highest_energy_escaping_neutron_zoom.png");
  zoomCanvas->SaveAs("highest_energy_escaping_neutron_zoom.pdf");

  std::cout << "Selected " << eventIds.size() << " events. Highest escaping neutron: event "
            << bestEvent << ", track " << bestTrack << ", birth KE " << bestBirth
            << " GeV, exit/loss energy " << bestExit << " GeV.\n";
  std::cout << "Created neutron_energy_steps_first20.* and highest_energy_escaping_neutron_zoom.*\n";
}
