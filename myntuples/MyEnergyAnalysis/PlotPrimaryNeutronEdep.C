// Plot direct primary-neutron energy deposition versus initial kinetic energy.
//
// Only NeutronTree rows with MotherTrackID == 0 are used.  DirectEdep contains
// SimEnergyDeposit energy assigned to the original neutron TrackID itself;
// deposits from secondary neutrons and all other descendants are excluded.
//
// Usage:
//   root -l -b -q 'PlotPrimaryNeutronEdep.C("your_output.root")'

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
TTree *findNeutronTree(TFile &file)
{
  if (auto *tree = dynamic_cast<TTree *>(file.Get("MyEnergyAnalysis/NeutronTree")))
    return tree;
  return dynamic_cast<TTree *>(file.Get("NeutronTree"));
}

void styleGraph(TGraph &graph, int color)
{
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(0.75);
  graph.SetMarkerColor(color);
  graph.SetLineColor(color);
}

void saveScatter(const std::vector<double> &initialMeV,
                 const std::vector<double> &depositMeV,
                 const char *title,
                 const char *output)
{
  TCanvas canvas("canvas", "canvas", 1100, 800);
  canvas.SetLeftMargin(0.13);
  canvas.SetBottomMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.09);
  canvas.SetGrid();

  TGraph graph(static_cast<int>(initialMeV.size()), initialMeV.data(), depositMeV.data());
  graph.SetTitle(title);
  graph.GetXaxis()->SetTitle("Initial primary-neutron kinetic energy [MeV]");
  graph.GetYaxis()->SetTitle("Energy deposited directly by primary-neutron track [MeV]");
  graph.GetXaxis()->SetTitleSize(0.045);
  graph.GetYaxis()->SetTitleSize(0.045);
  graph.GetXaxis()->SetLabelSize(0.038);
  graph.GetYaxis()->SetLabelSize(0.038);
  graph.GetYaxis()->SetTitleOffset(1.35);
  styleGraph(graph, kBlue + 1);
  graph.Draw("AP");
  canvas.SaveAs(output);
}
} // namespace

void PlotPrimaryNeutronEdep(const char *filename)
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);

  TFile file(filename, "READ");
  if (file.IsZombie())
  {
    std::cerr << "ERROR: could not open " << filename << '\n';
    return;
  }

  TTree *tree = findNeutronTree(file);
  if (!tree)
  {
    std::cerr << "ERROR: NeutronTree was not found in " << filename << '\n';
    return;
  }
  if (!tree->GetBranch("DirectEdep"))
  {
    std::cerr << "ERROR: NeutronTree has no DirectEdep branch. "
                 "Rerun MyEnergyAnalysis with the updated analyzer.\n";
    return;
  }

  int mother = 0;
  int nInteractions = 0;
  double birthKEGeV = 0.0;
  double directEdepGeV = 0.0;
  tree->SetBranchAddress("MotherTrackID", &mother);
  tree->SetBranchAddress("BirthKE", &birthKEGeV);
  tree->SetBranchAddress("DirectEdep", &directEdepGeV);
  tree->SetBranchAddress("NInteractions", &nInteractions);

  std::vector<double> allInitialMeV, allDepositMeV;
  std::vector<double> interactingInitialMeV, interactingDepositMeV;

  const Long64_t entries = tree->GetEntries();
  for (Long64_t i = 0; i < entries; ++i)
  {
    tree->GetEntry(i);
    if (mother != 0)
      continue;

    const double initialMeV = birthKEGeV * 1000.0;
    const double depositMeV = directEdepGeV * 1000.0;
    allInitialMeV.push_back(initialMeV);
    allDepositMeV.push_back(depositMeV);
    if (nInteractions > 0)
    {
      interactingInitialMeV.push_back(initialMeV);
      interactingDepositMeV.push_back(depositMeV);
    }
  }

  if (allInitialMeV.empty())
  {
    std::cerr << "ERROR: no primary neutron rows (MotherTrackID == 0) were found.\n";
    return;
  }

  saveScatter(allInitialMeV, allDepositMeV,
              "Primary neutrons;Initial primary-neutron kinetic energy [MeV];Direct primary-track energy deposit [MeV]",
              "primary_neutron_edep_vs_initial_ke.png");

  if (!interactingInitialMeV.empty())
    saveScatter(interactingInitialMeV, interactingDepositMeV,
                "Interacting primary neutrons;Initial primary-neutron kinetic energy [MeV];Direct primary-track energy deposit [MeV]",
                "interacting_primary_neutron_edep_vs_initial_ke.png");

  std::cout << "Primary neutrons plotted: " << allInitialMeV.size() << '\n'
            << "Interacting primary neutrons plotted: " << interactingInitialMeV.size() << '\n'
            << "Created primary_neutron_edep_vs_initial_ke.png\n";
  if (!interactingInitialMeV.empty())
    std::cout << "Created interacting_primary_neutron_edep_vs_initial_ke.png\n";
}
