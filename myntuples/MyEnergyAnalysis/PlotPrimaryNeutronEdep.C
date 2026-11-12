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
#include "TPaveText.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
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
  const double maximumMeV = *std::max_element(depositMeV.begin(), depositMeV.end());
  double displayScale = 1.0;
  const char *displayUnit = "MeV";
  if (maximumMeV < 0.1 && maximumMeV >= 1.0e-4)
  {
    displayScale = 1.0e3;
    displayUnit = "keV";
  }
  else if (maximumMeV < 1.0e-4)
  {
    displayScale = 1.0e6;
    displayUnit = "eV";
  }

  std::vector<double> displayDeposit(depositMeV.size());
  std::transform(depositMeV.begin(), depositMeV.end(), displayDeposit.begin(),
                 [displayScale](double value) { return value * displayScale; });
  const auto nonzero = std::count_if(displayDeposit.begin(), displayDeposit.end(),
                                     [](double value) { return value > 0.0; });
  const double maximumDisplay = maximumMeV * displayScale;
  const double yMaximum = maximumDisplay > 0.0 ? 1.18 * maximumDisplay : 1.0;

  TCanvas canvas("canvas", "canvas", 1100, 800);
  canvas.SetLeftMargin(0.13);
  canvas.SetBottomMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.09);
  canvas.SetGrid();

  TGraph graph(static_cast<int>(initialMeV.size()), initialMeV.data(), displayDeposit.data());
  graph.SetTitle(title);
  graph.GetXaxis()->SetTitle("Initial primary-neutron kinetic energy [MeV]");
  std::string yTitle = std::string("Energy deposited directly by primary-neutron track [") +
                       displayUnit + "]";
  graph.GetYaxis()->SetTitle(yTitle.c_str());
  graph.GetXaxis()->SetTitleSize(0.045);
  graph.GetYaxis()->SetTitleSize(0.045);
  graph.GetXaxis()->SetLabelSize(0.038);
  graph.GetYaxis()->SetLabelSize(0.038);
  graph.GetYaxis()->SetTitleOffset(1.35);
  graph.SetMinimum(0.0);
  graph.SetMaximum(yMaximum);
  styleGraph(graph, kBlue + 1);
  graph.Draw("AP");

  TPaveText counts(0.64, 0.74, 0.93, 0.88, "NDC");
  counts.SetFillColor(kWhite);
  counts.SetFillStyle(1001);
  counts.SetBorderSize(1);
  counts.SetTextAlign(12);
  counts.SetTextFont(42);
  counts.SetTextSize(0.032);
  std::ostringstream totalText, nonzeroText;
  totalText << "Entries: " << initialMeV.size();
  nonzeroText << "Nonzero direct deposits: " << nonzero;
  counts.AddText(totalText.str().c_str());
  counts.AddText(nonzeroText.str().c_str());
  if (maximumDisplay == 0.0)
    counts.AddText("Maximum direct deposit: 0");
  counts.Draw();

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
              "Primary neutrons",
              "primary_neutron_edep_vs_initial_ke.png");

  if (!interactingInitialMeV.empty())
    saveScatter(interactingInitialMeV, interactingDepositMeV,
                "Interacting primary neutrons",
                "interacting_primary_neutron_edep_vs_initial_ke.png");

  std::cout << "Primary neutrons plotted: " << allInitialMeV.size() << '\n'
            << "Interacting primary neutrons plotted: " << interactingInitialMeV.size() << '\n'
            << "Created primary_neutron_edep_vs_initial_ke.png\n";
  if (!interactingInitialMeV.empty())
    std::cout << "Created interacting_primary_neutron_edep_vs_initial_ke.png\n";
}
