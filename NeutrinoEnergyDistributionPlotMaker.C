//Imported ROOT Classes | Documentation in the ROOT Reference Guide: https://root.cern.ch/doc/master/index.html
#include<iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TString.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TAttAxis.h>
#include <TStyle.h>
using namespace std;


//Read in myntuples file
TFile *f = new TFile("/exp/dune/app/users/mophoven/EReco/srcs/myntuples/myntuples/MyEnergyAnalysis/test5files.root"); //Insert path to root file produced from MyEnergyAnalysis_module.cc
TTree *t1 = (TTree*)f->Get("MyEnergyAnalysis/MyTree");

//Defining variables

double True_LepE; //Lepton total Energy

//Kinetic energy of daughter particles
double Gen_numu_E; //Primary Neutrino Energy
double eP; //Proton Energy
double eN; //Neutron Energy
double ePi0; //Pion0 Energy
double ePip; //Pion+ Energy
double ePim; //Pion- Energy
double eOther; // Kaon+, Kaon-, Kaon0, Anti-Kaon, Kaon0-long, Kaon0-short, Gamma, IsHadron(pdg)

//Particle counter for each event
int nLep(0);
int nN(0);
int nP(0);
int nPip(0);
int nPim(0);
int nPi0(0);
double nOther(0);

int CCNC_truth; //charged current interaction detector 
int Mode_truth; //Neutrino Scattering Mechanism detector: 0=Quasi-elastic/elastic, 1=resonance production, 2=deep inelastic scattering 3=coherent production

//For plotting the energy distribution of the primary neutrino amoungst daughter particles
double Percent_energy_Lep;
double Percent_energy_N;
double Percent_energy_P;
double Percent_energy_Pip;
double Percent_energy_Pim;
double Percent_energy_Pi0;
double Percent_energy_Other;
  //For plotting the ratio of deposited energy to the kinetic energy against the primary neutrino energy
double Ratio_deposited_energy_mu;
double Ratio_deposited_energy_N;
double Ratio_deposited_energy_P;
double Ratio_deposited_energy_Pip;
double Ratio_deposited_energy_Pim;
double Ratio_deposited_energy_Pi0;
double Ratio_deposited_energy_Other;
  //Simulated deposited energy (muon,neutron,proton,pion+,pion-,pion0,Other)
double Sim_mu_Edep_b2;
double Sim_n_Edep_b2;
double Sim_p_Edep_b2;
double Sim_pip_Edep_b2;
double Sim_pim_Edep_b2;
double Sim_pi0_Edep_b2;
double Sim_Other_Edep_b2;
// Generator level outgoing lepton vtx
double Lepvtx_x, Lepvtx_y, Lepvtx_z;  


TFile outFile;

//Function will create sequence of plots with raw simulated data.
void Raw(){
    cout << "no action taken" << endl; 

    //Defining Path to Branches
    //Kinetic Energy of Particles
    t1->SetBranchAddress("eN", &eN); //neutron
    t1->SetBranchAddress("eP", &eP); //proton
    t1->SetBranchAddress("ePi0", &ePi0); //pion0
    t1->SetBranchAddress("ePip", &ePip); //pion+
    t1->SetBranchAddress("ePim", &ePim); //pion-
    t1->SetBranchAddress("eOther", &eOther); // Kaon+, Kaon-, Kaon0, Anti-Kaon, Kaon0-long, Kaon0-short, Gamma, IsHadron(pdg)
    t1->SetBranchAddress("True_LepE", &True_LepE); //lepton (muon)
    //Primary Neutrino Energy
    t1->SetBranchAddress("Gen_numu_E", &Gen_numu_E);
    //Charged Current Interaction
    t1->SetBranchAddress("CCNC_truth", &CCNC_truth);
    //Neutrino Scattering Mechanism
    t1->SetBranchAddress("Mode_truth", &Mode_truth);
    //Deposited Energy of Particles
    t1->SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2); //muon
    t1->SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2); //neutron
    t1->SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2); //proton
    t1->SetBranchAddress("Sim_pip_Edep_b2", &Sim_pip_Edep_b2); //pion+
    t1->SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2); //pion-
    t1->SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2); //pion0
    t1->SetBranchAddress("Sim_Other_Edep_b2", &Sim_Other_Edep_b2); // Kaon+, Kaon-, Kaon0, Anti-Kaon, Kaon0-long, Kaon0-short, Gamma, IsHadron(pdg)
    //number of particles per event
    t1->SetBranchAddress("nN",&nN);
    t1->SetBranchAddress("nP",&nP);
    t1->SetBranchAddress("nPip", &nPip);
    t1->SetBranchAddress("nPim", &nPim);
    t1->SetBranchAddress("nPi0", &nPi0);
    t1->SetBranchAddress("nOther",&nOther);
    //for fiducial volume cut: Generator level outgoing lepton vtx
    t1->SetBranchAddress("Lepvtx_x",&Lepvtx_x);
    t1->SetBranchAddress("Lepvtx_y",&Lepvtx_y);
    t1->SetBranchAddress("Lepvtx_z",&Lepvtx_z);


    //2D plots for fraction of kinetic energy/neutrino energy over neutrino energy
    TH2F *plot_2D_Lep_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by Lepton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_N_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by Neutron", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_P_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by Proton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pip_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by #pi^{+}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pim_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by #pi^{-}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pi0_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by #pi^{0}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Other_raw = new TH2F("Statistics", "Raw - Fraction of Neutrino Energy Received by Other", 50, 0, 7, 50, 0, 1);

    //2D plots for ratio of deposited energy to kinetic energy over neutrino energy
    TH2F *plot_2D_mu_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Lepton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Neutron Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Proton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Pion+ Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Pion- Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Pion0 Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_Edep_raw = new TH2F("Statistics", "Raw - Fraction of Energy of Other Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);

    //2D plots for ratio of deposited energy to kinetic energy over kinetic energy of particle
    TH2F *plot_2D_mu_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Lepton Deposited in Detector Lepton varying Lepton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Neutron Deposited in Detector varying Neutron Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Proton Deposited in Detector varying Proton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Pion+ Deposited in Detector varying Pion+ Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Pion- Deposited in Detector varying Pion- Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Pion0 Deposited in Detector varying Pion0 Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_EdepKE_raw = new TH2F("Statistics", "Fraction of Energy of Other Deposited in Detector varying Other Energy", 50, 0, 7, 50, 0, 1.4);


    // Loop over all events
    int nentries = 0; // Total input events
    nentries = t1->GetEntries();
    cout<< "nentries: " << nentries<<endl;


    for ( int ientry = 0; ientry < nentries; ientry++ )
    {
    t1->GetEntry(ientry);

    //fiducial volume cut
    if (abs(Lepvtx_x) > 310 ) continue;
    if (abs(Lepvtx_y) > 550) continue;
    if (50 > Lepvtx_z) continue;
    if (Lepvtx_z > 1244) continue;

    if (CCNC_truth == 0) { // if 0 then charged current interaction, if 1 then neutral current interaction    

    //Converting Energy Values from MeV to GeV
    Sim_mu_Edep_b2 = Sim_mu_Edep_b2/1000;
    Sim_n_Edep_b2 = Sim_n_Edep_b2/1000;
    Sim_p_Edep_b2 = Sim_p_Edep_b2/1000;
    Sim_pip_Edep_b2 = Sim_pip_Edep_b2/1000;
    Sim_pim_Edep_b2 = Sim_pim_Edep_b2/1000;
    Sim_pi0_Edep_b2 = Sim_pi0_Edep_b2/1000;
    Sim_Other_Edep_b2 = Sim_Other_Edep_b2/1000;    

    //multiply masses by number of particle
    //eP = eP + 0.938272; //adding proton mass to KE to get Energy in GeV
    //eN = eN + 0.939565; //adding neutron mass to KE to get Energy in GeV
    ePi0 = ePi0 + (0.134977)*(nPi0); //adding pion 0 mass to KE to get Energy in GeV
    ePip = ePip + (0.139571)*(nPip); //adding pion+ mass to KE to get Energy in Gev
    ePim = ePim + (0.139571)*(nPim); //adding pion- mass to KE to get Energy in Gev


    Percent_energy_Lep = True_LepE/Gen_numu_E;
    Percent_energy_N = eN/Gen_numu_E;
    Percent_energy_P = eP/Gen_numu_E;
    Percent_energy_Pip = ePip/Gen_numu_E;
    Percent_energy_Pim = ePim/Gen_numu_E;
    Percent_energy_Pi0 = ePi0/Gen_numu_E;
    Percent_energy_Other = eOther/Gen_numu_E;

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;

    //Energy Distribution Plots from Primary Neutrino Amongst Particles
    plot_2D_N_raw->Fill(Gen_numu_E, Percent_energy_N);
    plot_2D_P_raw->Fill(Gen_numu_E, Percent_energy_P);
    plot_2D_Pip_raw->Fill(Gen_numu_E, Percent_energy_Pip);
    plot_2D_Pim_raw->Fill(Gen_numu_E, Percent_energy_Pim);
    plot_2D_Pi0_raw->Fill(Gen_numu_E, Percent_energy_Pi0);
    plot_2D_Other_raw->Fill(Gen_numu_E, Percent_energy_Other);

    plot_2D_Lep_raw->Fill(Gen_numu_E, Percent_energy_Lep); //(x,y)
    plot_2D_Pi0_EdepKE_raw->Fill(ePi0, Ratio_deposited_energy_Pi0); //need pi0 to use all energy not only KE so moving above redefinition

    //Redefining Pion Energies to KE
    ePi0 = ePi0; // KE in Gev
    ePip = ePip; // KE in Gev
    ePim = ePim; // KE in Gev

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;

    //Energy Deposited in Detector Amongst Particles Across Primary Neutrino Energies Plots
    plot_2D_mu_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_mu); //(x,y)
    plot_2D_N_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_N);
    plot_2D_P_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_P);
    plot_2D_Pip_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_Pip);
    plot_2D_Pim_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_Pim);
    plot_2D_Pi0_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_Pi0);
    plot_2D_Other_Edep_raw->Fill(Gen_numu_E, Ratio_deposited_energy_Other);

    //Energy Deposited in Detector Amongst Particles Across Particle Energies Plots
    plot_2D_mu_EdepKE_raw->Fill(True_LepE, Ratio_deposited_energy_mu); //(x,y)
    plot_2D_N_EdepKE_raw->Fill(eN, Ratio_deposited_energy_N);
    plot_2D_P_EdepKE_raw->Fill(eP, Ratio_deposited_energy_P);
    plot_2D_Pip_EdepKE_raw->Fill(ePip, Ratio_deposited_energy_Pip);
    plot_2D_Pim_EdepKE_raw->Fill(ePim, Ratio_deposited_energy_Pim);
    plot_2D_Other_EdepKE_raw->Fill(eOther, Ratio_deposited_energy_Other);
    }//end CCNC_truth
    }// end ientry
    //Create Canvas,Draw Plot, Label Axes, Save
  
    //True Energy to Primary Neutrino Energy (Percent Energy) against Primary Neutrino Energy
    TCanvas *c1 = new TCanvas("c1", "c1",1800, 1350);
    c1->cd();
    plot_2D_Lep_raw->Draw("COLZ");
    plot_2D_Lep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Lep_raw->GetYaxis()->SetTitle("Lepton Total E/Primary #nu Energy "); //Set title of y-axis
    c1->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyLepton2D.jpg"); //Set file name
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 1350);
    c2->cd();
    plot_2D_N_raw->Draw("COLZ");
    plot_2D_N_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_raw->GetYaxis()->SetTitle("Neutron KE/Primary #nu Energy "); //Set title of y-axis
    c2->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyNeutron2D.jpg"); //Set file name
    
    
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 1350);
    c3->cd();  
    plot_2D_P_raw->Draw("COLZ");
    plot_2D_P_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_raw->GetYaxis()->SetTitle("Proton KE/Primary #nu Energy "); //Set title of y-axis
    c3->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 1350);
    c4->cd();  
    plot_2D_Pip_raw->Draw("COLZ");
    plot_2D_Pip_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_raw->GetYaxis()->SetTitle("Pion+ total E/Primary #nu Energy "); //Set title of y-axis
    c4->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 1350);
    c5->cd();  
    plot_2D_Pim_raw->Draw("COLZ");
    plot_2D_Pim_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_raw->GetYaxis()->SetTitle("Pion- total E/Primary #nu Energy "); //Set title of y-axis
    c5->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 1350);
    c6->cd();  
    plot_2D_Pi0_raw->Draw("COLZ");
    plot_2D_Pi0_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_raw->GetYaxis()->SetTitle("Pion0 total E/Primary #nu Energy "); //Set title of y-axis
    c6->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c7 = new TCanvas("c7", "c7", 1800, 1350);
    c7->cd();  
    plot_2D_Other_raw->Draw("COLZ");
    plot_2D_Other_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_raw->GetYaxis()->SetTitle("Other Energy/Primary #nu Energy "); //Set title of y-axis
    c7->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/PercentEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Primary Neutrino Energy

    TCanvas *c8 = new TCanvas("c8", "c8", 1800, 1350);
    c8->cd();  
    plot_2D_mu_Edep_raw->Draw("COLZ");
    plot_2D_mu_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_mu_Edep_raw->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c8->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergymuon2D.jpg"); //Set file name

    TCanvas *c9 = new TCanvas("c9", "c9", 1800, 1350);
    c9->cd();
    plot_2D_N_Edep_raw->Draw("COLZ");
    plot_2D_N_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_Edep_raw->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c9->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyNeutron2D.jpg"); //Set file name
    
    TCanvas *c10 = new TCanvas("c10", "c10", 1800, 1350);
    c10->cd();  
    plot_2D_P_Edep_raw->Draw("COLZ");
    plot_2D_P_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_Edep_raw->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c10->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c11 = new TCanvas("c11", "c11", 1800, 1350);
    c11->cd();  
    plot_2D_Pip_Edep_raw->Draw("COLZ");
    plot_2D_Pip_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_Edep_raw->GetYaxis()->SetTitle("Pion+ Edep/Pi+ KE "); //Set title of y-axis
    c11->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c12 = new TCanvas("c12", "c12", 1800, 1350);
    c12->cd();  
    plot_2D_Pim_Edep_raw->Draw("COLZ");
    plot_2D_Pim_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_Edep_raw->GetYaxis()->SetTitle("Pion- Edep/Pi- KE "); //Set title of y-axis
    c12->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c13 = new TCanvas("c13", "c13", 1800, 1350);
    c13->cd();  
    plot_2D_Pi0_Edep_raw->Draw("COLZ");
    plot_2D_Pi0_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_Edep_raw->GetYaxis()->SetTitle("Pion0 Edep/Pi0 KE "); //Set title of y-axis
    c13->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c14 = new TCanvas("c14", "c14", 1800, 1350);
    c14->cd();  
    plot_2D_Other_Edep_raw->Draw("COLZ");
    plot_2D_Other_Edep_raw->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_Edep_raw->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c14->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Kinetic Energy

    TCanvas *c15 = new TCanvas("c15", "c15", 1800, 1350);
    c15->cd();  
    plot_2D_mu_EdepKE_raw->Draw("COLZ");
    plot_2D_mu_EdepKE_raw->GetXaxis()->SetTitle("Lepton Total E (GeV)"); //Set title of x-axis
    plot_2D_mu_EdepKE_raw->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c15->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEmuon2D.jpg"); //Set file name

    TCanvas *c16 = new TCanvas("c16", "c16", 1800, 1350);
    c16->cd();
    plot_2D_N_EdepKE_raw->Draw("COLZ");
    plot_2D_N_EdepKE_raw->GetXaxis()->SetTitle("Neutron KE (GeV)"); //Set title of x-axis
    plot_2D_N_EdepKE_raw->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c16->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKENeutron2D.jpg"); //Set file name
    
    TCanvas *c17 = new TCanvas("c17", "c17", 1800, 1350);
    c17->cd();  
    plot_2D_P_EdepKE_raw->Draw("COLZ");
    plot_2D_P_EdepKE_raw->GetXaxis()->SetTitle("Proton KE (GeV)"); //Set title of x-axis
    plot_2D_P_EdepKE_raw->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c17->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEProton2D.jpg"); //Set file name
    
    TCanvas *c18 = new TCanvas("c18", "c18", 1800, 1350);
    c18->cd();  
    plot_2D_Pip_EdepKE_raw->Draw("COLZ");
    plot_2D_Pip_EdepKE_raw->GetXaxis()->SetTitle("Pi+ KE (GeV)"); //Set title of x-axis
    plot_2D_Pip_EdepKE_raw->GetYaxis()->SetTitle("Pi+ Edep/Pi+ KE "); //Set title of y-axis
    c18->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEPion+2D.jpg"); //Set file name  
    
    TCanvas *c19 = new TCanvas("c19", "c19", 1800, 1350);
    c19->cd();  
    plot_2D_Pim_EdepKE_raw->Draw("COLZ");
    plot_2D_Pim_EdepKE_raw->GetXaxis()->SetTitle("Pi- KE (GeV)"); //Set title of x-axis
    plot_2D_Pim_EdepKE_raw->GetYaxis()->SetTitle("Pi- Edep/Pi- KE "); //Set title of y-axis
    c19->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEPion-2D.jpg"); //Set file name  
    
    TCanvas *c20 = new TCanvas("c20", "c20", 1800, 1350);
    c20->cd();  
    plot_2D_Pi0_EdepKE_raw->Draw("COLZ");
    plot_2D_Pi0_EdepKE_raw->GetXaxis()->SetTitle("Pi0 total energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_EdepKE_raw->GetYaxis()->SetTitle("Pi0 Edep/Pi+ total energy "); //Set title of y-axis
    c20->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEPion02D.jpg"); //Set file name  
    
    TCanvas *c21 = new TCanvas("c21", "c21", 1800, 1350);
    c21->cd();  
    plot_2D_Other_EdepKE_raw->Draw("COLZ");
    plot_2D_Other_EdepKE_raw->GetXaxis()->SetTitle("Other KE (GeV)"); //Set title of x-axis
    plot_2D_Other_EdepKE_raw->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c21->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/raw/DepositedEnergyKEOther2D.jpg"); //Set file name
    

outFile.Close();
}



//Function will create sequence of plots with simulated data column normalized to 1.
void Column_Normalized(){
    cout << "plots will be column normalized" << endl;

    //Defining Path to Branches
    //Kinetic Energy of Particles
    t1->SetBranchAddress("eN", &eN); //neutron
    t1->SetBranchAddress("eP", &eP); //proton
    t1->SetBranchAddress("ePi0", &ePi0); //pion0
    t1->SetBranchAddress("ePip", &ePip); //pion+
    t1->SetBranchAddress("ePim", &ePim); //pion-
    t1->SetBranchAddress("eOther", &eOther);
    t1->SetBranchAddress("True_LepE", &True_LepE); //lepton (muon)
    //Primary Neutrino Energy
    t1->SetBranchAddress("Gen_numu_E", &Gen_numu_E);
    //Charged Current Interaction
    t1->SetBranchAddress("CCNC_truth", &CCNC_truth);
    //Neutrino Scattering Mechanism
    t1->SetBranchAddress("Mode_truth", &Mode_truth);
    //Deposited Energy of Particles
    t1->SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2); //muon
    t1->SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2); //neutron
    t1->SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2); //proton
    t1->SetBranchAddress("Sim_pip_Edep_b2", &Sim_pip_Edep_b2); //pion+
    t1->SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2); //pion-
    t1->SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2); //pion0
    t1->SetBranchAddress("Sim_Other_Edep_b2", &Sim_Other_Edep_b2); 
    //number of particles per event
    t1->SetBranchAddress("nN",&nN);
    t1->SetBranchAddress("nP",&nP);
    t1->SetBranchAddress("nPip", &nPip);
    t1->SetBranchAddress("nPim", &nPim);
    t1->SetBranchAddress("nPi0", &nPi0);
    t1->SetBranchAddress("nOther",&nOther);
    //for fiducial volume cut: Generator level outgoing lepton vtx
    t1->SetBranchAddress("Lepvtx_x",&Lepvtx_x);
    t1->SetBranchAddress("Lepvtx_y",&Lepvtx_y);
    t1->SetBranchAddress("Lepvtx_z",&Lepvtx_z);


    //2D plots for fraction of kinetic energy/neutrino energy over neutrino energy
    TH2F *plot_2D_Lep_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Lepton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_N_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Neutron", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_P_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Proton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pip_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Pion+", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pim_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Pion-", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pi0_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Pion0", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Other_cn = new TH2F("Statistics", "cn Fraction of Neutrino Energy Received by Other", 50, 0, 7, 50, 0, 1);

    //2D plots for ratio of deposited energy to kinetic energy over neutrino energy
    TH2F *plot_2D_mu_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Lepton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Neutron Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Proton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion+ Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion- Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion0 Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_Edep_cn = new TH2F("Statistics", "cn Fraction of Energy of Other Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);

    //2D plots for ratio of deposited energy to kinetic energy over kinetic energy of particle
    TH2F *plot_2D_mu_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Lepton Deposited in Detector varying Lepton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Neutron Deposited in Detector varying Neutron Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Proton Deposited in Detector varying Proton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion+ Deposited in Detector varying Pion+ Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion- Deposited in Detector varying Pion- Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Pion0 Deposited in Detector varying Pion0 Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_EdepKE_cn = new TH2F("Statistics", "cn Fraction of Energy of Other Deposited in Detector varying Other Energy", 50, 0, 7, 50, 0, 1.4);


    // Loop over all events
    int nentries = 0; // Total input events
    nentries = t1->GetEntries();
    cout<< "nentries: " << nentries<<endl;


    for ( int ientry = 0; ientry < nentries; ientry++ )
    {
    t1->GetEntry(ientry);

    //fiducial volume cut
    if (abs(Lepvtx_x) > 310 ) continue;
    if (abs(Lepvtx_y) > 550) continue;
    if (50 > Lepvtx_z) continue;
    if (Lepvtx_z > 1244) continue;

    if (CCNC_truth == 0) { // if 0 then charged current interaction, if 1 then neutral current interaction    

    //Converting Energy Values from MeV to GeV
    Sim_mu_Edep_b2 = Sim_mu_Edep_b2/1000;
    Sim_n_Edep_b2 = Sim_n_Edep_b2/1000;
    Sim_p_Edep_b2 = Sim_p_Edep_b2/1000;
    Sim_pip_Edep_b2 = Sim_pip_Edep_b2/1000;
    Sim_pim_Edep_b2 = Sim_pim_Edep_b2/1000;
    Sim_pi0_Edep_b2 = Sim_pi0_Edep_b2/1000;
    Sim_Other_Edep_b2 = Sim_Other_Edep_b2/1000;    

    //multiply masses by number of particle
    //eP = eP + 0.938272; //adding proton mass to KE to get Energy in GeV
    //eN = eN + 0.939565; //adding neutron mass to KE to get Energy in GeV
    ePi0 = ePi0 + (0.134977)*(nPi0); //adding pion 0 mass to KE to get Energy in GeV
    ePip = ePip + (0.139571)*(nPip); //adding pion+ mass to KE to get Energy in Gev
    ePim = ePim + (0.139571)*(nPim); //adding pion- mass to KE to get Energy in Gev


    Percent_energy_Lep = True_LepE/Gen_numu_E;
    Percent_energy_N = eN/Gen_numu_E;
    Percent_energy_P = eP/Gen_numu_E;
    Percent_energy_Pip = ePip/Gen_numu_E;
    Percent_energy_Pim = ePim/Gen_numu_E;
    Percent_energy_Pi0 = ePi0/Gen_numu_E;
    Percent_energy_Other = eOther/Gen_numu_E;

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;

    //Energy Distribution Plots from Primary Neutrino Amongst Particles
    plot_2D_Lep_cn->Fill(Gen_numu_E, Percent_energy_Lep); //(x,y)
    plot_2D_N_cn->Fill(Gen_numu_E, Percent_energy_N);
    plot_2D_P_cn->Fill(Gen_numu_E, Percent_energy_P);
    plot_2D_Pip_cn->Fill(Gen_numu_E, Percent_energy_Pip);
    plot_2D_Pim_cn->Fill(Gen_numu_E, Percent_energy_Pim);
    plot_2D_Pi0_cn->Fill(Gen_numu_E, Percent_energy_Pi0);
    plot_2D_Other_cn->Fill(Gen_numu_E, Percent_energy_Other);

    plot_2D_Pi0_EdepKE_cn->Fill(ePi0, Ratio_deposited_energy_Pi0); //need pi0 to use all energy not only KE so moving above redefinition

    //Redefining Pion Energies to KE
    ePi0 = ePi0; // KE in Gev
    ePip = ePip; // KE in Gev
    ePim = ePim; // KE in Gev

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;

    //Energy Deposited in Detector Amongst Particles Across Primary Neutrino Energies Plots
    plot_2D_mu_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_mu); //(x,y)
    plot_2D_N_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_N);
    plot_2D_P_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_P);
    plot_2D_Pip_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_Pip);
    plot_2D_Pim_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_Pim);
    plot_2D_Pi0_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_Pi0);
    plot_2D_Other_Edep_cn->Fill(Gen_numu_E, Ratio_deposited_energy_Other);

    //Energy Deposited in Detector Amongst Particles Across Particle Energies Plots
    plot_2D_mu_EdepKE_cn->Fill(True_LepE, Ratio_deposited_energy_mu); //(x,y)
    plot_2D_N_EdepKE_cn->Fill(eN, Ratio_deposited_energy_N);
    plot_2D_P_EdepKE_cn->Fill(eP, Ratio_deposited_energy_P);
    plot_2D_Pip_EdepKE_cn->Fill(ePip, Ratio_deposited_energy_Pip);
    plot_2D_Pim_EdepKE_cn->Fill(ePim, Ratio_deposited_energy_Pim);
    plot_2D_Other_EdepKE_cn->Fill(eOther, Ratio_deposited_energy_Other);
    }//end CCNC_truth
    }// end ientry
  
    //Column Normalizing Plots to 1 
    cout << "Column Normalizing to 1..." << endl;
    //Percent Energy plots
    for (int ix=1; ix<=plot_2D_Lep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Lep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Lep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Lep_cn->GetNbinsY(); iy++){
        plot_2D_Lep_cn->SetBinContent(ix, iy, plot_2D_Lep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_cn->GetNbinsY(); iy++){
        plot_2D_N_cn->SetBinContent(ix, iy, plot_2D_N_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_cn->GetNbinsY(); iy++){
        plot_2D_P_cn->SetBinContent(ix, iy, plot_2D_P_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_cn->GetNbinsY(); iy++){
        plot_2D_Pip_cn->SetBinContent(ix, iy, plot_2D_Pip_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_cn->GetNbinsY(); iy++){
        plot_2D_Pim_cn->SetBinContent(ix, iy, plot_2D_Pim_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_cn->GetNbinsY(); iy++){
        plot_2D_Pi0_cn->SetBinContent(ix, iy, plot_2D_Pi0_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_cn->GetNbinsY(); iy++){
        plot_2D_Other_cn->SetBinContent(ix, iy, plot_2D_Other_cn->GetBinContent(ix, iy)/colnorm);
    }
    }
    //Deposited Energy Plots
    for (int ix=1; ix<=plot_2D_mu_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_mu_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_mu_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_mu_Edep_cn->GetNbinsY(); iy++){
        plot_2D_mu_Edep_cn->SetBinContent(ix, iy, plot_2D_mu_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_Edep_cn->GetNbinsY(); iy++){
        plot_2D_N_Edep_cn->SetBinContent(ix, iy, plot_2D_N_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_Edep_cn->GetNbinsY(); iy++){
        plot_2D_P_Edep_cn->SetBinContent(ix, iy, plot_2D_P_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_Edep_cn->GetNbinsY(); iy++){
        plot_2D_Pip_Edep_cn->SetBinContent(ix, iy, plot_2D_Pip_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_Edep_cn->GetNbinsY(); iy++){
        plot_2D_Pim_Edep_cn->SetBinContent(ix, iy, plot_2D_Pim_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_Edep_cn->GetNbinsY(); iy++){
        plot_2D_Pi0_Edep_cn->SetBinContent(ix, iy, plot_2D_Pi0_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_Edep_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_Edep_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_Edep_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_Edep_cn->GetNbinsY(); iy++){
        plot_2D_Other_Edep_cn->SetBinContent(ix, iy, plot_2D_Other_Edep_cn->GetBinContent(ix, iy)/colnorm);
    }
    }  
    //Deposited Energy/Kinetic Energy Plots
    for (int ix=1; ix<=plot_2D_mu_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_mu_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_mu_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_mu_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_mu_EdepKE_cn->SetBinContent(ix, iy, plot_2D_mu_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_N_EdepKE_cn->SetBinContent(ix, iy, plot_2D_N_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_P_EdepKE_cn->SetBinContent(ix, iy, plot_2D_P_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_Pip_EdepKE_cn->SetBinContent(ix, iy, plot_2D_Pip_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_Pim_EdepKE_cn->SetBinContent(ix, iy, plot_2D_Pim_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_Pi0_EdepKE_cn->SetBinContent(ix, iy, plot_2D_Pi0_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_EdepKE_cn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_EdepKE_cn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_EdepKE_cn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_EdepKE_cn->GetNbinsY(); iy++){
        plot_2D_Other_EdepKE_cn->SetBinContent(ix, iy, plot_2D_Other_EdepKE_cn->GetBinContent(ix, iy)/colnorm);
    }
    }
    
    //Create Canvas,Draw Plot, Label Axes, Save
    //Kinetic Energy to Primary Neutrino Energy (Percent Energy) against Primary Neutrino Energy
    TCanvas *c1 = new TCanvas("c1", "c1",1800, 1350);
    c1->cd();
    gStyle->SetOptStat(11);
    plot_2D_Lep_cn->Draw("COLZ");
    plot_2D_Lep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Lep_cn->GetYaxis()->SetTitle("Lepton Total E/Primary #nu Energy "); //Set title of y-axis
    c1->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyLepton2D.jpg"); //Set file name
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 1350);
    c2->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_cn->Draw("COLZ");
    plot_2D_N_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_cn->GetYaxis()->SetTitle("Neutron KE/Primary #nu Energy "); //Set title of y-axis
    c2->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/fnal_data/columnnormalized/PercentEnergyNeutron2D.jpg"); //Set file name
    
    
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 1350);
    c3->cd();  
    gStyle->SetOptStat(11);
    plot_2D_P_cn->Draw("COLZ");
    plot_2D_P_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_cn->GetYaxis()->SetTitle("Proton KE/Primary #nu Energy "); //Set title of y-axis
    c3->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 1350);
    c4->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_cn->Draw("COLZ");
    plot_2D_Pip_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_cn->GetYaxis()->SetTitle("Pion+ total E/Primary #nu Energy "); //Set title of y-axis
    c4->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 1350);
    c5->cd(); 
    gStyle->SetOptStat(11); 
    plot_2D_Pim_cn->Draw("COLZ");
    plot_2D_Pim_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_cn->GetYaxis()->SetTitle("Pion- total E/Primary #nu Energy "); //Set title of y-axis
    c5->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 1350);
    c6->cd(); 
    gStyle->SetOptStat(11); 
    plot_2D_Pi0_cn->Draw("COLZ");
    plot_2D_Pi0_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_cn->GetYaxis()->SetTitle("Pion0 total E/Primary #nu Energy "); //Set title of y-axis
    c6->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c7 = new TCanvas("c7", "c7", 1800, 1350);
    c7->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_cn->Draw("COLZ");
    plot_2D_Other_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_cn->GetYaxis()->SetTitle("Other Energy/Primary #nu Energy "); //Set title of y-axis
    c7->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/PercentEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Primary Neutrino Energy

    TCanvas *c8 = new TCanvas("c8", "c8", 1800, 1350);
    c8->cd();  
    gStyle->SetOptStat(11);
    plot_2D_mu_Edep_cn->Draw("COLZ");
    plot_2D_mu_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_mu_Edep_cn->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c8->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergymuon2D.jpg"); //Set file name

    TCanvas *c9 = new TCanvas("c9", "c9", 1800, 1350);
    c9->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_Edep_cn->Draw("COLZ");
    plot_2D_N_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_Edep_cn->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c9->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyNeutron2D.jpg"); //Set file name
    
    TCanvas *c10 = new TCanvas("c10", "c10", 1800, 1350);
    c10->cd();  
    gStyle->SetOptStat(11);
    plot_2D_P_Edep_cn->Draw("COLZ");
    plot_2D_P_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_Edep_cn->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c10->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c11 = new TCanvas("c11", "c11", 1800, 1350);
    c11->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_Edep_cn->Draw("COLZ");
    plot_2D_Pip_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_Edep_cn->GetYaxis()->SetTitle("Pion+ Edep/Pi+ KE "); //Set title of y-axis
    c11->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c12 = new TCanvas("c12", "c12", 1800, 1350);
    c12->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pim_Edep_cn->Draw("COLZ");
    plot_2D_Pim_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_Edep_cn->GetYaxis()->SetTitle("Pion- Edep/Pi- KE "); //Set title of y-axis
    c12->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c13 = new TCanvas("c13", "c13", 1800, 1350);
    c13->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pi0_Edep_cn->Draw("COLZ");
    plot_2D_Pi0_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_Edep_cn->GetYaxis()->SetTitle("Pion0 Edep/Pi0 KE "); //Set title of y-axis
    c13->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c14 = new TCanvas("c14", "c14", 1800, 1350);
    c14->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_Edep_cn->Draw("COLZ");
    plot_2D_Other_Edep_cn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_Edep_cn->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c14->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Kinetic Energy

    TCanvas *c15 = new TCanvas("c15", "c15", 1800, 1350);
    c15->cd();  
    gStyle->SetOptStat(11);
    plot_2D_mu_EdepKE_cn->Draw("COLZ");
    plot_2D_mu_EdepKE_cn->GetXaxis()->SetTitle("Lepton Total E (GeV)"); //Set title of x-axis
    plot_2D_mu_EdepKE_cn->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c15->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEmuon2D.jpg"); //Set file name

    TCanvas *c16 = new TCanvas("c16", "c16", 1800, 1350);
    c16->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_EdepKE_cn->Draw("COLZ");
    plot_2D_N_EdepKE_cn->GetXaxis()->SetTitle("Neutron KE (GeV)"); //Set title of x-axis
    plot_2D_N_EdepKE_cn->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c16->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKENeutron2D.jpg"); //Set file name
    
    TCanvas *c17 = new TCanvas("c17", "c17", 1800, 1350);
    c17->cd();  
    gStyle->SetOptStat(11);
    plot_2D_P_EdepKE_cn->Draw("COLZ");
    plot_2D_P_EdepKE_cn->GetXaxis()->SetTitle("Proton KE (GeV)"); //Set title of x-axis
    plot_2D_P_EdepKE_cn->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c17->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEProton2D.jpg"); //Set file name
    
    TCanvas *c18 = new TCanvas("c18", "c18", 1800, 1350);
    c18->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_EdepKE_cn->Draw("COLZ");
    plot_2D_Pip_EdepKE_cn->GetXaxis()->SetTitle("Pi+ KE (GeV)"); //Set title of x-axis
    plot_2D_Pip_EdepKE_cn->GetYaxis()->SetTitle("Pi+ Edep/Pi+ KE "); //Set title of y-axis
    c18->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEPion+2D.jpg"); //Set file name  
    
    TCanvas *c19 = new TCanvas("c19", "c19", 1800, 1350);
    c19->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pim_EdepKE_cn->Draw("COLZ");
    plot_2D_Pim_EdepKE_cn->GetXaxis()->SetTitle("Pi- KE (GeV)"); //Set title of x-axis
    plot_2D_Pim_EdepKE_cn->GetYaxis()->SetTitle("Pi- Edep/Pi- KE "); //Set title of y-axis
    c19->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEPion-2D.jpg"); //Set file name  
    
    TCanvas *c20 = new TCanvas("c20", "c20", 1800, 1350);
    c20->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pi0_EdepKE_cn->Draw("COLZ");
    plot_2D_Pi0_EdepKE_cn->GetXaxis()->SetTitle("Pi0 total energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_EdepKE_cn->GetYaxis()->SetTitle("Pi0 Edep/Pi+ total energy "); //Set title of y-axis
    c20->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEPion02D.jpg"); //Set file name  
    
    TCanvas *c21 = new TCanvas("c21", "c21", 1800, 1350);
    c21->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_EdepKE_cn->Draw("COLZ");
    plot_2D_Other_EdepKE_cn->GetXaxis()->SetTitle("Other KE (GeV)"); //Set title of x-axis
    plot_2D_Other_EdepKE_cn->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c21->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/columnnormalized/DepositedEnergyKEOther2D.jpg"); //Set file name
    

outFile.Close();
}



//Function will create sequence of plots with simulated data zero-suppressed (i.e. empty data points are not plotted).
void Zero_Suppressed(){
    cout << "plots will be zero suppressed" << endl;

    //Defining Path to Branches
    //Kinetic Energy of Particles
    t1->SetBranchAddress("eN", &eN); //neutron
    t1->SetBranchAddress("eP", &eP); //proton
    t1->SetBranchAddress("ePi0", &ePi0); //pion0
    t1->SetBranchAddress("ePip", &ePip); //pion+
    t1->SetBranchAddress("ePim", &ePim); //pion-
    t1->SetBranchAddress("eOther", &eOther);
    t1->SetBranchAddress("True_LepE", &True_LepE); //lepton (muon)
    //Primary Neutrino Energy
    t1->SetBranchAddress("Gen_numu_E", &Gen_numu_E);
    //Charged Current Interaction
    t1->SetBranchAddress("CCNC_truth", &CCNC_truth);
    //Neutrino Scattering Mechanism
    t1->SetBranchAddress("Mode_truth", &Mode_truth);
    //Deposited Energy of Particles
    t1->SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2); //muon
    t1->SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2); //neutron
    t1->SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2); //proton
    t1->SetBranchAddress("Sim_pip_Edep_b2", &Sim_pip_Edep_b2); //pion+
    t1->SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2); //pion-
    t1->SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2); //pion0
    t1->SetBranchAddress("Sim_Other_Edep_b2", &Sim_Other_Edep_b2); 
    //number of particles per event
    t1->SetBranchAddress("nLep",&nLep);
    t1->SetBranchAddress("nN",&nN);
    t1->SetBranchAddress("nP",&nP);
    t1->SetBranchAddress("nPip", &nPip);
    t1->SetBranchAddress("nPim", &nPim);
    t1->SetBranchAddress("nPi0", &nPi0);
    t1->SetBranchAddress("nOther",&nOther);
    //for fiducial volume cut: Generator level outgoing lepton vtx
    t1->SetBranchAddress("Lepvtx_x",&Lepvtx_x);
    t1->SetBranchAddress("Lepvtx_y",&Lepvtx_y);
    t1->SetBranchAddress("Lepvtx_z",&Lepvtx_z);


    //2D plots for fraction of kinetic energy/neutrino energy over neutrino energy
    TH2F *plot_2D_Lep_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Lepton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_N_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Neutron", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_P_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Proton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pip_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Pion+", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pim_zs= new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Pion-", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pi0_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Pion0", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Other_zs = new TH2F("Statistics", "zs Fraction of Neutrino Energy Received by Other", 50, 0, 7, 50, 0, 1);

    //2D plots for ratio of deposited energy to kinetic energy over neutrino energy
    TH2F *plot_2D_mu_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Lepton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Neutron Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Proton Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion+ Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion- Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion0 Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_Edep_zs = new TH2F("Statistics", "zs Fraction of Energy of Other Deposited in Detector varying Neutrino Energy", 50, 0, 7, 50, 0, 1.4);

    //2D plots for ratio of deposited energy to kinetic energy over kinetic energy of particle
    TH2F *plot_2D_mu_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Lepton Deposited in Detector varying Lepton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Neutron Deposited in Detector varying Neutron Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Proton Deposited in Detector varying Proton Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion+ Deposited in Detector varying Pion+ Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion- Deposited in Detector varying Pion- Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Pion0 Deposited in Detector varying Pion0 Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_EdepKE_zs = new TH2F("Statistics", "zs Fraction of Energy of Other Deposited in Detector varying Other Energy", 50, 0, 7, 50, 0, 1.4);


    // Loop over all events
    int nentries = 0; // Total input events
    nentries = t1->GetEntries();
    cout<< "nentries: " << nentries<<endl;

    cout << "zero-suppressing plots..." << endl;

    for ( int ientry = 0; ientry < nentries; ientry++ )
    {
    t1->GetEntry(ientry);

    //fiducial volume cut
    if (abs(Lepvtx_x) > 310 ) continue;
    if (abs(Lepvtx_y) > 550) continue;
    if (50 > Lepvtx_z) continue;
    if (Lepvtx_z > 1244) continue;

    if (CCNC_truth == 0) { // if 0 then charged current interaction, if 1 then neutral current interaction    

    //Converting Energy Values from MeV to GeV
    Sim_mu_Edep_b2 = Sim_mu_Edep_b2/1000;
    Sim_n_Edep_b2 = Sim_n_Edep_b2/1000;
    Sim_p_Edep_b2 = Sim_p_Edep_b2/1000;
    Sim_pip_Edep_b2 = Sim_pip_Edep_b2/1000;
    Sim_pim_Edep_b2 = Sim_pim_Edep_b2/1000;
    Sim_pi0_Edep_b2 = Sim_pi0_Edep_b2/1000;
    Sim_Other_Edep_b2 = Sim_Other_Edep_b2/1000;    

    //multiply masses by number of particle
    //eP = eP + 0.938272; //adding proton mass to KE to get Energy in GeV
    //eN = eN + 0.939565; //adding neutron mass to KE to get Energy in GeV
    ePi0 = ePi0 + (0.134977)*(nPi0); //adding pion 0 mass to KE to get Energy in GeV
    ePip = ePip + (0.139571)*(nPip); //adding pion+ mass to KE to get Energy in Gev
    ePim = ePim + (0.139571)*(nPim); //adding pion- mass to KE to get Energy in Gev


    Percent_energy_Lep = True_LepE/Gen_numu_E;
    Percent_energy_N = eN/Gen_numu_E;
    Percent_energy_P = eP/Gen_numu_E;
    Percent_energy_Pip = ePip/Gen_numu_E;
    Percent_energy_Pim = ePim/Gen_numu_E;
    Percent_energy_Pi0 = ePi0/Gen_numu_E;
    Percent_energy_Other = eOther/Gen_numu_E;

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;


    //Energy Distribution Plots from Primary Neutrino Amongst Particles
    //0-suppression
    if (nLep !=0){plot_2D_Lep_zs->Fill(Gen_numu_E, Percent_energy_Lep);} //(x,y)
    if (nN != 0){plot_2D_N_zs->Fill(Gen_numu_E, Percent_energy_N);}
    if (nP != 0){plot_2D_P_zs->Fill(Gen_numu_E, Percent_energy_P);}
    if (nPip != 0){plot_2D_Pip_zs->Fill(Gen_numu_E, Percent_energy_Pip);}
    if (nPim != 0){plot_2D_Pim_zs->Fill(Gen_numu_E, Percent_energy_Pim);}
    if (nPi0 != 0){plot_2D_Pi0_zs->Fill(Gen_numu_E, Percent_energy_Pi0);}
    if (nOther !=0){plot_2D_Other_zs->Fill(Gen_numu_E, Percent_energy_Other);}

    if (nPi0 != 0){plot_2D_Pi0_EdepKE_zs->Fill(ePi0, Ratio_deposited_energy_Pi0);} //need pi0 to use all energy not only KE so moving above redefinition

    //Redefining Pion Energies to KE
    ePi0 = ePi0; // KE in GeV
    ePip = ePip; // KE in Gev
    ePim = ePim; // KE in Gev

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;

    //Energy Deposited in Detector Amongst Particles Across Primary Neutrino Energies Plots
    //0-suppression  
    if (nLep !=0){plot_2D_mu_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_mu);}
    if (nN != 0){plot_2D_N_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_N);}
    if (nP != 0){plot_2D_P_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_P);}
    if (nPip != 0){plot_2D_Pip_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_Pip);}
    if (nPim != 0){plot_2D_Pim_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_Pim);}
    if (nPi0 != 0){plot_2D_Pi0_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_Pi0);}
    if (nOther !=0){plot_2D_Other_Edep_zs->Fill(Gen_numu_E, Ratio_deposited_energy_Other);}

    //Energy Deposited in Detector Amongst Particles Across Particle Energies Plots
    //0-suppression 
    if (nLep !=0){plot_2D_mu_EdepKE_zs->Fill(Gen_numu_E, Ratio_deposited_energy_mu);}
    if (nN != 0){plot_2D_N_EdepKE_zs->Fill(eN, Ratio_deposited_energy_N);}
    if (nP != 0){plot_2D_P_EdepKE_zs->Fill(eP, Ratio_deposited_energy_P);}
    if (nPip != 0){plot_2D_Pip_EdepKE_zs->Fill(ePip, Ratio_deposited_energy_Pip);}
    if (nPim != 0){plot_2D_Pim_EdepKE_zs->Fill(ePim, Ratio_deposited_energy_Pim);}
    if (nOther !=0){plot_2D_Other_EdepKE_zs->Fill(eOther, Ratio_deposited_energy_Other);} 
    }//end CCNC_truth
    }// end ientry
    //Create Canvas,Draw Plot, Label Axes, Save
  
    //Kinetic Energy to Primary Neutrino Energy (Percent Energy) against Primary Neutrino Energy
    TCanvas *c1 = new TCanvas("c1", "c1",1800, 1350);
    c1->cd();
    plot_2D_Lep_zs->Draw("COLZ");
    plot_2D_Lep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Lep_zs->GetYaxis()->SetTitle("Lepton Total E/Primary #nu Energy "); //Set title of y-axis
    c1->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyLepton2D.jpg"); //Set file name
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 1350);
    c2->cd();
    plot_2D_N_zs->Draw("COLZ");
    plot_2D_N_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_zs->GetYaxis()->SetTitle("Neutron KE/Primary #nu Energy y"); //Set title of y-axis
    c2->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyNeutron2D.jpg"); //Set file name
    
    
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 1350);
    c3->cd();  
    plot_2D_P_zs->Draw("COLZ");
    plot_2D_P_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_zs->GetYaxis()->SetTitle("Proton KE/Primary #nu Energy "); //Set title of y-axis
    c3->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 1350);
    c4->cd();  
    plot_2D_Pip_zs->Draw("COLZ");
    plot_2D_Pip_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_zs->GetYaxis()->SetTitle("Pion+ total E/Primary #nu Energy "); //Set title of y-axis
    c4->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 1350);
    c5->cd();  
    plot_2D_Pim_zs->Draw("COLZ");
    plot_2D_Pim_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_zs->GetYaxis()->SetTitle("Pion- total E/Primary #nu Energy "); //Set title of y-axis
    c5->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 1350);
    c6->cd();  
    plot_2D_Pi0_zs->Draw("COLZ");
    plot_2D_Pi0_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_zs->GetYaxis()->SetTitle("Pion0 total E/Primary #nu Energy "); //Set title of y-axis
    c6->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c7 = new TCanvas("c7", "c7", 1800, 1350);
    c7->cd();  
    plot_2D_Other_zs->Draw("COLZ");
    plot_2D_Other_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_zs->GetYaxis()->SetTitle("Other Energy/Primary #nu Energy "); //Set title of y-axis
    c7->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/PercentEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Primary Neutrino Energy

    TCanvas *c8 = new TCanvas("c8", "c8", 1800, 1350);
    c8->cd();  
    plot_2D_mu_Edep_zs->Draw("COLZ");
    plot_2D_mu_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_mu_Edep_zs->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c8->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergymuon2D.jpg"); //Set file name

    TCanvas *c9 = new TCanvas("c9", "c9", 1800, 1350);
    c9->cd();
    plot_2D_N_Edep_zs->Draw("COLZ");
    plot_2D_N_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_Edep_zs->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c9->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyNeutron2D.jpg"); //Set file name
    
    TCanvas *c10 = new TCanvas("c10", "c10", 1800, 1350);
    c10->cd();  
    plot_2D_P_Edep_zs->Draw("COLZ");
    plot_2D_P_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_Edep_zs->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c10->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyProton2D.jpg"); //Set file name
    
    TCanvas *c11 = new TCanvas("c11", "c11", 1800, 1350);
    c11->cd();  
    plot_2D_Pip_Edep_zs->Draw("COLZ");
    plot_2D_Pip_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_Edep_zs->GetYaxis()->SetTitle("Pion+ Edep/Pi+ KE "); //Set title of y-axis
    c11->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyPion+2D.jpg"); //Set file name  
    
    TCanvas *c12 = new TCanvas("c12", "c12", 1800, 1350);
    c12->cd();  
    plot_2D_Pim_Edep_zs->Draw("COLZ");
    plot_2D_Pim_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_Edep_zs->GetYaxis()->SetTitle("Pion- Edep/Pi- KE "); //Set title of y-axis
    c12->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyPion-2D.jpg"); //Set file name  
    
    TCanvas *c13 = new TCanvas("c13", "c13", 1800, 1350);
    c13->cd();  
    plot_2D_Pi0_Edep_zs->Draw("COLZ");
    plot_2D_Pi0_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_Edep_zs->GetYaxis()->SetTitle("Pion0 Edep/Pi0 KE "); //Set title of y-axis
    c13->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyPion02D.jpg"); //Set file name  
    
    TCanvas *c14 = new TCanvas("c14", "c14", 1800, 1350);
    c14->cd();  
    plot_2D_Other_Edep_zs->Draw("COLZ");
    plot_2D_Other_Edep_zs->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_Edep_zs->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c14->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyOther2D.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Kinetic Energy

    TCanvas *c15 = new TCanvas("c15", "c15", 1800, 1350);
    c15->cd();  
    plot_2D_mu_EdepKE_zs->Draw("COLZ");
    plot_2D_mu_EdepKE_zs->GetXaxis()->SetTitle("Lepton Total E (GeV)"); //Set title of x-axis
    plot_2D_mu_EdepKE_zs->GetYaxis()->SetTitle("Lepton Edep/Lepton Total E "); //Set title of y-axis
    c15->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEmuon2D.jpg"); //Set file name

    TCanvas *c16 = new TCanvas("c16", "c16", 1800, 1350);
    c16->cd();
    plot_2D_N_EdepKE_zs->Draw("COLZ");
    plot_2D_N_EdepKE_zs->GetXaxis()->SetTitle("Neutron KE (GeV)"); //Set title of x-axis
    plot_2D_N_EdepKE_zs->GetYaxis()->SetTitle("Neutron Edep/Neutron KE "); //Set title of y-axis
    c16->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKENeutron2D.jpg"); //Set file name
    
    TCanvas *c17 = new TCanvas("c17", "c17", 1800, 1350);
    c17->cd();  
    plot_2D_P_EdepKE_zs->Draw("COLZ");
    plot_2D_P_EdepKE_zs->GetXaxis()->SetTitle("Proton KE (GeV)"); //Set title of x-axis
    plot_2D_P_EdepKE_zs->GetYaxis()->SetTitle("Proton Edep/Proton KE "); //Set title of y-axis
    c17->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEProton2D.jpg"); //Set file name
    
    TCanvas *c18 = new TCanvas("c18", "c18", 1800, 1350);
    c18->cd();  
    plot_2D_Pip_EdepKE_zs->Draw("COLZ");
    plot_2D_Pip_EdepKE_zs->GetXaxis()->SetTitle("Pi+ KE (GeV)"); //Set title of x-axis
    plot_2D_Pip_EdepKE_zs->GetYaxis()->SetTitle("Pi+ Edep/Pi+ KE "); //Set title of y-axis
    c18->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEPion+2D.jpg"); //Set file name  
    
    TCanvas *c19 = new TCanvas("c19", "c19", 1800, 1350);
    c19->cd();  
    plot_2D_Pim_EdepKE_zs->Draw("COLZ");
    plot_2D_Pim_EdepKE_zs->GetXaxis()->SetTitle("Pi- KE (GeV)"); //Set title of x-axis
    plot_2D_Pim_EdepKE_zs->GetYaxis()->SetTitle("Pi- Edep/Pi- KE "); //Set title of y-axis
    c19->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEPion-2D.jpg"); //Set file name  
    
    TCanvas *c20 = new TCanvas("c20", "c20", 1800, 1350);
    c20->cd();  
    plot_2D_Pi0_EdepKE_zs->Draw("COLZ");
    plot_2D_Pi0_EdepKE_zs->GetXaxis()->SetTitle("Pi0 total energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_EdepKE_zs->GetYaxis()->SetTitle("Pi0 Edep/Pi+ total energy "); //Set title of y-axis
    c20->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEPion02D.jpg"); //Set file name  
    
    TCanvas *c21 = new TCanvas("c21", "c21", 1800, 1350);
    c21->cd();  
    plot_2D_Other_EdepKE_zs->Draw("COLZ");
    plot_2D_Other_EdepKE_zs->GetXaxis()->SetTitle("Other KE (GeV)"); //Set title of x-axis
    plot_2D_Other_EdepKE_zs->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c21->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zerosuppressed/DepositedEnergyKEOther2D.jpg"); //Set file name
    

outFile.Close();
}
//Function will create sequence of plots with simulated data that is zero-suppressed and column normalized to 1.
void ZSCN(){
    cout << "plots will be column normalized and zero suppressed" << endl;

    //Defining Path to Branches
    //Kinetic Energy of Particles
    t1->SetBranchAddress("eN", &eN); //neutron
    t1->SetBranchAddress("eP", &eP); //proton
    t1->SetBranchAddress("ePi0", &ePi0); //pion0
    t1->SetBranchAddress("ePip", &ePip); //pion+
    t1->SetBranchAddress("ePim", &ePim); //pion-
    t1->SetBranchAddress("eOther", &eOther);
    t1->SetBranchAddress("True_LepE", &True_LepE); //lepton (muon)
    //Primary Neutrino Energy
    t1->SetBranchAddress("Gen_numu_E", &Gen_numu_E);
    //Charged Current Interaction
    t1->SetBranchAddress("CCNC_truth", &CCNC_truth);
    //Neutrino Scattering Mechanism
    t1->SetBranchAddress("Mode_truth", &Mode_truth);
    //Deposited Energy of Particles
    t1->SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2); //muon
    t1->SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2); //neutron
    t1->SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2); //proton
    t1->SetBranchAddress("Sim_pip_Edep_b2", &Sim_pip_Edep_b2); //pion+
    t1->SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2); //pion-
    t1->SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2); //pion0
    t1->SetBranchAddress("Sim_Other_Edep_b2", &Sim_Other_Edep_b2); 
    //number of particles per event
    t1->SetBranchAddress("nLep",&nLep);
    t1->SetBranchAddress("nN",&nN);
    t1->SetBranchAddress("nP",&nP);
    t1->SetBranchAddress("nPip", &nPip);
    t1->SetBranchAddress("nPim", &nPim);
    t1->SetBranchAddress("nPi0", &nPi0);
    t1->SetBranchAddress("nOther",&nOther);
    //for fiducial volume cut: Generator level outgoing lepton vtx
    t1->SetBranchAddress("Lepvtx_x",&Lepvtx_x);
    t1->SetBranchAddress("Lepvtx_y",&Lepvtx_y);
    t1->SetBranchAddress("Lepvtx_z",&Lepvtx_z);


    //2D plots for fraction of kinetic energy/neutrino energy over neutrino energy
    TH2F *plot_2D_Lep_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by #mu", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_N_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by Neutron", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_P_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by Proton", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pip_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by #pi^{+}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pim_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by #pi^{-}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Pi0_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by #pi^{0}", 50, 0, 7, 50, 0, 1);
    TH2F *plot_2D_Other_zscn = new TH2F("Statistics", "RP | Fraction of #nu Energy Received by Other", 50, 0, 7, 50, 0, 1);

    //2D plots for ratio of deposited energy to kinetic energy over neutrino energy
    TH2F *plot_2D_mu_Edep_zscn = new TH2F("Statistics", "RP | #mu Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_Edep_zscn = new TH2F("Statistics", "RP | Neutron Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_Edep_zscn = new TH2F("Statistics", "RP | Proton Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_Edep_zscn = new TH2F("Statistics", "RP | #pi^{+} Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_Edep_zscn = new TH2F("Statistics", "RP | #pi^{-} Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_Edep_zscn = new TH2F("Statistics", "RP | #pi^{0} Deposited Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_Edep_zscn = new TH2F("Statistics", "RP | Other Deposited Energy", 50, 0, 7, 50, 0, 1.4);

    //2D plots for ratio of deposited energy to kinetic energy over kinetic energy of particle
    TH2F *plot_2D_mu_EdepKE_zscn = new TH2F("Statistics", "RP | #mu Energy Deposited Across Total Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_N_EdepKE_zscn = new TH2F("Statistics", "RP | Neutron Energy Deposited Across KE", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_P_EdepKE_zscn = new TH2F("Statistics", "RP | Proton Energy Deposited Across KE", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pip_EdepKE_zscn = new TH2F("Statistics", "RP | #pi^{+} Energy Deposited Across KE", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pim_EdepKE_zscn = new TH2F("Statistics", "RP | #pi^{-} Energy Deposited Across KE", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Pi0_EdepKE_zscn = new TH2F("Statistics", "RP | #pi^{0} Energy Deposited Across Total Energy", 50, 0, 7, 50, 0, 1.4);
    TH2F *plot_2D_Other_EdepKE_zscn = new TH2F("Statistics", "RP | Other Energy Deposited Across KE", 50, 0, 7, 50, 0, 1.4);


    // Loop over all events
    int nentries = 0; // Total input events
    int pi0_count_Edep_over_one = 0; //counting number of pi0 with deposited energy larger than total energy
    int pi0_count_Edep_under_one = 0; //a check 
    int pi0_count_Edep_equals_one = 0; //a check
    int fiducial_volume_cut_counter = 0; //Counts how many events get cut out of plots
    int charged_current_counter = 0; //counts how many events are charged current 
    nentries = t1->GetEntries();
    cout<< "nentries: " << nentries<<endl;

    cout << "Zero-Suppressing plots..." << endl;

    cout << "Filter Neutrino Scattering Mode: N/A" << endl;

    for ( int ientry = 0; ientry < nentries; ientry++ )
    {
    t1->GetEntry(ientry);

    //fiducial volume cut
    if (abs(Lepvtx_x) > 310 ){  fiducial_volume_cut_counter++; continue;} //logic for continue : if an event falls within this condition, it is skipped
    if (abs(Lepvtx_y) > 550){ fiducial_volume_cut_counter++; continue;}
    if (50 > Lepvtx_z){  fiducial_volume_cut_counter++; continue;}
    if (Lepvtx_z > 1244){  fiducial_volume_cut_counter++; continue;}

    

    if (CCNC_truth == 0) { charged_current_counter++; // if 0 then charged current interaction, if 1 then neutral current interaction    
    //if (Mode_truth == 1) { // if 0=QE/El, 1=res, 2=DIS, 3=coherent production

    //Converting Deposited Energy Values from MeV to GeV
    Sim_mu_Edep_b2 = Sim_mu_Edep_b2/1000;
    Sim_n_Edep_b2 = Sim_n_Edep_b2/1000;
    Sim_p_Edep_b2 = Sim_p_Edep_b2/1000;
    Sim_pip_Edep_b2 = Sim_pip_Edep_b2/1000;
    Sim_pim_Edep_b2 = Sim_pim_Edep_b2/1000;
    Sim_pi0_Edep_b2 = Sim_pi0_Edep_b2/1000;
    Sim_Other_Edep_b2 = Sim_Other_Edep_b2/1000;    

    //Convert KE to total E: multiply masses by number of particle and add to KE
    //eP = eP + 0.938272; //adding proton mass to KE to get Energy in GeV
    //eN = eN + 0.939565; //adding neutron mass to KE to get Energy in GeV
    ePi0 = ePi0 + (0.134977)*(nPi0); //adding pion 0 mass to KE to get total Energy in GeV
    ePip = ePip + (0.139571)*(nPip); //adding pion+ mass to KE to get total Energy in Gev
    ePim = ePim + (0.139571)*(nPim); //adding pion- mass to KE to get total Energy in Gev



    Percent_energy_Lep = True_LepE/Gen_numu_E; //total over neutrino E
    Percent_energy_N = eN/Gen_numu_E; //KE over neutrino E
    Percent_energy_P = eP/Gen_numu_E; //KE over neutrino E
    Percent_energy_Pip = ePip/Gen_numu_E; // total E over neutrino E
    Percent_energy_Pim = ePim/Gen_numu_E; // total E over neutrino E
    Percent_energy_Pi0 = ePi0/Gen_numu_E; // total E over neutrino E
    Percent_energy_Other = eOther/Gen_numu_E; //KE over neutrino E

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE; //Edep over total E
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN; //Edep over KE
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP; //Edep over KE
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip; //Edep over total E
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim; //Edep over total E
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0; //Edep over total E
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther; //Edep over KE


    //Energy Distribution Plots from Primary Neutrino Amongst Particles
    //0-supression
    if (nLep !=0){plot_2D_Lep_zscn->Fill(Gen_numu_E, Percent_energy_Lep);} // (x,y) | y-axis: total E over neutrino E
    if (nN != 0){plot_2D_N_zscn->Fill(Gen_numu_E, Percent_energy_N);} // y-axis: KE over neutrino E
    if (nP != 0){plot_2D_P_zscn->Fill(Gen_numu_E, Percent_energy_P);} // y-axis: KE over neutrino E
    if (nPip != 0){plot_2D_Pip_zscn->Fill(Gen_numu_E, Percent_energy_Pip);} // y-axis: total E over neutrino E
    if (nPim != 0){plot_2D_Pim_zscn->Fill(Gen_numu_E, Percent_energy_Pim);} // y-axis: total E over neutrino E
    if (nPi0 != 0){plot_2D_Pi0_zscn->Fill(Gen_numu_E, Percent_energy_Pi0);} // y-axis: total E over neutrino E
    if (nOther !=0){plot_2D_Other_zscn->Fill(Gen_numu_E, Percent_energy_Other);} // y-axis: KE over neutrino E

    //Energy Deposited in Detector Amongst Particles Across Particle Energies Plots
    if (nPi0 != 0){plot_2D_Pi0_EdepKE_zscn->Fill(ePi0, Ratio_deposited_energy_Pi0);} //need pi0 to use all energy not only KE so moving above redefinition
    
    if (Sim_pi0_Edep_b2 > ePi0) {
        pi0_count_Edep_over_one++;
    }

    if (Sim_pi0_Edep_b2 < ePi0) {
        pi0_count_Edep_under_one++;
    }

    if (Sim_pi0_Edep_b2 == ePi0) {
        pi0_count_Edep_equals_one ++;
    }   

    //Redefining Pion Total Energies to KE
    ePi0 = ePi0; // KE in Gev
    ePip = ePip; // KE in Gev
    ePim = ePim; // KE in Gev

    Ratio_deposited_energy_mu = Sim_mu_Edep_b2/True_LepE;
    Ratio_deposited_energy_N = Sim_n_Edep_b2/eN;
    Ratio_deposited_energy_P = Sim_p_Edep_b2/eP;
    Ratio_deposited_energy_Pip = Sim_pip_Edep_b2/ePip;
    Ratio_deposited_energy_Pim = Sim_pim_Edep_b2/ePim;
    Ratio_deposited_energy_Pi0 = Sim_pi0_Edep_b2/ePi0;
    Ratio_deposited_energy_Other = Sim_Other_Edep_b2/eOther;



    //Energy Deposited in Detector Amongst Particles Across Primary Neutrino Energies Plots
      //0-supression  
    if (nLep !=0){plot_2D_mu_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_mu);}
    if (nN != 0){plot_2D_N_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_N);}
    if (nP != 0){plot_2D_P_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_P);}
    if (nPip != 0){plot_2D_Pip_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_Pip);}
    if (nPim != 0){plot_2D_Pim_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_Pim);}
    if (nPi0 != 0){plot_2D_Pi0_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_Pi0);}
    if (nOther !=0){plot_2D_Other_Edep_zscn->Fill(Gen_numu_E, Ratio_deposited_energy_Other);}

    //Energy Deposited in Detector Amongst Particles Across Particle Energies Plots
    //0-supression 
    if (nLep !=0){plot_2D_mu_EdepKE_zscn->Fill(True_LepE, Ratio_deposited_energy_mu);}
    if (nN != 0){plot_2D_N_EdepKE_zscn->Fill(eN, Ratio_deposited_energy_N);}
    if (nP != 0){plot_2D_P_EdepKE_zscn->Fill(eP, Ratio_deposited_energy_P);}
    if (nPip != 0){plot_2D_Pip_EdepKE_zscn->Fill(ePip, Ratio_deposited_energy_Pip);}
    if (nPim != 0){plot_2D_Pim_EdepKE_zscn->Fill(ePim, Ratio_deposited_energy_Pim);}
    if (nOther !=0){plot_2D_Other_EdepKE_zscn->Fill(eOther, Ratio_deposited_energy_Other);}
    //}//end Mode_truth
    }//end CCNC_truth
    }// end ientry
  
    cout << "pi0 with Edep larger than total E: " << pi0_count_Edep_over_one << endl;
    cout << "pi0 with Edep less than total E: " << pi0_count_Edep_under_one << endl;
    cout << "pi0 with Edep equal to total E: " << pi0_count_Edep_equals_one << endl;
    cout << "number of events cut due to outside fiducial volume: " << fiducial_volume_cut_counter << endl;
    cout << "number of charged current events: " << charged_current_counter << endl;
    
    int total_pi0 = 0;
    total_pi0 = pi0_count_Edep_over_one + pi0_count_Edep_under_one + pi0_count_Edep_equals_one;

    cout << "total_pi0: " << total_pi0 << endl;

    //Column Normalizing Plots to 1 
    cout << "Column Normalizing to 1..." << endl;
    //Percent Energy plots
    for (int ix=1; ix<=plot_2D_Lep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Lep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Lep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Lep_zscn->GetNbinsY(); iy++){
        plot_2D_Lep_zscn->SetBinContent(ix, iy, plot_2D_Lep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_zscn->GetNbinsY(); iy++){
        plot_2D_N_zscn->SetBinContent(ix, iy, plot_2D_N_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_zscn->GetNbinsY(); iy++){
        plot_2D_P_zscn->SetBinContent(ix, iy, plot_2D_P_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_zscn->GetNbinsY(); iy++){
        plot_2D_Pip_zscn->SetBinContent(ix, iy, plot_2D_Pip_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_zscn->GetNbinsY(); iy++){
        plot_2D_Pim_zscn->SetBinContent(ix, iy, plot_2D_Pim_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_zscn->GetNbinsY(); iy++){
        plot_2D_Pi0_zscn->SetBinContent(ix, iy, plot_2D_Pi0_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_zscn->GetNbinsY(); iy++){
        plot_2D_Other_zscn->SetBinContent(ix, iy, plot_2D_Other_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }
    //Deposited Energy Plots
    for (int ix=1; ix<=plot_2D_mu_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_mu_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_mu_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_mu_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_mu_Edep_zscn->SetBinContent(ix, iy, plot_2D_mu_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_N_Edep_zscn->SetBinContent(ix, iy, plot_2D_N_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_P_Edep_zscn->SetBinContent(ix, iy, plot_2D_P_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_Pip_Edep_zscn->SetBinContent(ix, iy, plot_2D_Pip_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_Pim_Edep_zscn->SetBinContent(ix, iy, plot_2D_Pim_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_Pi0_Edep_zscn->SetBinContent(ix, iy, plot_2D_Pi0_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_Edep_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_Edep_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_Edep_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_Edep_zscn->GetNbinsY(); iy++){
        plot_2D_Other_Edep_zscn->SetBinContent(ix, iy, plot_2D_Other_Edep_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }  
    //Deposited Energy/Kinetic Energy Plots
    for (int ix=1; ix<=plot_2D_mu_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_mu_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_mu_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_mu_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_mu_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_mu_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_N_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_N_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_N_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_N_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_N_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_N_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_P_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_P_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_P_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_P_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_P_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_P_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pip_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pip_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pip_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pip_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_Pip_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_Pip_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pim_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pim_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pim_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pim_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_Pim_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_Pim_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Pi0_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Pi0_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Pi0_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Pi0_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_Pi0_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_Pi0_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    for (int ix=1; ix<=plot_2D_Other_EdepKE_zscn->GetNbinsX(); ix++){
    double colnorm =0.;
    for (int iy=1; iy<=plot_2D_Other_EdepKE_zscn->GetNbinsY(); iy++){
        colnorm += plot_2D_Other_EdepKE_zscn->GetBinContent(ix,iy);
    }
    for (int iy=1; iy<= plot_2D_Other_EdepKE_zscn->GetNbinsY(); iy++){
        plot_2D_Other_EdepKE_zscn->SetBinContent(ix, iy, plot_2D_Other_EdepKE_zscn->GetBinContent(ix, iy)/colnorm);
    }
    }

    //Create Canvas,Draw Plot, Label Axes, Save
   //Kinetic Energy to Primary Neutrino Energy (Percent Energy) against Primary Neutrino Energy
    TCanvas *c1 = new TCanvas("c1", "c1",1800, 1350);
    c1->cd();
    gStyle->SetOptStat(11);
    plot_2D_Lep_zscn->Draw("COLZ"); //Sets z-axis
    //plot_2D_Lep_zscn->GetZaxis()->SetRangeUser(0,1.0); //Sets z-axis range
    plot_2D_Lep_zscn->SetTitleSize(0.06);
    plot_2D_Lep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Lep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Lep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Lep_zscn->GetYaxis()->SetTitle("Muon Total E/Primary #nu Energy "); //Set title of y-axis
    c1->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyLepton.jpg"); //Set file name
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 1350);
    c2->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_zscn->Draw("COLZ"); 
    plot_2D_N_zscn->SetTitleSize(0.06);
    plot_2D_N_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_N_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_N_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_zscn->GetYaxis()->SetTitle("Neutrons KE/Primary #nu Energy "); //Set title of y-axis
    c2->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyNeutron2D_RP.jpg"); //Set file name
    
    
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 1350);
    c3->cd();
    gStyle->SetOptStat(11);  
    plot_2D_P_zscn->Draw("COLZ");
    plot_2D_P_zscn->SetTitleSize(0.06);
    plot_2D_P_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_P_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_P_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_zscn->GetYaxis()->SetTitle("Protons KE/Primary #nu Energy "); //Set title of y-axis
    c3->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyProton2D_RP.jpg"); //Set file name
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 1350);
    c4->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_zscn->Draw("COLZ");
    plot_2D_Pip_zscn->SetTitleSize(0.06);
    plot_2D_Pip_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pip_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pip_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_zscn->GetYaxis()->SetTitle("Pions+ Total E/Primary #nu Energy "); //Set title of y-axis
    c4->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyPion+2D_RP.jpg"); //Set file name  
    
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 1350);
    c5->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pim_zscn->Draw("COLZ");
    plot_2D_Pim_zscn->SetTitleSize(0.06);
    plot_2D_Pim_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pim_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pim_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_zscn->GetYaxis()->SetTitle("Pions- Total E/Primary #nu Energy "); //Set title of y-axis
    c5->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyPion-2D_RP.jpg"); //Set file name  
    
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 1350);
    c6->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pi0_zscn->Draw("COLZ");
    plot_2D_Pi0_zscn->SetTitleSize(0.06);
    plot_2D_Pi0_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_zscn->GetYaxis()->SetTitle("Pions0 Total E/Primary #nu Energy "); //Set title of y-axis
    c6->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyPion02D_RP.jpg"); //Set file name  
    
    TCanvas *c7 = new TCanvas("c7", "c7", 1800, 1350);
    c7->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_zscn->Draw("COLZ");
    plot_2D_Other_zscn->SetTitleSize(0.06);
    plot_2D_Other_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Other_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Other_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_zscn->GetYaxis()->SetTitle("Other Energy/Primary #nu Energy "); //Set title of y-axis
    c7->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/PercentEnergyOther2D_RP.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Primary Neutrino Energy

    TCanvas *c8 = new TCanvas("c8", "c8", 1800, 1350);
    c8->cd();  
    gStyle->SetOptStat(11);
    plot_2D_mu_Edep_zscn->Draw("COLZ");
    plot_2D_mu_Edep_zscn->SetTitleSize(0.06);
    plot_2D_mu_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_mu_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_mu_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_mu_Edep_zscn->GetYaxis()->SetTitle("Muons Edep/Muons Total E "); //Set title of y-axis
    c8->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergymuon2D_RP.jpg"); //Set file name

    TCanvas *c9 = new TCanvas("c9", "c9", 1800, 1350);
    c9->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_Edep_zscn->Draw("COLZ");
    plot_2D_N_Edep_zscn->SetTitleSize(0.06);
    plot_2D_N_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_N_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_N_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_N_Edep_zscn->GetYaxis()->SetTitle("Neutrons Edep/Neutrons KE "); //Set title of y-axis
    c9->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyNeutron2D_RP.jpg"); //Set file name
    
    TCanvas *c10 = new TCanvas("c10", "c10", 1800, 1350);
    c10->cd();  
    gStyle->SetOptStat(11);
    plot_2D_P_Edep_zscn->Draw("COLZ");
    plot_2D_P_Edep_zscn->SetTitleSize(0.06);
    plot_2D_P_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_P_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_P_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_P_Edep_zscn->GetYaxis()->SetTitle("Protons Edep/Protons KE "); //Set title of y-axis
    c10->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyProton2D_RP.jpg"); //Set file name
    
    TCanvas *c11 = new TCanvas("c11", "c11", 1800, 1350);
    c11->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_Edep_zscn->Draw("COLZ");
    plot_2D_Pip_Edep_zscn->SetTitleSize(0.06);
    plot_2D_Pip_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pip_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pip_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pip_Edep_zscn->GetYaxis()->SetTitle("Pions+ Edep/Pions+ KE "); //Set title of y-axis
    c11->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyPion+2D_RP.jpg"); //Set file name  
    
    TCanvas *c12 = new TCanvas("c12", "c12", 1800, 1350);
    c12->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pim_Edep_zscn->Draw("COLZ");
    plot_2D_Pim_Edep_zscn->SetTitleSize(0.06);
    plot_2D_Pim_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pim_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pim_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pim_Edep_zscn->GetYaxis()->SetTitle("Pions- Edep/Pions- KE "); //Set title of y-axis
    c12->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyPion-2D_RP.jpg"); //Set file name  
    
    TCanvas *c13 = new TCanvas("c13", "c13", 1800, 1350);
    c13->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pi0_Edep_zscn->Draw("COLZ");
    plot_2D_Pi0_Edep_zscn->SetTitleSize(0.06);
    plot_2D_Pi0_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_Edep_zscn->GetYaxis()->SetTitle("Pions0 Edep/Pions0 KE "); //Set title of y-axis
    c13->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyPion02D_RP.jpg"); //Set file name  
    
    TCanvas *c14 = new TCanvas("c14", "c14", 1800, 1350);
    c14->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_Edep_zscn->Draw("COLZ");
    plot_2D_Other_Edep_zscn->SetTitleSize(0.06);
    plot_2D_Other_Edep_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Other_Edep_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Other_Edep_zscn->GetXaxis()->SetTitle("Primary Neutrino Energy (GeV)"); //Set title of x-axis
    plot_2D_Other_Edep_zscn->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c14->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyOther2D_RP.jpg"); //Set file name

    //Ratio Deposited Energy to Kinetic Energy against Kinetic Energy

    TCanvas *c15 = new TCanvas("c15", "c15", 1800, 1350);
    c15->cd();  
    gStyle->SetOptStat(11);
    plot_2D_mu_EdepKE_zscn->Draw("COLZ");
    plot_2D_mu_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_mu_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_mu_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_mu_EdepKE_zscn->GetXaxis()->SetTitle("Muon Total E (GeV)"); //Set title of x-axis
    plot_2D_mu_EdepKE_zscn->GetYaxis()->SetTitle("Muon Edep/Muon Total E "); //Set title of y-axis
    c15->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEmuon2D_RP.jpg"); //Set file name

    TCanvas *c16 = new TCanvas("c16", "c16", 1800, 1350);
    c16->cd();
    gStyle->SetOptStat(11);
    plot_2D_N_EdepKE_zscn->Draw("COLZ");
    plot_2D_N_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_N_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_N_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_N_EdepKE_zscn->GetXaxis()->SetTitle("Neutrons KE (GeV)"); //Set title of x-axis
    plot_2D_N_EdepKE_zscn->GetYaxis()->SetTitle("Neutrons Edep/Neutron KE "); //Set title of y-axis
    c16->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKENeutron2D_RP.jpg"); //Set file name
    
    TCanvas *c17 = new TCanvas("c17", "c17", 1800, 1350);
    c17->cd();  
    gStyle->SetOptStat(11);
    plot_2D_P_EdepKE_zscn->Draw("COLZ");
    plot_2D_P_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_P_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_P_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_P_EdepKE_zscn->GetXaxis()->SetTitle("Protons KE (GeV)"); //Set title of x-axis
    plot_2D_P_EdepKE_zscn->GetYaxis()->SetTitle("Protons Edep/Proton KE "); //Set title of y-axis
    c17->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEProton2D_RP.jpg"); //Set file name
    
    TCanvas *c18 = new TCanvas("c18", "c18", 1800, 1350);
    c18->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pip_EdepKE_zscn->Draw("COLZ");
    plot_2D_Pip_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_Pip_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pip_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pip_EdepKE_zscn->GetXaxis()->SetTitle("Pi+ KE (GeV)"); //Set title of x-axis
    plot_2D_Pip_EdepKE_zscn->GetYaxis()->SetTitle("Pions+ Edep/Pions+ KE "); //Set title of y-axis
    c18->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEPion+2D_RP.jpg"); //Set file name  
    
    TCanvas *c19 = new TCanvas("c19", "c19", 1800, 1350);
    c19->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pim_EdepKE_zscn->Draw("COLZ");
    plot_2D_Pim_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_Pim_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pim_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pim_EdepKE_zscn->GetXaxis()->SetTitle("Pi- KE (GeV)"); //Set title of x-axis
    plot_2D_Pim_EdepKE_zscn->GetYaxis()->SetTitle("Pions- Edep/Pions- KE "); //Set title of y-axis
    c19->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEPion-2D_RP.jpg"); //Set file name  
    
    TCanvas *c20 = new TCanvas("c20", "c20", 1800, 1350);
    c20->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Pi0_EdepKE_zscn->Draw("COLZ");
    plot_2D_Pi0_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_Pi0_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Pi0_EdepKE_zscn->GetXaxis()->SetTitle("Pi0 total energy (GeV)"); //Set title of x-axis
    plot_2D_Pi0_EdepKE_zscn->GetYaxis()->SetTitle("Pions0 Edep/Pions0 total energy "); //Set title of y-axis
    c20->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEPion02D_RP.jpg"); //Set file name  
    
    TCanvas *c21 = new TCanvas("c21", "c21", 1800, 1350);
    c21->cd();  
    gStyle->SetOptStat(11);
    plot_2D_Other_EdepKE_zscn->Draw("COLZ");
    plot_2D_Other_EdepKE_zscn->SetTitleSize(0.06);
    plot_2D_Other_EdepKE_zscn->GetXaxis()->SetTitleSize(0.05);
    plot_2D_Other_EdepKE_zscn->GetYaxis()->SetTitleSize(0.05);
    plot_2D_Other_EdepKE_zscn->GetXaxis()->SetTitle("Other KE (GeV)"); //Set title of x-axis
    plot_2D_Other_EdepKE_zscn->GetYaxis()->SetTitle("Other Edep/Other KE "); //Set title of y-axis
    c21->SaveAs("/exp/dune/app/users/mfucci/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/NNBurrito/fnal_data/zscn/DepositedEnergyKEOther2D_RP.jpg"); //Set file name
    

outFile.Close();
}

//Main Function
//Calls on specific functions as defined above based on the value of plot_filter
//Set plot_filter to "Raw", "Column_Normalized", "Zero_Suppressed", or " ZSCN" (Zero Suppressed + Column Normalized) to produce the respective plots
void NeutrinoEnergyDistributionPlotsMaker_FNAL_DataSource(TString plot_filter = "ZSCN") {

if (plot_filter == "Raw"){
    Raw();
}

 else if (plot_filter == "Column_Normalized"){
    Column_Normalized();
}

 else if (plot_filter == "Zero_Suppressed"){
    Zero_Suppressed();
}

else if (plot_filter == "ZSCN"){
    ZSCN();
}

else {
    cout << "Set plot_filter to an acceptable string: Raw, Column_Normalized, Zero_Suppressed, or ZSCN." << endl;

}
}

