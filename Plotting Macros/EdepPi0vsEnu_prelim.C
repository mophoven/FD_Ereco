void EdepPi0vsEnu_prelim(){
TFile* rootFile = TFile::Open("test5files.root");
        if (!rootFile){
        std::cerr << "Error: Could not open file" << std::endl;
        return;
        }
        //rootFile -> cd ("MyEnergyAnalysis");
TTree* ntuple = dynamic_cast<TTree*>(rootFile->Get("MyEnergyAnalysis/MyTree"));
        if (!ntuple){
        std::cerr << "Error: Could not find nTuple" << std::endl;
        return;
        }
TH2F* PeNu = new TH2F("EdepPi0_prelim", "Energy Deposited by pi0 vs True Neutrino Energy", 40, 0, 7, 40, 0, 3);


        double Gen_numu_EVal, Sim_pi0_Edep_b2Val, ePi0Val, dep_E_ratioVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2Val);
        ntuple -> SetBranchAddress("ePi0", &ePi0Val);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);
                        if(ePi0Val==0 || Sim_pi0_Edep_b2Val == 0) {
                                continue;
                        }
                dep_E_ratioVal = (Sim_pi0_Edep_b2Val/1000) / ePi0Val;

                std::cout << dep_E_ratioVal << "   " << Gen_numu_EVal << std::endl;

                PeNu->Fill(Gen_numu_EVal, dep_E_ratioVal);
        }
TCanvas* canvas = new TCanvas();

PeNu->Draw("colz");

PeNu->SetTitle("Energy Deposited by Pi0 vs True Neutrino Energy");
PeNu->GetXaxis()->SetTitle("Generated Neutrino Energy (GeV)");
PeNu->GetYaxis()->SetTitle("Deposted pi0 Energy / True Pi0 Energy");

canvas->Draw("colz");
}
