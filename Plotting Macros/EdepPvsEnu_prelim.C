void EdepPvsEnu_prelim(){
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
TH2F* P = new TH2F("EdepP_prelim", "Energy Deposited by Protons vs True Neutrino Energy", 40, 0, 7, 40, 0, 3);


        double Gen_numu_EVal, Sim_p_Edep_b2Val, ePVal, dep_E_ratioVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2Val);
        ntuple -> SetBranchAddress("eP", &ePVal);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);
                        if(ePVal==0 || Sim_p_Edep_b2Val == 0) {
                                continue;
                        }
                dep_E_ratioVal = (Sim_p_Edep_b2Val/1000) / ePVal;

                std::cout << dep_E_ratioVal << "   " << Gen_numu_EVal << std::endl;

                P->Fill(Gen_numu_EVal, dep_E_ratioVal);
        }
TCanvas* canvas = new TCanvas();

P->Draw("colz");

P->SetTitle("Energy Deposited by Protons vs True Neutrino Energy");
P->GetXaxis()->SetTitle("Generated Neutrino Energy (GeV)");
P->GetYaxis()->SetTitle("Deposted Proton Energy / True Proton Energy");

canvas->Draw("colz");
}
