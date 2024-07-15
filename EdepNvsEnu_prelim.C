void EdepNvsEnu_prelim(){
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
TH2F* N = new TH2F("EdepN_prelim", "Energy Deposited by Neutrons vs True Neutrino Energy", 40, 0, 7, 40, 0, 3);


        double Gen_numu_EVal, Sim_n_Edep_b2Val, eNVal, dep_E_ratioVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2Val);
        ntuple -> SetBranchAddress("eN", &eNVal);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);
                        if(eNVal==0 || Sim_n_Edep_b2Val == 0) {
                                continue;
                        }
                dep_E_ratioVal = (Sim_n_Edep_b2Val/1000) / eNVal;

                std::cout << dep_E_ratioVal << "   " << Gen_numu_EVal << std::endl;

                N->Fill(Gen_numu_EVal, dep_E_ratioVal);
        }
TCanvas* canvas = new TCanvas();

N->Draw("colz");

N->SetTitle("Energy Deposited by Neutrons vs True Neutrino Energy");
N->GetXaxis()->SetTitle("Generated Neutrino Energy (GeV)");
N->GetYaxis()->SetTitle("Deposted Neutron Energy / True Neutron Energy");

canvas->Draw("colz");
}
