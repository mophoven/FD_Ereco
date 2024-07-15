void EdepPimvsEnu_prelim(){
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
TH2F* Pim = new TH2F("EdepPi-_prelim", "Energy Deposited by pi- vs True Neutrino Energy", 40, 0, 7, 40, 0, 3);


        double Gen_numu_EVal, Sim_pim_Edep_b2Val, ePimVal, dep_E_ratioVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2Val);
        ntuple -> SetBranchAddress("ePim", &ePimVal);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);
                        if(ePimVal==0 || Sim_pim_Edep_b2Val == 0) {
                                continue;
                        }
                dep_E_ratioVal = (Sim_pim_Edep_b2Val/1000) / ePimVal;

                std::cout << dep_E_ratioVal << "   " << Gen_numu_EVal << std::endl;

                Pim->Fill(Gen_numu_EVal, dep_E_ratioVal);
        }
TCanvas* canvas = new TCanvas();

Pim->Draw("colz");

Pim->SetTitle("Energy Deposited by Pi- vs True Neutrino Energy");
Pim->GetXaxis()->SetTitle("Generated Neutrino Energy (GeV)");
Pim->GetYaxis()->SetTitle("Deposted pi- Energy / True pi- Energy");

canvas->Draw("colz");
}
