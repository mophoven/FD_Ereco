void EdepmuvsEnu_prelim(){
TFile* rootFile = TFile::Open("myntuple.root");
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
TH2F* PeNu = new TH2F("Edepmu", "Energy Deposited by muons vs True Neutrino Energy", 50, 0, 7, 50, 0, 10);


        double Gen_numu_EVal, Sim_mu_Edep_b2Val, Sim_EMuVal, dep_E_ratioVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2Val);
        //ntuple -> SetBranchAddress("Sim_nMu", &Sim_nMuVal);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);

          // Section to add cuts, here is an example of one to reject points if either value in the ratio is zero.
          
          /*if(Sim_mu_Edep_b2Val != 0 && Sim_nMuVal !=0) {
                                dep_N_ratioVal = Sim_mu_Edep_b2Val / Sim_nMuVal;
                        } else {
                                dep_N_ratioVal = 0;
                        }
                        */

                PeNu->Fill(Gen_numu_EVal, Sim_mu_Edep_b2Val);
        }
TCanvas* canvas = new TCanvas();

PeNu->Draw("colz");

PeNu->SetTitle("Energy Deposited by muons vs True Neutrino Energy");
PeNu->GetXaxis()->SetTitle("Beam Nu_e Energy (GeV)");
PeNu->GetYaxis()->SetTitle("Deposted mu Energy / mu Events (MeV)");

canvas->Draw("colz");

}
