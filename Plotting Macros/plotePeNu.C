void plotePeNu(){
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
TH2F* PeNu = new TH2F("PeNu", "Proton Energy vs. Generator Nu mu Energy", 50, 0, 50, 50, 0, 10);

        double Gen_numu_EVal, ePVal;
        ntuple -> SetBranchAddress("Gen_numu_E", &Gen_numu_EVal);
        ntuple -> SetBranchAddress("eP", &ePVal);

        for (Long64_t i = 0; i < ntuple->GetEntries(); ++i){
                ntuple->GetEntry(i);

                PeNu->Fill(Gen_numu_EVal, ePVal);
        }
TCanvas* canvas = new TCanvas();

PeNu->Draw("colz");

PeNu->SetTitle("Proton Energy vs. Generator Nu mu Energy");
PeNu->GetXaxis()->SetTitle("Generated Nu mu Energy (GeV)");
PeNu->GetYaxis()->SetTitle("Proton Energy (GeV)");

canvas->Draw("colz");

}
