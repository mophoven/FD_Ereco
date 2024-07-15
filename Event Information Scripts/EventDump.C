void EventDump(){
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

//Vec Vectors;
//ntuple->MakeClass("Vec", &Vectors);

int Event, Run, SubRun;

double eP, eN, ePip, ePim, ePi0, eOther;

int nLep, nP, nN, nPip, nPim, nPi0, nOther;

int Sim_nEle, Sim_nNue, Sim_nMu, Sim_nNuMu, Sim_nTau, Sim_nNutau, Sim_nPhoton, Sim_nPionNeutral, Sim_nPionCharged, Sim_nNeutron, Sim_nProton;

int CCNC_truth, Mode_truth;

double Gen_numu_E, Sim_numu_E, Sim_mu_Edep_b2, Sim_n_Edep_b2, Sim_p_Edep_b2, Sim_pip_Edep_b2, Sim_pim_Edep_b2, Sim_pi0_Edep_b2, Sim_Other_Edep_b2;

 //std::vector<int>* P_PDG;
 //std::vector<int>* P_mother;
 //std::vector<float> P_vtx_x;
 //std::vector<float> P_vtx_y;
 //std::vector<float> P_vtx_z;
 //std::vector<float> P_E;




        ntuple->SetBranchAddress("Event", &Event);
        ntuple->SetBranchAddress("Run", &Run);
        ntuple->SetBranchAddress("SubRun", &SubRun);
        ntuple->SetBranchAddress("eP", &eP);
        ntuple->SetBranchAddress("eN", &eN);
        ntuple->SetBranchAddress("ePip", &ePip);
        ntuple->SetBranchAddress("ePim", &ePim);
        ntuple->SetBranchAddress("ePi0", &ePi0);
        ntuple->SetBranchAddress("eOther", &eOther);
        ntuple->SetBranchAddress("nLep", &nLep);
        ntuple->SetBranchAddress("nP", &nP);
        ntuple->SetBranchAddress("nN", &nN);
        ntuple->SetBranchAddress("nPip", &nPip);
        ntuple->SetBranchAddress("nPim", &nPim);
        ntuple->SetBranchAddress("nPi0", &nPi0);
        ntuple->SetBranchAddress("nOther", &nOther);
        ntuple->SetBranchAddress("Sim_nEle", &Sim_nEle);
        ntuple->SetBranchAddress("Sim_nNue", &Sim_nNue);
        ntuple->SetBranchAddress("Sim_nMu", &Sim_nMu);
        ntuple->SetBranchAddress("Sim_nNumu", &Sim_nNuMu);
        ntuple->SetBranchAddress("Sim_nTau", &Sim_nTau);
        ntuple->SetBranchAddress("Sim_nNutau", &Sim_nNutau);
        ntuple->SetBranchAddress("Sim_nPhoton", &Sim_nPhoton);
        ntuple->SetBranchAddress("Sim_nPionNeutral", &Sim_nPionNeutral);
        ntuple->SetBranchAddress("Sim_nPionCharged", &Sim_nPionCharged);
        ntuple->SetBranchAddress("Sim_nNeutron", &Sim_nNeutron);
        ntuple->SetBranchAddress("Sim_nProton", &Sim_nProton);
        ntuple->SetBranchAddress("CCNC_truth", &CCNC_truth);
        ntuple->SetBranchAddress("Mode_truth", &Mode_truth);
        ntuple->SetBranchAddress("Gen_numu_E", &Gen_numu_E);
        ntuple->SetBranchAddress("Sim_numu_E", &Sim_numu_E);
        ntuple->SetBranchAddress("Sim_mu_Edep_b2", &Sim_mu_Edep_b2);
        ntuple->SetBranchAddress("Sim_n_Edep_b2", &Sim_n_Edep_b2);
        ntuple->SetBranchAddress("Sim_p_Edep_b2", &Sim_p_Edep_b2);
        ntuple->SetBranchAddress("Sim_pip_Edep_b2", &Sim_pip_Edep_b2);
        ntuple->SetBranchAddress("Sim_pim_Edep_b2", &Sim_pim_Edep_b2);
        ntuple->SetBranchAddress("Sim_pi0_Edep_b2", &Sim_pi0_Edep_b2);
        ntuple->SetBranchAddress("Sim_Other_Edep_b2", &Sim_Other_Edep_b2);
        //ntuple->SetBranchAddress("P_E", &P_E);
        //ntuple->SetBranchAddress("P_PDG", &Vec);
        //ntuple->SetBranchAddress("P_mother", &Vector);
        //ntuple->SetBranchAddress("P_vtx_x", &P_vtx_x);
        //ntuple->SetBranchAddress("P_vtx_y", &P_vtx_y);
        //ntuple->SetBranchAddress("P_vtx_z", &P_vtx_z);



for (Long64_t i = 0; i < 500; ++i){
                ntuple->GetEntry(i);
                        std::cout << i << std::endl;
                        std::cout << "Event Number: " << Event << std::endl;
                        std::cout << "Run: " << Run << std::endl;
                        std::cout << "SubRun: " << SubRun << std::endl;
                        std::cout << "Neutral (1) or Charged (0) Current Interaction? " << CCNC_truth << std::endl;
                        std::cout << "Interaction Type (0=QE, 1=RP, 2=DIS, 3=CP): " << Mode_truth << std::endl;
                        std::cout << "Generated Neutrino Energy: " << Gen_numu_E << std::endl;
                        std::cout << "Simulated Nuetrino Energy: " << Sim_numu_E << std::endl;
                        std::cout << "Number of Protons: " << nP << std::endl;
                        std::cout << "True Proton Energy: " << eP << std::endl;
                        std::cout << "Number of Neutrons: " << nN << std::endl;
                        std::cout << "True Neutron Energy: " << eN << std::endl;
                        std::cout << "Number of Pi+: " << nPip << std::endl;
                        std::cout << "True Pi+ Energy: " << ePip << std::endl;
                        std::cout << "Number of Pi-: " << nPim << std::endl;
                        std::cout << "True Pi- Energy: " << ePim << std::endl;
                        std::cout << "Number of Pi0: " << nPi0 << std::endl;
                        std::cout << "True Pi0 Energy: " << ePi0 << std::endl;
                        std::cout << "Number of Other Particles: " << nOther << std::endl;
                        std::cout << "True Energy of Other Particles: " << eOther <<std::endl;
                        std::cout << "Number of Simulated Electrons: " << Sim_nEle << std::endl;
                        std::cout << "Number of Simulated Electron Neutrinos: " << Sim_nNue << std::endl;
                        std::cout << "Number of Simulated Muons: " << Sim_nMu << std::endl;
                        std::cout << "Energy Deposited by Muons: " << Sim_mu_Edep_b2 << std::endl;
                        std::cout << "Number of Simulated Muon neutrinos: "  << Sim_nNuMu << std::endl;
                        std::cout << "Number of Simulated Tau Particles: " << Sim_nTau << std::endl;
                        std::cout << "Number of Simulated Tau Neutrinos: " << Sim_nNutau << std::endl;
                        std::cout << "Number of Simulated Photons: " << Sim_nPhoton << std::endl;
                        std::cout << "Number of Simulated Pi0: " << Sim_nPionNeutral << std::endl;
                        std::cout << "Energy Deposited by Pi0: " << Sim_pi0_Edep_b2 << std::endl;
                        std::cout << "Number of Simulated Charged Pions: " << Sim_nPionCharged << std::endl;
                        std::cout << "Energy Deposited by Pi+: " << Sim_pip_Edep_b2 << std::endl;
                        std::cout << "Energy Deposited by Pi-: " << Sim_pim_Edep_b2 << std::endl;
                        std::cout << "Number of Simulated Protons: " << Sim_nProton << std::endl;
                        std::cout << "Energy Deposited by Protons: " << Sim_p_Edep_b2 << std::endl;
                        std::cout << "Number of Simulated Neutrons: " << Sim_nNeutron << std::endl;
                        std::cout << "Enery Deposited by Neutrons: " << Sim_n_Edep_b2 << std::endl;
                        //for(const auto& element : Vec.P_PDGVal){
                                //std::cout << element << ", " << std::endl;
                                //}
                        //std::cout << "PDG of Each Particle: ";
                        //printVector(*P_PDG);
                        //std::cout << "Mother of each particle (-1 is no mother): ";
                        //printVector(*P_mother);
                        //std::cout << "x Position of each Particle: ";
                        //PrintVector(*P_vtx_x);
                        //std::cout << "y Position of each Particle: ";
                        //PrintVector(*P_vtx_y);
                        //std::cout << "z Position of each Particle: ";
                        //PrintVector(*P_vtx_z);
                        std::cout << "----------------------------------------" << std::endl;

                }
        ntuple->Print();
        }                      
