#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


void loopTree(TTree* tree){
        TTreeReader reader(tree);

        TTreeReaderValue<int> eventReader(reader, "Event");
        TTreeReaderValue<std::vector<int>> FSmotherReader(reader, "P_mother");
        TTreeReaderValue<std::vector<int>> FSPDGReader(reader, "P_PDG");
        TTreeReaderValue<int> typeReader (reader, "Mode_truth");
        TTreeReaderValue<int> ccncReader (reader, "CCNC_truth");
        TTreeReaderValue<std::vector<int>> simPDGReader (reader, "SimP_PDG_vec");
        TTreeReaderValue<std::vector<int>> simMomReader (reader, "SimP_Mom_vec");

                while (reader.Next()) {
                        const int& EventVal = *eventReader;
                        const std::vector<int>& motherVal = *FSmotherReader;
                        const std::vector<int>& PDGVal = *FSPDGReader;
                        const int& TypeVal = *typeReader;
                        const int& ccncVal = *ccncReader;
                        const std::vector<int>& simPDGVal = *simPDGReader;
                        const std::vector<int>& simMomVal = *simMomReader;

                        std::cout << "Event #: " << EventVal << std::endl;

                        std::cout << "Interaction Type: ";
                        if(TypeVal == 0){
                                std::cout << "Quasi-Elastic" << std::endl;
                        }
                        else if (TypeVal == 1){
                                std::cout << "Resonant Production" << std::endl;
                        }
                        else if (TypeVal == 2){
                                std::cout << "Deep Inelastic Scattering" << std::endl;
                        }
                        else if (TypeVal == 3){
                                std::cout << "Coherent Production" << std::endl;
                        }

                        std::cout << "Final State Mother Particle: ";

                for (auto motherEntry : motherVal){
                        std::cout << motherEntry  << " ";
                }
                std::cout << std::endl;

                std::cout << "Final State PDG Code: ";
                for (auto PDGEntry : PDGVal){
                        std::cout << PDGEntry << " ";
                }
                std::cout << std::endl;

                std::cout << "Geant PDG: ";
                for (auto GPDGEntry : simPDGVal){
                        std::cout << GPDGEntry << " ";
                }
                std::cout << std::endl;

                std::cout << "Geant Mother: ";
                for(auto GMomEntry : simMomVal){
                        std::cout << GMomEntry << " ";
                }
                std::cout << std::endl;
                std::cout << "-------------------------------------------" << std::endl;
                }
        }

void VecLooker(){
        TFile* file = TFile::Open("test5files.root");

/*      TDirectory* myDir = file->GetDirectory("MyEnergyAnalysis");
if (myDir) {
  myDir->cd();
}
*/
        TTree* tree = (TTree*)file->Get("MyEnergyAnalysis/MyTree");

        if (tree) {
                loopTree(tree);
        }
        else{
                std::cerr << "Error: Could not find Tree" << std::endl;
        }
file->Close();
}
