/**
 * @file   MyEnergyAnalysis_module.cc
 * @brief  A file to read and analyze art::Event records from a DUNE FD MC file,
 * @author Wei Shi (wei.shi.1@stonybrook.edu)
 *
 * Adapted from https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/_AnalysisExample_
 */

// Include headers: starting from LArSoft and going up the software
// layers (nusimdata, art, etc.), ending with C++ is standard.

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <cmath>
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

namespace {

  // This is a local namespace. Stuff declared here will not be
  // visible beyond this file. We will define functions at the end,
  // but we declare them here so that the module can freely use them.

  // Utility function to get the diagonal of the detector
  double DetectorDiagonal(geo::GeometryCore const& geom);

  // Sort MC particles based on its start momentum P(0)
  bool MomentumOrderMCParticle(const simb::MCParticle*, const simb::MCParticle*);

  // Ancestor Mother is primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle&, int, std::map<int, const simb::MCParticle*>);

  //Ancestor Mother is Neutron
  bool IsAncestorMotherNeutron(const simb::MCParticle&, std::vector<int>, std::map<int, const simb::MCParticle*>);

  //Ancestor Mother is Proton
  bool IsAncestorMotherProton(const simb::MCParticle&, std::vector<int>, std::map<int, const simb::MCParticle*>);

  //Ancestor Mother is Pi+
  bool IsAncestorMotherPip(const simb::MCParticle&, std::vector<int>, std::map<int, const simb::MCParticle*>);

  //Ancestor Mother is Pi-
  bool IsAncestorMotherPim(const simb::MCParticle&, std::vector<int>, std::map<int, const simb::MCParticle*>);

  //Ancestor Mother is pi0
  bool IsAncestorMotherPi0(const simb::MCParticle&, std::vector<int>, std::map<int, const simb::MCParticle*>);


} // local namespace

// An outside package call this module like lar::example::MyEnergyAnalysis

namespace lar {
  namespace example {

    // BEGIN MyEnergyAnalysis group
    // -----------------------------------------------
    // class definition
    //
    // This class produces a ROOT tree that contains information
    // from the generated/simulated and reconstructed particles.
    //
    // Configuration parameters
    // =========================
    //
    // - GenieGenModuleLabel (string, default: "generator"): tag of the input data
    //   product with the event generator information
    //
    // - SimulationLabel (string, default: "largeant"): tag of the input data
    //   product with the detector simulation information (typically an instance
    //   of the LArG4 module)
    //
    class MyEnergyAnalysis : public art::EDAnalyzer {
    public:

      // This structure describes the configuration parameters of the module.
      // Any missing or unknown parameters will generate a configuration error.

      struct Config {

        // Save some typing:
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // One Atom for each parameter
        fhicl::Atom<art::InputTag> GenieGenModuleLabel{
          Name("GenieGenModuleLabel"),
          Comment("tag of the input data product with the event generator "
                  "information")};

        fhicl::Atom<art::InputTag> SimulationLabel{
          Name("SimulationLabel"),
          Comment("tag of the input data product with the detector simulation "
                  "information")};

      }; // Config

      using Parameters = art::EDAnalyzer::Table<Config>;

      /// Constructor: configures the module (see the Config structure above)
      explicit MyEnergyAnalysis(Parameters const& config);

      // This method is called once, at the start of the job. In this
      // example, it will define the histograms and n-tuples we'll
      // write.
      virtual void beginJob() override;

      // This method is called once, at the start of each run. It's a
      // good place to read databases or files that may have
      // run-dependent information.
      virtual void beginRun(const art::Run& run) override;

      // The analysis routine, called once per event.
      virtual void analyze(const art::Event& event) override;

    private:

      // The parameters we will read from the .fcl file.
      art::InputTag fGenieGenModuleLabel;     // The name of the producer that generated particles e.g. GENIE
      art::InputTag fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector

      // The n-tuple to create
      TTree* fNtuple;

      // Event info
      int fEvent;  // number of the event being processed
      int fRun;    // number of the run being processed
      int fSubRun; // number of the sub-run being processed

      // Add nu information
      double eP, eN, ePip, ePim, ePi0, eOther;    // Energy of particles
      int nLep, nP, nN, nPip, nPim, nPi0, nOther;                            // number of particles
      double E_vis_true;                 // True vis energy [GeV]

      //
      // Variables related to geneator/simulation
      //
      int fSimPDG;                       // MCParticle PDG ID
      std::vector<int> fSimP_TrackID_vec;
      std::vector<int> EDep_TrackID_vec;
      std::vector<int> fSimP_PDG_vec;
      std::vector<int> fSimP_Mom_vec;
      std::vector<int> fSimP_SC_vec;
      std::vector<float> fSimP_vtx_x_vec;
      std::vector<float> fSimP_vtx_y_vec;
      std::vector<float> fSimP_vtx_z_vec;
      std::vector<float> fSimP_ptot_vec;
      std::vector<float> fSimP_px_vec;
      std::vector<float> fSimP_py_vec;
      std::vector<float> fSimP_pz_vec;
      std::vector<float> fSimP_E_vec;
      std::vector<float> fSimP_M_vec;
      std::vector<float> fSimP_Ek_vec;

      int fSimTrackID;                   // GEANT ID of the particle being processed
      int EDepTrackID;

      int primarylep_trkID;
      std::vector<int> neutron_trkID;
      std::vector<int> proton_trkID;
      std::vector<int> pip_trkID;
      std::vector<int> pim_trkID;
      std::vector<int> pi0_trkID;
      // Next is to code it in vectors
      // std::vector<int> pi0_trkID;

      int fSim_nEle;                     // No. of Sim electrons (e+/e-)
      int fSim_nNue;                     // No. of Sim electron neutrinos (nue and nuebar)
      int fSim_nMu;                      // No. of Sim muons (mu+/mu-)
      int fSim_nNumu;                    // No. of Sim muon neutrinos (numu and numubar)
      int fSim_nTau;                     // No. of Sim tau leptons (+/-)
      int fSim_nNutau;                   // No. of Sim tau neutrinos (nutau and nutaubar)
      int fSim_nPhoton;                  // No. of Sim photons
      int fSim_nPionNeutral;             // No. of Sim pi+/pi-
      int fSim_nPip;
      int fSim_nPim;                     // No. of Sim pi0
      int fSim_nNeutron;                 // No. of Sim neutrons
      int fSim_nProton;                  // No. of Sim protons
      double fSim_LepE, fSim_HadE;       // Energy of Sim lep and had

      int fCCNC_truth;     		 //0=CC 1=NC
      int fMode_truth;     		 //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      int fInteractionType; 		 // Interaction type
      double fNuvtxx_truth; 		 //Genie true neutrino interaction vertex x
      double fNuvtxy_truth; 		 //Genie true neutrino interaction vertex y
      double fNuvtxz_truth;  		 //Genie true neutrino interaction vertex z
      int fNuPDG;                        // Generator level neutrino PDG code
      int fLepPDG;                       // Generator level outgoing lepton PDG code
      double fLepMomX, fLepMomY, fLepMomZ;      // Generator level outgoing lepton momentum
      double fLepvtx_x, fLepvtx_y, fLepvtx_z;               // Generator level outgoing lepton vtx
      double fVis_LepE;                      // Generator level neutrino lepton energy [GeV]
      double fLepMass;                      // Generator level neutrino lepton mass [GeV]
      int fStatusCode;                    // Generator level neutrino lepton statuscode
      double fLepNuAngle;                // Angle b/w nu and lepton

      double fGen_numu_E;                // Energy of generator level neutrino [GeV]
      double fSim_numu_E;                // Energy of leading muon (anti) neutrino

      double fSim_mu_start_vx;           // x position of the muon trajectory start
      double fSim_mu_start_vy;           // y .....................................
      double fSim_mu_start_vz;           // z .....................................
      std::vector<double> fSim_mu_start_4position; // (x,y,z,t) of the muon trajectory start
      double fSim_mu_end_vx;             // x position of the muon trajectory end
      double fSim_mu_end_vy;             // y ...................................
      double fSim_mu_end_vz;             // z ...................................
      std::vector<double> fSim_mu_end_4position;   // ................................ end
      double fSim_mu_start_px;           // x momentum of the muon trajectory start
      double fSim_mu_start_py;           // y .....................................
      double fSim_mu_start_pz;           // z .....................................
      double fSim_mu_start_E;            // Energy of leading mu
      std::vector<double> fSim_mu_start_4mommenta; // (Px,Py,Pz,E) of the muon trajectory start
      double fSim_mu_end_px;             // x momentum of the muon trajectory end
      double fSim_mu_end_py;             // y ...................................
      double fSim_mu_end_pz;             // z ...................................
      double fSim_mu_end_E;              // Energy of leading mu
      std::vector<double> fSim_mu_end_4mommenta;   // ................................... end
      double fSim_mu_track_length;       // leading mu track length

// For now these just store the first Particle created's Energy, since each particle would be far more complex

// Each Particle Processed's vertexes, momenta, and four vectors

      double fSim_P_start_vx;
      double fSim_P_start_vy;           
      double fSim_P_start_vz;           
      std::vector<double> fSim_P_start_4position;
      double fSim_P_end_vx;             
      double fSim_P_end_vy;             
      double fSim_P_end_vz;          
      std::vector<double> fSim_P_end_4position;   
      double fSim_P_start_px;           
      double fSim_P_start_py;         
      double fSim_P_start_pz;          
      double fSim_P_start_E;           
      std::vector<double> fSim_P_start_4mommenta; 
      double fSim_P_end_px;             
      double fSim_P_end_py;            
      double fSim_P_end_pz;          
      double fSim_P_end_E;             
      std::vector<double> fSim_P_end_4mommenta;  
      double fSim_P_track_length;    

      double fSim_N_start_vx;           
      double fSim_N_start_vy;           
      double fSim_N_start_vz;           
      std::vector<double> fSim_N_start_4position;
      double fSim_N_end_vx;             
      double fSim_N_end_vy;             
      double fSim_N_end_vz;          
      std::vector<double> fSim_N_end_4position;   
      double fSim_N_start_px;           
      double fSim_N_start_py;         
      double fSim_N_start_pz;          
      double fSim_N_start_E;           
      std::vector<double> fSim_N_start_4mommenta; 
      double fSim_N_end_px;             
      double fSim_N_end_py;            
      double fSim_N_end_pz;          
      double fSim_N_end_E;             
      std::vector<double> fSim_N_end_4mommenta;  
      double fSim_N_track_length;      

      double fSim_Pi0_start_vx;           
      double fSim_Pi0_start_vy;           
      double fSim_Pi0_start_vz;           
      std::vector<double> fSim_Pi0_start_4position;
      double fSim_Pi0_end_vx;             
      double fSim_Pi0_end_vy;             
      double fSim_Pi0_end_vz;          
      std::vector<double> fSim_Pi0_end_4position;   
      double fSim_Pi0_start_px;           
      double fSim_Pi0_start_py;         
      double fSim_Pi0_start_pz;          
      double fSim_Pi0_start_E;           
      std::vector<double> fSim_Pi0_start_4mommenta; 
      double fSim_Pi0_end_px;             
      double fSim_Pi0_end_py;            
      double fSim_Pi0_end_pz;          
      double fSim_Pi0_end_E;             
      std::vector<double> fSim_Pi0_end_4mommenta;  
      double fSim_Pi0_track_length;      

      double fSim_pip_start_vx;           
      double fSim_pip_start_vy;           
      double fSim_pip_start_vz;           
      std::vector<double> fSim_pip_start_4position;
      double fSim_pip_end_vx;             
      double fSim_pip_end_vy;             
      double fSim_pip_end_vz;          
      std::vector<double> fSim_pip_end_4position;   
      double fSim_pip_start_px;           
      double fSim_pip_start_py;         
      double fSim_pip_start_pz;          
      double fSim_pip_start_E;           
      std::vector<double> fSim_pip_start_4mommenta; 
      double fSim_pip_end_px;             
      double fSim_pip_end_py;            
      double fSim_pip_end_pz;          
      double fSim_pip_end_E;             
      std::vector<double> fSim_pip_end_4mommenta;  
      double fSim_pip_track_length;      

      double fSim_pim_start_vx;           
      double fSim_pim_start_vy;           
      double fSim_pim_start_vz;           
      std::vector<double> fSim_pim_start_4position;
      double fSim_pim_end_vx;             
      double fSim_pim_end_vy;             
      double fSim_pim_end_vz;          
      std::vector<double> fSim_pim_end_4position;   
      double fSim_pim_start_px;           
      double fSim_pim_start_py;         
      double fSim_pim_start_pz;          
      double fSim_pim_start_E;           
      std::vector<double> fSim_pim_start_4mommenta; 
      double fSim_pim_end_px;             
      double fSim_pim_end_py;            
      double fSim_pim_end_pz;          
      double fSim_pim_end_E;             
      std::vector<double> fSim_pim_end_4mommenta;  
      double fSim_pim_track_length;      


      std::vector<double> fSim_primary_end_energy;	//Primary particle in interaction's final energy
      std::vector<double> fSim_daughter_begin_energy;	//Sum of daughter particle's energy per interaction 

      double fSim_mu_Edep_b2;                // [MeV]
      double fSim_n_Edep_b2;                // [MeV] Energy Deposit of neutron
      double fSim_p_Edep_b2;                // [MeV] Energy Deposit of proton
      double fSim_pip_Edep_b2;                // [MeV] Energy Deposity of Pion+
      double fSim_pim_Edep_b2;                // [MeV] Energy Deposity of Pion-
      double fSim_pi0_Edep_b2;                // [MeV] Energy Deposity of Pion0
      double fSim_Other_Edep_b2;                // [MeV] Energy Deposity of eOther ; includes kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)

      // Two ways (a, b) to access collection plane +
      // Two ways (1, 2) of get E deposit for sim::IDE
      // Method b
      //double fSim_hadronic_Edep_b1;
      double fSim_hadronic_Edep_b2;
      //double fSim_hadronic_Edep_NonCollectionPlane_b2;  // [MeV]
      //double fSim_hadronic_Edep_b2_debug;          // [MeV]
      int fSim_n_hadronic_Edep_b;        // Number of hadronic energy deposits
      std::vector<float> fSim_hadronic_hit_x_b;
      std::vector<float> fSim_hadronic_hit_y_b;
      std::vector<float> fSim_hadronic_hit_z_b;
      //std::vector<float> fSim_hadronic_hit_Edep_b1;
      std::vector<float> fSim_hadronic_hit_Edep_b2;

      std::vector<std::string> fP_int_class_string;
      std::vector<unsigned long long> fP_int_class;

      //
      // Other variables that will be shared between different methods.
      //
      geo::GeometryCore const* fGeometryService; // pointer to Geometry provider
      double fElectronsToGeV;                    // conversion factor for no. of ionization electrons to energy deposited in GeV

      // True info for each particle generated
      std::vector<int> fP_PDG;                         // PDG code for each particle
      std::vector<int> fP_TrackID;                     // TrackID for each particle
      int fP_num;                         // Number of types of particle
      std::vector<int> fP_StatusCode;                  // Status code for each particle, https://internal.dunescience.org/doxygen/GENIEGen__module_8cc_source.html
      std::vector<float> fP_vtx_x;                    // Position: x component for each particle
      std::vector<float> fP_vtx_y;                    // Position: y component for each particle
      std::vector<float> fP_vtx_z;                    // Position: z component for each particle
      std::vector<float> fP_ptot;                     // Total momentum for each particle
      std::vector<float> fP_px;                       // Momentum: x component for each particle
      std::vector<float> fP_py;                       // Momentum: y component for each particle
      std::vector<float> fP_pz;                       // Momentum: z component for each particle
      std::vector<float> fP_E;                        // Energy for each particle [GeV]
      std::vector<float> fP_mass;                     // Mass for each particle [GeV/c^2]
      std::vector<float> fP_Ek;                       // Kinetic Energy for each particle [GeV]
      std::vector<int>   fP_mother;                    // Find the parent of the produced particle. -1 means this particle has no mother

      // True info for energy
      double fTrue_HadE;                              // True had E by adding all fP_E (!=lepton)
      double fTrue_LepE;                              // True Lep E by adding all fP_E (==lepton)
      double fVis_HadE;                               // Visible had E

    }; // class MyEnergyAnalysis

    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor

    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const& config)
      : EDAnalyzer(config)
      , fGenieGenModuleLabel(config().GenieGenModuleLabel())
      , fSimulationProducerLabel(config().SimulationLabel())
    {
      // Get a pointer to the geometry service provider.
      fGeometryService = lar::providerFrom<geo::Geometry>();

      // Tell beforehand all the data the module is going to read ("consumes") or
      // might read ("may_consume").
      consumes<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
      consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::beginJob()
    {
      // Get the detector length
      const double detectorLength = DetectorDiagonal(*fGeometryService);
      std::cout << "Detector length=" << detectorLength << " cm" << std::endl;

      // Access art's TFileService, which will handle creating and writing
      // histograms and n-tuples for us.
      art::ServiceHandle<art::TFileService const> tfs;

      // Define n-tuples
      fNtuple = tfs->make<TTree>("MyTree", "MyTree");

      fNtuple->Branch("Event",                    &fEvent,                  "Event/I");
      fNtuple->Branch("SubRun",                   &fSubRun,                 "SubRun/I");
      fNtuple->Branch("Run",                      &fRun,                    "Run/I");
      // Add true nu information
      fNtuple->Branch("Vis_LepE",          &fVis_LepE,         "Vis_LepE/D");
      fNtuple->Branch("LepMass",           &fLepMass,          "LepMass/D");

      fNtuple->Branch("eP",            &eP,            "eP/D");
      fNtuple->Branch("eN",            &eN,            "eN/D");
      fNtuple->Branch("ePip",          &ePip,          "ePip/D");
      fNtuple->Branch("ePim",          &ePim,          "ePim/D");
      fNtuple->Branch("ePi0",          &ePi0,          "ePi0/D");
      fNtuple->Branch("eOther",        &eOther,        "eOther/D");
      fNtuple->Branch("nLep",          &nLep,           "nLep/I");
      fNtuple->Branch("nP",            &nP,            "nP/I");
      fNtuple->Branch("nN",            &nN,            "nN/I");
      fNtuple->Branch("nPip",          &nPip,          "nPip/I");
      fNtuple->Branch("nPim",          &nPim,          "nPim/I");
      fNtuple->Branch("nPi0",          &nPi0,          "nPi0/I");
      fNtuple->Branch("nOther",        &nOther,        "nOther/D");
      fNtuple->Branch("E_vis_true",    &E_vis_true,    "E_vis_true/D");

      // GEN neutrino E
      fNtuple->Branch("Gen_numu_E",               &fGen_numu_E,             "Gen_numu_E/D");
      fNtuple->Branch("CCNC_truth",               &fCCNC_truth,             "CCNC_truth/I");
      fNtuple->Branch("Mode_truth",               &fMode_truth,             "Mode_truth/I");
      fNtuple->Branch("InteractionType",          &fInteractionType,        "InteractionType/I");
      fNtuple->Branch("Nuvtxx_truth",             &fNuvtxx_truth,           "Nuvtxx_truth/D");
      fNtuple->Branch("Nuvtxy_truth",             &fNuvtxy_truth,           "Nuvtxy_truth/D");
      fNtuple->Branch("Nuvtxz_truth",             &fNuvtxz_truth,           "Nuvtxz_truth/D");
      // Generator level PDG code
      fNtuple->Branch("LepPDG",        	          &fLepPDG,     	          "LepPDG/I");
      fNtuple->Branch("neuPDG",         	        &fNuPDG,             	    "neuPDG/I");
      fNtuple->Branch("LepNuAngle",               &fLepNuAngle,             "LepNuAngle/D");
      fNtuple->Branch("LepMomX",                  &fLepMomX,                "LepMomX/D");
      fNtuple->Branch("LepMomY",                  &fLepMomY,                "LepMomY/D");
      fNtuple->Branch("LepMomZ",                  &fLepMomZ,                "LepMomZ/D");
      fNtuple->Branch("Lepvtx_x",                 &fLepvtx_x,               "Lepvtx_x/D");
      fNtuple->Branch("Lepvtx_y",                 &fLepvtx_y,               "Lepvtx_y/D");
      fNtuple->Branch("Lepvtx_z",                 &fLepvtx_z,               "Lepvtx_z/D");
      fNtuple->Branch("StatusCode",         	    &fStatusCode,             "StatusCode/I");

      // Simulation branches Sim*
      fNtuple->Branch("SimP_TrackID_vec",              &fSimP_TrackID_vec);
      fNtuple->Branch("SimP_PDG_vec",                  &fSimP_PDG_vec);
      fNtuple->Branch("SimP_Mom_vec",                  &fSimP_Mom_vec);
      fNtuple->Branch("SimP_SC_vec",                   &fSimP_SC_vec);
      fNtuple->Branch("SimP_vtx_x_vec",                &fSimP_vtx_x_vec);
      fNtuple->Branch("SimP_vtx_y_vec",                &fSimP_vtx_y_vec);
      fNtuple->Branch("SimP_vtx_z_vec",                &fSimP_vtx_z_vec);
      fNtuple->Branch("SimP_ptot_vec",                 &fSimP_ptot_vec);
      fNtuple->Branch("SimP_px_vec",                   &fSimP_px_vec);
      fNtuple->Branch("SimP_py_vec",                   &fSimP_py_vec);
      fNtuple->Branch("SimP_pz_vec",                   &fSimP_pz_vec);
      fNtuple->Branch("SimP_E_vec",                    &fSimP_E_vec);
      fNtuple->Branch("SimP_M_vec",                    &fSimP_M_vec);
      fNtuple->Branch("SimP_Ek_vec",                   &fSimP_Ek_vec);

      fNtuple->Branch("Sim_nEle",                 &fSim_nEle,               "Sim_nEle/I");
      fNtuple->Branch("Sim_nNue",                 &fSim_nNue,               "Sim_nNue/I");
      fNtuple->Branch("Sim_nMu",                  &fSim_nMu,                "Sim_nMu/I");
      fNtuple->Branch("Sim_nNumu",                &fSim_nNumu,              "Sim_nNumu/I");
      fNtuple->Branch("Sim_nTau",                 &fSim_nTau,               "Sim_nTau/I");
      fNtuple->Branch("Sim_nNutau",               &fSim_nNutau,             "Sim_nNutau/I");
      fNtuple->Branch("Sim_nPhoton",              &fSim_nPhoton,            "Sim_nPhoton/I");
      fNtuple->Branch("Sim_nPionNeutral",         &fSim_nPionNeutral,       "Sim_nPionNeutral/I");
      fNtuple->Branch("Sim_nPip",                 &fSim_nPip,               "Sim_nPip/I");
      fNtuple->Branch("Sim_nPim",                 &fSim_nPim,                "Sim_nPim/I");
      fNtuple->Branch("Sim_nNeutron",             &fSim_nNeutron,           "Sim_nNeutron/I");
      fNtuple->Branch("Sim_nProton",              &fSim_nProton,            "Sim_nProton/I");
      fNtuple->Branch("Sim_LepE",                 &fSim_LepE,               "Sim_LepE/D");
      fNtuple->Branch("Sim_HadE",                 &fSim_HadE,               "Sim_HadE/D");

      // GEANT level neutrino E
      fNtuple->Branch("Sim_numu_E",               &fSim_numu_E,             "Sim_numu_E/D");
      // muon position
      fNtuple->Branch("Sim_mu_start_vx",          &fSim_mu_start_vx,        "Sim_mu_start_vx/D");
      fNtuple->Branch("Sim_mu_start_vy",          &fSim_mu_start_vy,        "Sim_mu_start_vy/D");
      fNtuple->Branch("Sim_mu_start_vz",          &fSim_mu_start_vz,        "Sim_mu_start_vz/D");
      fNtuple->Branch("Sim_mu_end_vx",            &fSim_mu_end_vx,          "Sim_mu_end_vx/D");
      fNtuple->Branch("Sim_mu_end_vy",            &fSim_mu_end_vy,          "Sim_mu_end_vy/D");
      fNtuple->Branch("Sim_mu_end_vz",            &fSim_mu_end_vz,          "Sim_mu_end_vz/D");
      fNtuple->Branch("Sim_mu_start_px",          &fSim_mu_start_px,        "Sim_mu_start_px/D");
      fNtuple->Branch("Sim_mu_start_py",          &fSim_mu_start_py,        "Sim_mu_start_py/D");
      fNtuple->Branch("Sim_mu_start_pz",          &fSim_mu_start_pz,        "Sim_mu_start_pz/D");
      fNtuple->Branch("Sim_mu_start_E",           &fSim_mu_start_E,         "Sim_mu_start_E/D");
      fNtuple->Branch("Sim_mu_end_px",            &fSim_mu_end_px,          "Sim_mu_end_px/D");
      fNtuple->Branch("Sim_mu_end_py",            &fSim_mu_end_py,          "Sim_mu_end_py/D");
      fNtuple->Branch("Sim_mu_end_pz",            &fSim_mu_end_pz,          "Sim_mu_end_pz/D");
      fNtuple->Branch("Sim_mu_end_E",             &fSim_mu_end_E,           "Sim_mu_end_E/D");
      fNtuple->Branch("Sim_mu_start_4position",   &fSim_mu_start_4position);
      fNtuple->Branch("Sim_mu_end_4position",     &fSim_mu_end_4position);
      fNtuple->Branch("Sim_mu_start_4mommenta",   &fSim_mu_start_4mommenta);
      fNtuple->Branch("Sim_mu_end_4mommenta",     &fSim_mu_end_4mommenta);
      fNtuple->Branch("Sim_mu_track_length",      &fSim_mu_track_length,    "Sim_mu_track_length/D");

      fNtuple->Branch("Sim_P_start_vx",          &fSim_P_start_vx,        "Sim_P_start_vx/D");
      fNtuple->Branch("Sim_P_start_vy",          &fSim_P_start_vy,        "Sim_P_start_vy/D");
      fNtuple->Branch("Sim_P_start_vz",          &fSim_P_start_vz,        "Sim_P_start_vz/D");
      fNtuple->Branch("Sim_P_end_vx",            &fSim_P_end_vx,          "Sim_P_end_vx/D");
      fNtuple->Branch("Sim_P_end_vy",            &fSim_P_end_vy,          "Sim_P_end_vy/D");
      fNtuple->Branch("Sim_P_end_vz",            &fSim_P_end_vz,          "Sim_P_end_vz/D");
      fNtuple->Branch("Sim_P_start_px",          &fSim_P_start_px,        "Sim_P_start_px/D");
      fNtuple->Branch("Sim_P_start_py",          &fSim_P_start_py,        "Sim_P_start_py/D");
      fNtuple->Branch("Sim_P_start_pz",          &fSim_P_start_pz,        "Sim_P_start_pz/D");
      fNtuple->Branch("Sim_P_start_E",           &fSim_P_start_E,         "Sim_P_start_E/D");
      fNtuple->Branch("Sim_P_end_px",            &fSim_P_end_px,          "Sim_P_end_px/D");
      fNtuple->Branch("Sim_P_end_py",            &fSim_P_end_py,          "Sim_P_end_py/D");
      fNtuple->Branch("Sim_P_end_pz",            &fSim_P_end_pz,          "Sim_P_end_pz/D");
      fNtuple->Branch("Sim_P_end_E",             &fSim_P_end_E,           "Sim_P_end_E/D");
      fNtuple->Branch("Sim_P_start_4position",   &fSim_P_start_4position);
      fNtuple->Branch("Sim_P_end_4position",     &fSim_P_end_4position);
      fNtuple->Branch("Sim_P_start_4mommenta",   &fSim_P_start_4mommenta);
      fNtuple->Branch("Sim_P_end_4mommenta",     &fSim_P_end_4mommenta);
      fNtuple->Branch("Sim_P_track_length",      &fSim_P_track_length,    "Sim_P_track_length/D");

      fNtuple->Branch("Sim_N_start_vx",          &fSim_N_start_vx,        "Sim_N_start_vx/D");
      fNtuple->Branch("Sim_N_start_vy",          &fSim_N_start_vy,        "Sim_N_start_vy/D");
      fNtuple->Branch("Sim_N_start_vz",          &fSim_N_start_vz,        "Sim_N_start_vz/D");
      fNtuple->Branch("Sim_N_end_vx",            &fSim_N_end_vx,          "Sim_N_end_vx/D");
      fNtuple->Branch("Sim_N_end_vy",            &fSim_N_end_vy,          "Sim_N_end_vy/D");
      fNtuple->Branch("Sim_N_end_vz",            &fSim_N_end_vz,          "Sim_N_end_vz/D");
      fNtuple->Branch("Sim_N_start_px",          &fSim_N_start_px,        "Sim_N_start_px/D");
      fNtuple->Branch("Sim_N_start_py",          &fSim_N_start_py,        "Sim_N_start_py/D");
      fNtuple->Branch("Sim_N_start_pz",          &fSim_N_start_pz,        "Sim_N_start_pz/D");
      fNtuple->Branch("Sim_N_start_E",           &fSim_N_start_E,         "Sim_N_start_E/D");
      fNtuple->Branch("Sim_N_end_px",            &fSim_N_end_px,          "Sim_N_end_px/D");
      fNtuple->Branch("Sim_N_end_py",            &fSim_N_end_py,          "Sim_N_end_py/D");
      fNtuple->Branch("Sim_N_end_pz",            &fSim_N_end_pz,          "Sim_N_end_pz/D");
      fNtuple->Branch("Sim_N_end_E",             &fSim_N_end_E,           "Sim_N_end_E/D");
      fNtuple->Branch("Sim_N_start_4position",   &fSim_N_start_4position);
      fNtuple->Branch("Sim_N_end_4position",     &fSim_N_end_4position);
      fNtuple->Branch("Sim_N_start_4mommenta",   &fSim_N_start_4mommenta);
      fNtuple->Branch("Sim_N_end_4mommenta",     &fSim_N_end_4mommenta);
      fNtuple->Branch("Sim_N_track_length",      &fSim_N_track_length,    "Sim_N_track_length/D");

      fNtuple->Branch("Sim_Pi0_start_vx",          &fSim_Pi0_start_vx,        "Sim_Pi0_start_vx/D");
      fNtuple->Branch("Sim_Pi0_start_vy",          &fSim_Pi0_start_vy,        "Sim_Pi0_start_vy/D");
      fNtuple->Branch("Sim_Pi0_start_vz",          &fSim_Pi0_start_vz,        "Sim_Pi0_start_vz/D");
      fNtuple->Branch("Sim_Pi0_end_vx",            &fSim_Pi0_end_vx,          "Sim_Pi0_end_vx/D");
      fNtuple->Branch("Sim_Pi0_end_vy",            &fSim_Pi0_end_vy,          "Sim_Pi0_end_vy/D");
      fNtuple->Branch("Sim_Pi0_end_vz",            &fSim_Pi0_end_vz,          "Sim_Pi0_end_vz/D");
      fNtuple->Branch("Sim_Pi0_start_px",          &fSim_Pi0_start_px,        "Sim_Pi0_start_px/D");
      fNtuple->Branch("Sim_Pi0_start_py",          &fSim_Pi0_start_py,        "Sim_Pi0_start_py/D");
      fNtuple->Branch("Sim_Pi0_start_pz",          &fSim_Pi0_start_pz,        "Sim_Pi0_start_pz/D");
      fNtuple->Branch("Sim_Pi0_start_E",           &fSim_Pi0_start_E,         "Sim_Pi0_start_E/D");
      fNtuple->Branch("Sim_Pi0_end_px",            &fSim_Pi0_end_px,          "Sim_Pi0_end_px/D");
      fNtuple->Branch("Sim_Pi0_end_py",            &fSim_Pi0_end_py,          "Sim_Pi0_end_py/D");
      fNtuple->Branch("Sim_Pi0_end_pz",            &fSim_Pi0_end_pz,          "Sim_Pi0_end_pz/D");
      fNtuple->Branch("Sim_Pi0_end_E",             &fSim_Pi0_end_E,           "Sim_Pi0_end_E/D");
      fNtuple->Branch("Sim_Pi0_start_4position",   &fSim_Pi0_start_4position);
      fNtuple->Branch("Sim_Pi0_end_4position",     &fSim_Pi0_end_4position);
      fNtuple->Branch("Sim_Pi0_start_4mommenta",   &fSim_Pi0_start_4mommenta);
      fNtuple->Branch("Sim_Pi0_end_4mommenta",     &fSim_Pi0_end_4mommenta);
      fNtuple->Branch("Sim_Pi0_track_length",      &fSim_Pi0_track_length,    "Sim_Pi0_track_length/D");

      fNtuple->Branch("Sim_pip_start_vx",          &fSim_pip_start_vx,        "Sim_pip_start_vx/D");
      fNtuple->Branch("Sim_pip_start_vy",          &fSim_pip_start_vy,        "Sim_pip_start_vy/D");
      fNtuple->Branch("Sim_pip_start_vz",          &fSim_pip_start_vz,        "Sim_pip_start_vz/D");
      fNtuple->Branch("Sim_pip_end_vx",            &fSim_pip_end_vx,          "Sim_pip_end_vx/D");
      fNtuple->Branch("Sim_pip_end_vy",            &fSim_pip_end_vy,          "Sim_pip_end_vy/D");
      fNtuple->Branch("Sim_pip_end_vz",            &fSim_pip_end_vz,          "Sim_pip_end_vz/D");
      fNtuple->Branch("Sim_pip_start_px",          &fSim_pip_start_px,        "Sim_pip_start_px/D");
      fNtuple->Branch("Sim_pip_start_py",          &fSim_pip_start_py,        "Sim_pip_start_py/D");
      fNtuple->Branch("Sim_pip_start_pz",          &fSim_pip_start_pz,        "Sim_pip_start_pz/D");
      fNtuple->Branch("Sim_pip_start_E",           &fSim_pip_start_E,         "Sim_pip_start_E/D");
      fNtuple->Branch("Sim_pip_end_px",            &fSim_pip_end_px,          "Sim_pip_end_px/D");
      fNtuple->Branch("Sim_pip_end_py",            &fSim_pip_end_py,          "Sim_pip_end_py/D");
      fNtuple->Branch("Sim_pip_end_pz",            &fSim_pip_end_pz,          "Sim_pip_end_pz/D");
      fNtuple->Branch("Sim_pip_end_E",             &fSim_pip_end_E,           "Sim_pip_end_E/D");
      fNtuple->Branch("Sim_pip_start_4position",   &fSim_pip_start_4position);
      fNtuple->Branch("Sim_pip_end_4position",     &fSim_pip_end_4position);
      fNtuple->Branch("Sim_pip_start_4mommenta",   &fSim_pip_start_4mommenta);
      fNtuple->Branch("Sim_pip_end_4mommenta",     &fSim_pip_end_4mommenta);
      fNtuple->Branch("Sim_pip_track_length",      &fSim_pip_track_length,    "Sim_pip_track_length/D");

      fNtuple->Branch("Sim_pim_start_vx",          &fSim_pim_start_vx,        "Sim_pim_start_vx/D");
      fNtuple->Branch("Sim_pim_start_vy",          &fSim_pim_start_vy,        "Sim_pim_start_vy/D");
      fNtuple->Branch("Sim_pim_start_vz",          &fSim_pim_start_vz,        "Sim_pim_start_vz/D");
      fNtuple->Branch("Sim_pim_end_vx",            &fSim_pim_end_vx,          "Sim_pim_end_vx/D");
      fNtuple->Branch("Sim_pim_end_vy",            &fSim_pim_end_vy,          "Sim_pim_end_vy/D");
      fNtuple->Branch("Sim_pim_end_vz",            &fSim_pim_end_vz,          "Sim_pim_end_vz/D");
      fNtuple->Branch("Sim_pim_start_px",          &fSim_pim_start_px,        "Sim_pim_start_px/D");
      fNtuple->Branch("Sim_pim_start_py",          &fSim_pim_start_py,        "Sim_pim_start_py/D");
      fNtuple->Branch("Sim_pim_start_pz",          &fSim_pim_start_pz,        "Sim_pim_start_pz/D");
      fNtuple->Branch("Sim_pim_start_E",           &fSim_pim_start_E,         "Sim_pim_start_E/D");
      fNtuple->Branch("Sim_pim_end_px",            &fSim_pim_end_px,          "Sim_pim_end_px/D");
      fNtuple->Branch("Sim_pim_end_py",            &fSim_pim_end_py,          "Sim_pim_end_py/D");
      fNtuple->Branch("Sim_pim_end_pz",            &fSim_pim_end_pz,          "Sim_pim_end_pz/D");
      fNtuple->Branch("Sim_pim_end_E",             &fSim_pim_end_E,           "Sim_pim_end_E/D");
      fNtuple->Branch("Sim_pim_start_4position",   &fSim_pim_start_4position);
      fNtuple->Branch("Sim_pim_end_4position",     &fSim_pim_end_4position);
      fNtuple->Branch("Sim_pim_start_4mommenta",   &fSim_pim_start_4mommenta);
      fNtuple->Branch("Sim_pim_end_4mommenta",     &fSim_pim_end_4mommenta);
      fNtuple->Branch("Sim_pim_track_length",      &fSim_pim_track_length,    "Sim_pim_track_length/D");

      fNtuple->Branch("Sim_primary_end_energy",    &fSim_primary_end_energy);
      fNtuple->Branch("Sim_daughter_begin_energy", &fSim_daughter_begin_energy);

      fNtuple->Branch("Sim_mu_Edep_b2",           &fSim_mu_Edep_b2,         "Sim_mu_Edep_b2/D");
      fNtuple->Branch("Sim_n_Edep_b2",            &fSim_n_Edep_b2,          "Sim_n_Edep_b2/D");
      fNtuple->Branch("Sim_p_Edep_b2",            &fSim_p_Edep_b2,          "Sim_p_Edep_b2/D");
      fNtuple->Branch("Sim_pip_Edep_b2",          &fSim_pip_Edep_b2,        "Sim_pip_Edep_b2/D");
      fNtuple->Branch("Sim_pim_Edep_b2",          &fSim_pim_Edep_b2,        "Sim_pim_Edep_b2/D");
      fNtuple->Branch("Sim_pi0_Edep_b2",          &fSim_pi0_Edep_b2,        "Sim_pi0_Edep_b2/D");
      fNtuple->Branch("Sim_Other_Edep_b2",        &fSim_Other_Edep_b2,      "Sim_Other_Edep_b2/D");

      fNtuple->Branch("Sim_hadronic_Edep_b2",     &fSim_hadronic_Edep_b2,   "Sim_hadronic_Edep_b2/D");
      fNtuple->Branch("Sim_n_hadronic_Edep_b",    &fSim_n_hadronic_Edep_b,  "Sim_n_hadronic_Edep_b/I");
      fNtuple->Branch("Sim_hadronic_hit_x_b",     &fSim_hadronic_hit_x_b);
      fNtuple->Branch("Sim_hadronic_hit_y_b",     &fSim_hadronic_hit_y_b);
      fNtuple->Branch("Sim_hadronic_hit_z_b",     &fSim_hadronic_hit_z_b);
      fNtuple->Branch("Sim_hadronic_hit_Edep_b2", &fSim_hadronic_hit_Edep_b2);

      // True info for each particle
      fNtuple->Branch("P_num",                    &fP_num,                "P_num/I");
      fNtuple->Branch("P_mother",                 &fP_mother);
      fNtuple->Branch("P_TrackID",                &fP_TrackID);
      fNtuple->Branch("P_PDG",        	          &fP_PDG);
      fNtuple->Branch("P_StatusCode",             &fP_StatusCode);
      fNtuple->Branch("P_vtx_x",                  &fP_vtx_x);
      fNtuple->Branch("P_vtx_y",                  &fP_vtx_y);
      fNtuple->Branch("P_vtx_z",                  &fP_vtx_z);
      fNtuple->Branch("P_ptot",                   &fP_ptot);
      fNtuple->Branch("P_px",                     &fP_px);
      fNtuple->Branch("P_py",                     &fP_py);
      fNtuple->Branch("P_pz",                     &fP_pz);
      fNtuple->Branch("P_E",                      &fP_E);
      fNtuple->Branch("P_mass",                   &fP_mass);
      fNtuple->Branch("P_Ek",                     &fP_Ek);

      // Reconstruction branches
      fNtuple->Branch("True_HadE",                &fTrue_HadE,            "True_HadE/D");
      fNtuple->Branch("True_LepE",                &fTrue_LepE,            "True_LepE/D");
      fNtuple->Branch("Vis_HadE",                 &fVis_HadE,             "Vis_HadE/D");
      
      fNtuple->Branch("P_int_class_string",       &fP_int_class_string);
      fNtuple->Branch("P_int_class",		  &fP_int_class);

    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::beginRun(const art::Run& /*run*/)
    {
      // Conversion factor for no. of ionization electrons to energy deposited in GeV
      // The ultimate source of this conversion factor is
      // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
      art::ServiceHandle<sim::LArG4Parameters const> larParameters;
      fElectronsToGeV = 1. / larParameters->GeVToElectrons();
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::analyze(const art::Event& event)
    {
      // Fetching basic event information.
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      // Initialize
      fGen_numu_E                = 0.;
      fCCNC_truth     		       = -9999.;
      fMode_truth     		       = -9999.;
      fInteractionType 		       = -9999.;
      fNuvtxx_truth 		         = -9999.;
      fNuvtxy_truth         		 = -9999.;
      fNuvtxz_truth  	        	 = -9999.;
      fSim_numu_E                = 0.;
	    
      fSim_mu_start_vx           = -9999.;
      fSim_mu_start_vy           = -9999.;
      fSim_mu_start_vz           = -9999.;
      fSim_mu_end_vx             = -9999.;
      fSim_mu_end_vy             = -9999.;
      fSim_mu_end_vz             = -9999.;
      fSim_mu_start_px           = -9999.;
      fSim_mu_start_py           = -9999.;
      fSim_mu_start_pz           = -9999.;
      fSim_mu_start_E            = -9999.;
      fSim_mu_end_px             = -9999.;
      fSim_mu_end_py             = -9999.;
      fSim_mu_end_pz             = -9999.;
      fSim_mu_end_E              = -9999.;
      fSim_mu_track_length       = -9999.;

      fSim_P_start_vx           = -9999.;
      fSim_P_start_vy           = -9999.;
      fSim_P_start_vz           = -9999.;
      fSim_P_end_vx             = -9999.;
      fSim_P_end_vy             = -9999.;
      fSim_P_end_vz             = -9999.;
      fSim_P_start_px           = -9999.;
      fSim_P_start_py           = -9999.;
      fSim_P_start_pz           = -9999.;
      fSim_P_start_E            = -9999.;
      fSim_P_end_px             = -9999.;
      fSim_P_end_py             = -9999.;
      fSim_P_end_pz             = -9999.;
      fSim_P_end_E              = -9999.;
      fSim_P_track_length       = -9999.;
	    
      fSim_N_start_vx           = -9999.;
      fSim_N_start_vy           = -9999.;
      fSim_N_start_vz           = -9999.;
      fSim_N_end_vx             = -9999.;
      fSim_N_end_vy             = -9999.;
      fSim_N_end_vz             = -9999.;
      fSim_N_start_px           = -9999.;
      fSim_N_start_py           = -9999.;
      fSim_N_start_pz           = -9999.;
      fSim_N_start_E            = -9999.;
      fSim_N_end_px             = -9999.;
      fSim_N_end_py             = -9999.;  
      fSim_N_end_pz             = -9999.;
      fSim_N_end_E              = -9999.;
      fSim_N_track_length       = -9999.;

      fSim_Pi0_start_vx           = -9999.;
      fSim_Pi0_start_vy           = -9999.;
      fSim_Pi0_start_vz           = -9999.;
      fSim_Pi0_end_vx             = -9999.;
      fSim_Pi0_end_vy             = -9999.;
      fSim_Pi0_end_vz             = -9999.;
      fSim_Pi0_start_px           = -9999.;
      fSim_Pi0_start_py           = -9999.;
      fSim_Pi0_start_pz           = -9999.;
      fSim_Pi0_start_E            = -9999.;
      fSim_Pi0_end_px             = -9999.;
      fSim_Pi0_end_py             = -9999.;
      fSim_Pi0_end_pz             = -9999.;
      fSim_Pi0_end_E              = -9999.;
      fSim_Pi0_track_length       = -9999.;

      fSim_pip_start_vx           = -9999.;
      fSim_pip_start_vy           = -9999.;
      fSim_pip_start_vz           = -9999.;
      fSim_pip_end_vx             = -9999.;
      fSim_pip_end_vy             = -9999.;
      fSim_pip_end_vz             = -9999.;
      fSim_pip_start_px           = -9999.;
      fSim_pip_start_py           = -9999.;
      fSim_pip_start_pz           = -9999.;
      fSim_pip_start_E            = -9999.;
      fSim_pip_end_px             = -9999.;
      fSim_pip_end_py             = -9999.;
      fSim_pip_end_pz             = -9999.;
      fSim_pip_end_E              = -9999.;
      fSim_pip_track_length       = -9999.;

      fSim_pim_start_vx           = -9999.;
      fSim_pim_start_vy           = -9999.;
      fSim_pim_start_vz           = -9999.;
      fSim_pim_end_vx             = -9999.;
      fSim_pim_end_vy             = -9999.;
      fSim_pim_end_vz             = -9999.;
      fSim_pim_start_px           = -9999.;
      fSim_pim_start_py           = -9999.;
      fSim_pim_start_pz           = -9999.;
      fSim_pim_start_E            = -9999.;
      fSim_pim_end_px             = -9999.;
      fSim_pim_end_py             = -9999.;
      fSim_pim_end_pz             = -9999.;
      fSim_pim_end_E              = -9999.;
      fSim_pim_track_length       = -9999.;
	    
      fSim_LepE                  = 0.;
      fSim_HadE                  = 0.;

    //Initialize track ID
      primarylep_trkID           = -1;
      neutron_trkID.clear();
      proton_trkID.clear();
      pip_trkID.clear();
      pim_trkID.clear();
      pi0_trkID.clear();

      // Initialize true info
      fLepNuAngle = -9999.;
      fLepMomX    = -9999.;
      fLepMomY    = -9999.;
      fLepMomZ    = -9999.;
      fLepvtx_x    = -9999.;
      fLepvtx_y    = -9999.;
      fLepvtx_z    = -9999.;
      fVis_LepE    = -9999.;
      fLepMass     = -9999.;

      fP_num        = 0;
      fP_PDG.clear();
      fP_mother.clear();
      fP_TrackID.clear();
      fP_StatusCode.clear();
      fP_vtx_x.clear();
      fP_vtx_y.clear();
      fP_vtx_z.clear();
      fP_ptot.clear();
      fP_px.clear();
      fP_py.clear();
      fP_pz.clear();
      fP_E.clear();
      fP_mass.clear();
      fP_Ek.clear();

      fSimP_TrackID_vec.clear();
      EDep_TrackID_vec.clear();
      fSimP_PDG_vec.clear();
      fSimP_Mom_vec.clear();
      fSimP_SC_vec.clear();
      fSimP_vtx_x_vec.clear();
      fSimP_vtx_y_vec.clear();
      fSimP_vtx_z_vec.clear();
      fSimP_ptot_vec.clear();
      fSimP_px_vec.clear();
      fSimP_py_vec.clear();
      fSimP_pz_vec.clear();
      fSimP_E_vec.clear();
      fSimP_M_vec.clear();
      fSimP_Ek_vec.clear();

      fSim_mu_Edep_b2        = 0.;
      fSim_n_Edep_b2         = 0.;
      fSim_p_Edep_b2         = 0.;
      fSim_pip_Edep_b2       = 0.;
      fSim_pim_Edep_b2       = 0.;
      fSim_pi0_Edep_b2       = 0.;
      fSim_Other_Edep_b2     = 0.;
      fSim_hadronic_Edep_b2  = 0.;

      fSim_mu_start_4position.clear();
      fSim_mu_end_4position.clear();
      fSim_mu_start_4mommenta.clear();
      fSim_mu_end_4mommenta.clear();

      fSim_P_start_4position.clear();
      fSim_P_end_4position.clear();
      fSim_P_start_4mommenta.clear();
      fSim_P_end_4mommenta.clear();

      fSim_N_start_4position.clear();
      fSim_N_end_4position.clear();
      fSim_N_start_4mommenta.clear();
      fSim_N_end_4mommenta.clear();

      fSim_Pi0_start_4position.clear();
      fSim_Pi0_end_4position.clear();
      fSim_Pi0_start_4mommenta.clear();
      fSim_Pi0_end_4mommenta.clear();

      fSim_pim_start_4position.clear();
      fSim_pim_end_4position.clear();
      fSim_pim_start_4mommenta.clear();
      fSim_pim_end_4position.clear();

      fSim_pip_start_4position.clear();
      fSim_pip_end_4position.clear();
      fSim_pip_start_4mommenta.clear();
      fSim_pip_end_4mommenta.clear();

      fSim_primary_end_energy.clear();
      fSim_daughter_begin_energy.clear();

      fP_int_class_string.clear();
      fP_int_class.clear();

      /*
        for (int i = 0; i < 4; i++) {
        fSim_mu_start_4position[i] = -9999.;
        fSim_mu_end_4position[i]   = -9999.;
        fSim_mu_start_4mommenta[i] = -9999.;
        fSim_mu_end_4mommenta[i]   = -9999.;

	fSim_P_start_4position[i] = -9999.;
	fSim_P_end_4position[i] = -9999.;
	fSim_P_start_4mommenta[i] = -9999.;
        fSim_P_end_4mommenta[i] = -9999.;

	fSim_N_start_4position[i] = -9999.;
	fSim_N_end_4position[i] = -9999.;
	fSim_N_start_4mommenta[i] = -9999.;
	fSim_N_end_4mommenta[i] = -9999.;

	fSim_Pi0_start_4position[i] = -9999.;
	fSim_Pi0_end_4position[i] = -9999.;
	fSim_Pi0_start_4mommenta[i] = -9999.;
	fSim_Pi0_end_4mommenta[i] = -9999.;

	fSim_pip_start_4position[i] = -9999.;
	fSim_pip_end_4position[i] = -9999.;
	fSim_pip_start_4mommenta[i] = -9999.;
	fSim_pip_end_4mommenta[i] = -9999.;

	fSim_pim_start_4position[i] = -9999.;
	fSim_pim_end_4position[i] = -9999.;
	fSim_pim_start_4mommenta[i] = -9999.;
	fSim_pim_end_4momenta[i] = -9999.;
	
      }
*/
      fSim_hadronic_hit_x_b.clear();
      fSim_hadronic_hit_y_b.clear();
      fSim_hadronic_hit_z_b.clear();

      fSim_hadronic_hit_Edep_b2.clear();

      // LArSoft data products: https://larsoft.org/important-concepts-in-larsoft/data-products/

      //
      // Process generator level info
      //

      // c.f. https://github.com/DUNE/dunetpc/blob/master/dune/FDSensOpt/CAFMaker_module.cc#L720
      //      https://github.com/DUNE/dunetpc/blob/master/dune/FDSensOpt/NueAna_module.cc#L639
      art::Handle<std::vector<simb::MCTruth>> mctruthListHandle; // Generator level truth
      std::vector<art::Ptr<simb::MCTruth>> mclist;
      if ( event.getByLabel(fGenieGenModuleLabel, mctruthListHandle) ) art::fill_ptr_vector(mclist, mctruthListHandle);

      // There could be more than one MCTruth, e.g., you might have multiple neutrino interactions per spill,
      // in which case you’d run GENIE multiple times and have one MCTruth per interaction.
      // Or you might want one MCTruth information for the GENIE event and another that overlays cosmic simulation or data onto the same event
      if ( mclist.size() )
      {
        fGen_numu_E = mclist[0]->GetNeutrino().Nu().E(); // true neutrino energy
        fCCNC_truth   = mclist[0]->GetNeutrino().CCNC(); // CC or NC interaction
      	fMode_truth   = mclist[0]->GetNeutrino().Mode(); // Interaction mode (QE/1-pi/DIS...)
        fInteractionType = mclist[0]->GetNeutrino().InteractionType(); // Interaction type
        fNuvtxx_truth = mclist[0]->GetNeutrino().Nu().Vx(); //Genie true neutrino interaction vertex x
	fNuvtxy_truth = mclist[0]->GetNeutrino().Nu().Vy(); //Genie true neutrino interaction vertex y
      	fNuvtxz_truth = mclist[0]->GetNeutrino().Nu().Vz(); //Genie true neutrino interaction vertex z
      	fNuPDG    = mclist[0]->GetNeutrino().Nu().PdgCode(); // Generator level neutrino PDG code
	fLepPDG     = mclist[0]->GetNeutrino().Lepton().PdgCode(); // Generator level lepton PDG code
        fLepMomX    = mclist[0]->GetNeutrino().Lepton().Momentum().X(); // Generator level lepton momentum x
        fLepMomY    = mclist[0]->GetNeutrino().Lepton().Momentum().Y();  // Generator level lepton momentum y
        fLepMomZ    = mclist[0]->GetNeutrino().Lepton().Momentum().Z();  // Generator level lepton momentum z
        fLepvtx_x       = mclist[0]->GetNeutrino().Lepton().Vx();   // Generator level lepton vtx x
        fLepvtx_y       = mclist[0]->GetNeutrino().Lepton().Vy();   // Generator level lepton vtx y
        fLepvtx_z       = mclist[0]->GetNeutrino().Lepton().Vz();   // Generator level lepton vtx z
        fLepMass        = mclist[0]->GetNeutrino().Lepton().Mass();
        fVis_LepE       = mclist[0]->GetNeutrino().Lepton().Momentum().T() - fLepMass; // Generator level neutrino lepton kinetic energy
        fStatusCode     = mclist[0]->GetNeutrino().Lepton().StatusCode();  // Generator level neutrino lepton statuscode
        fLepNuAngle = mclist[0]->GetNeutrino().Nu().Momentum().Vect().Angle(mclist[0]->GetNeutrino().Lepton().Momentum().Vect()); // Angle b/w nu and lepton
      }
      // Is evt vtx GetNeutrino().Nu().Vx()?


      // Add true particle counts

      eP = 0.;
      eN = 0.;
      ePip = 0.;
      ePim = 0.;
      ePi0 = 0.;
      eOther = 0.;

      nLep = 0;
      nP = 0;
      nN = 0;
      nPip = 0;
      nPim = 0;
      nPi0 = 0;
      nOther = 0;

      fP_num = mclist[0]->NParticles();
      // std::cout << "fP_num: " << fP_num << "\n\n";

      // Initialize
      fTrue_HadE = 0.;
      fTrue_LepE = 0.;
      fVis_HadE = 0.;
      double proton_mass = 0.93827; // GeV

      // Choose CC event only
      if ( fCCNC_truth ==0 )
      {
        for ( int p = 0; p < mclist[0]->NParticles(); p++ )
        {
          fP_TrackID.push_back(mclist[0]->GetParticle(p).TrackId());
          fP_PDG.push_back(mclist[0]->GetParticle(p).PdgCode());
          fP_mother.push_back(mclist[0]->GetParticle(p).Mother());
          fP_StatusCode.push_back(mclist[0]->GetParticle(p).StatusCode());
          fP_vtx_x.push_back(mclist[0]->GetParticle(p).Vx());
          fP_vtx_y.push_back(mclist[0]->GetParticle(p).Vy());
          fP_vtx_z.push_back(mclist[0]->GetParticle(p).Vz());
          fP_ptot.push_back(mclist[0]->GetParticle(p).P());
          fP_px.push_back(mclist[0]->GetParticle(p).Px());
          fP_py.push_back(mclist[0]->GetParticle(p).Py());
          fP_pz.push_back(mclist[0]->GetParticle(p).Pz());
          fP_E.push_back(mclist[0]->GetParticle(p).E());
          fP_mass.push_back(mclist[0]->GetParticle(p).Mass());
          fP_Ek.push_back(fP_E.at(p) - fP_mass.at(p));

          // Stable Final State
          // The sum of true energy of hadrons and leptons should be true nu energy minus binding energy
          // Paper related to the binding energy: https://link.springer.com/article/10.1140/epjc/s10052-019-6750-3
          // Calculate true vis had E
          if ( fP_StatusCode.at(p) == 1 ) // Stable Final State
          {

            // Calculate true Lep E
            if ( abs(fP_PDG.at(p)) == 13 )
            {
              fTrue_LepE += fP_E.at(p);
            }

            if ( abs(fP_PDG.at(p)) <= 999 && abs(fP_PDG.at(p)) >=100 ) // kPdgMeson
            {
              fTrue_HadE += fP_E.at(p);
            }
            else if (fP_PDG.at(p) == 2212 || fP_PDG.at(p) == 2112) // kPdgProton or kPdgNeutron
            {
              fTrue_HadE += fP_Ek.at(p);
            }
            else if ( fP_PDG.at(p) <= 9999 && fP_PDG.at(p) >= 1000  ) // kPdgBaryon except proton and neutron
            {
              fTrue_HadE += fP_Ek.at(p) + ( fP_mass.at(p) - proton_mass );
            }
            else if ( fP_PDG.at(p) >= -9999 && fP_PDG.at(p) <= -1000  ) // kPdgAntiBaryon except proton and neutron, antihyperon
            {
              fTrue_HadE += fP_Ek.at(p) + 2*fP_mass.at(p) + ( fP_mass.at(p) - proton_mass );
            }
            else if (fP_PDG.at(p) == 22) // kPdgGamma
            {
              fTrue_HadE += fP_E.at(p);
            }

            if ( abs(fP_PDG.at(p)) == 13 ) // kPdgMuon
            {
              nLep++;
            }
            if ( fP_PDG.at(p) == 2212 ) // kPdgProton
            {
              eP += fP_Ek.at(p);
              nP++;
            }
            else if ( fP_PDG.at(p) == 2112 ) // kPdgNeutron
            {
              eN += fP_Ek.at(p);
              nN++;
            }
            else if ( fP_PDG.at(p) == 211 ) // kPdgPiP
            {
              ePip += fP_Ek.at(p);
              nPip++;
            }
            else if ( fP_PDG.at(p) == -211 ) // kPdgPiM
            {
              ePim += fP_Ek.at(p);
              nPim++;
            }
            else if ( fP_PDG.at(p) == 111 ) // kPdgPi0
            {
              ePi0 += fP_Ek.at(p);
              nPi0++;
            }
            else if ( fP_PDG.at(p) == 321 || fP_PDG.at(p) == -321 || fP_PDG.at(p) == 311 || fP_PDG.at(p) == -311 || fP_PDG.at(p) == 130 || fP_PDG.at(p) == 310 || fP_PDG.at(p) == 22 || (fP_PDG.at(p)>=100 && fP_PDG.at(p)<=9999) || (fP_PDG.at(p)>=-9999 && fP_PDG.at(p)<=-100)) // kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)
            {
              eOther += fP_Ek.at(p);
              nOther++;
            }
          } //end kIStHadronInTheNucleus
        } // end mclist[0]->NParticles() loop

        // True visible energy:
        double pi0_mass = 0.134977; // GeV
        fVis_HadE = eP + ePip + ePim + ePi0 + eOther + nPi0 * pi0_mass;
        E_vis_true = fVis_LepE + fVis_HadE; // KE of leptons and hadrons
        // neutron will not deposit, so it cannot be counted in the E_vis_true
        // VisTrue_NDFD = LepE + HadE,
        // HadE = eP + ePip + ePim + ePi0 + (0.135 * nipi0) + eother

      } // end CC events selection

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // Get all the simulated channels for the event. These channels
      // include the energy deposited for each simulated track.
      auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

      // Create a map pf MCParticle to its track ID, to be used for hadronic part later
      std::map<int, const simb::MCParticle*> particleMap;
      // Create a map of energy deposits to its track ID
      std::map<int, double> EDepMap;

      //
      // Process Sim MCparticles info
      //

      art::Handle<std::vector<simb::MCParticle>> particleHandle; // GEANT 4 level truth

      // Then fill the vector with all the objects
      if ( !event.getByLabel(fSimulationProducerLabel, particleHandle) ) {
        // If no MCParticles in an event, throw an exception to force this module to stop.
        throw cet::exception("MyEnergyAnalysis") << " No simb::MCParticle objects in this event - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Store specific particles
      std::vector<const simb::MCParticle*> SimParticles;
      std::vector<const simb::MCParticle*> SimElectrons;
      std::vector<const simb::MCParticle*> SimNues;
      std::vector<const simb::MCParticle*> SimMuons;
      std::vector<const simb::MCParticle*> SimNumus;
      std::vector<const simb::MCParticle*> SimTaus;
      std::vector<const simb::MCParticle*> SimNutaus;
      std::vector<const simb::MCParticle*> SimPhotons;
      std::vector<const simb::MCParticle*> SimNeutralPions;
      std::vector<const simb::MCParticle*> SimPip;
      std::vector<const simb::MCParticle*> SimPim;
      std::vector<const simb::MCParticle*> SimNeutrons;
      std::vector<const simb::MCParticle*> SimProtons;

      // Loop over the list of particles in the event
      // GENIE: primary process; GEANT4: primary+secondary
      for ( auto const& particle : (*particleHandle) ) {

        // For the methods you can call for MCParticle, see ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCParticle.h.
        fSimTrackID = particle.TrackId();
        fSimP_TrackID_vec.push_back(fSimTrackID);
        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[fSimTrackID] = &particle;

        // Only for primary particles in the event
        fSimPDG = particle.PdgCode();
        fSimP_PDG_vec.push_back(fSimPDG);
        fSimP_Mom_vec.push_back(particle.Mother());
        fSimP_SC_vec.push_back(particle.StatusCode());
        fSimP_vtx_x_vec.push_back(particle.Vx());
        fSimP_vtx_y_vec.push_back(particle.Vy());
        fSimP_vtx_z_vec.push_back(particle.Vz());
        fSimP_ptot_vec.push_back(particle.P());
        fSimP_px_vec.push_back(particle.Px());
        fSimP_py_vec.push_back(particle.Py());
        fSimP_pz_vec.push_back(particle.Pz());
        fSimP_E_vec.push_back(particle.E());
        fSimP_M_vec.push_back(particle.Mass());
        fSimP_Ek_vec.push_back(particle.E()-particle.Mass());

        // Take note of primary lepton track id, to be used later
        if ( particle.Process() == "primary" && abs(fSimPDG) == 13 ) {
          primarylep_trkID = fSimTrackID;
          if (false) std::cout << "primarylep_trkID: " << primarylep_trkID << std::endl; // the primary lep should always have trk id = 1
        }

        //Take note of neutron trackID
        if ( fSimPDG == 2112){
          neutron_trkID.push_back(fSimTrackID);
          //std::cout << "neutron_trkID: " << neutron_trkID << std::endl;
        }

        //Take note of proton trackID
        if ( fSimPDG == 2212){
          proton_trkID.push_back(fSimTrackID);
          //std::cout << "proton_trkID: " << proton_trkID << std::endl;
        }

        // Take note of pip track ID
        if ( fSimPDG == 211 ) {
          pip_trkID.push_back(fSimTrackID);
          //std::cout << "pip_trkID: " << pip_trkID << std::endl;
        }

        // Take note of primary pim track ID
        if ( fSimPDG == -211 ) {
          pim_trkID.push_back(fSimTrackID);
          //std::cout << "pim_trkID: " << pim_trkID << std::endl;
        }

        // Take note of primary pi0 track ID
        if ( fSimPDG == 111 ) {
          pi0_trkID.push_back(fSimTrackID);
          // Could add an counter here to see how many pi0 in the event. If no pi0s, when calculate energy deposit later you don't need to check pi0 at all
          if (pi0_trkID.size() == 0) {
            fSim_pi0_Edep_b2 = 0;
          }

          //std::cout << "pi0_trkID: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
        }

        //if ( fSimPDG == 22 && particle.Mother() == pi0_trkID) {

          //std::cout << "gamma: Mother trkid: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
        //}

        // Calculate sim_lepE and sim_hadE
        if ( particle.StatusCode() == 1 )
        {
         // Sim_LepE
         if ( abs(fSimPDG) == 13 ) fSim_LepE += particle.E();
         // Sim_HadE
         if ( abs(fSimPDG) <= 999 && abs(fSimPDG) >=100 ) // kPdgMeson
            {
              fSim_HadE += particle.E();
            }
            else if (fSimPDG == 2212 || fSimPDG == 2112) // kPdgProton or kPdgNeutron
            {
              fSim_HadE += particle.E()-particle.Mass();
            }
            else if ( fSimPDG <= 9999 && fSimPDG >= 1000  ) // kPdgBaryon except proton and neutron
            {
              fSim_HadE += particle.E()-particle.Mass() + ( particle.Mass() - proton_mass );
            }
            else if ( fSimPDG >= -9999 && fSimPDG <= -1000  ) // kPdgAntiBaryon except proton and neutron, antihyperon
            {
              fSim_HadE += particle.E()-particle.Mass() + 2*particle.Mass() + ( particle.Mass() - proton_mass );
            }
            else if (fSimPDG == 22) // kPdgGamma
            {
              fSim_HadE += particle.E();
            }
        }
		SimParticles.push_back(&particle);
        //if ( particle.Process() == "primary" ) {
          if ( abs(fSimPDG) == 11 )   SimElectrons.push_back(&particle);
          if ( abs(fSimPDG) == 12 )   SimNues.push_back(&particle);
          if ( abs(fSimPDG) == 13 )   SimMuons.push_back(&particle);
          if ( abs(fSimPDG) == 14 )   SimNumus.push_back(&particle);
          if ( abs(fSimPDG) == 15 )   SimTaus.push_back(&particle);
          if ( abs(fSimPDG) == 16 )   SimNutaus.push_back(&particle);
          if ( abs(fSimPDG) == 22 )   SimPhotons.push_back(&particle);
          if ( abs(fSimPDG) == 111 )  SimNeutralPions.push_back(&particle);
          if ( fSimPDG == 211 )       SimPip.push_back(&particle);
          if ( fSimPDG == -211 )      SimPim.push_back(&particle);
          if ( abs(fSimPDG) == 2112 ) SimNeutrons.push_back(&particle);
          if ( abs(fSimPDG) == 2212 ) SimProtons.push_back(&particle);
        //}

      } // end loop over all particles in the event.

      fSim_nEle         = SimElectrons.size();
      fSim_nNue         = SimNues.size();
      fSim_nMu          = SimMuons.size();
      fSim_nNumu        = SimNumus.size();
      fSim_nTau         = SimTaus.size();
      fSim_nNutau       = SimNutaus.size();
      fSim_nPhoton      = SimPhotons.size();
      fSim_nPionNeutral = SimNeutralPions.size();
      fSim_nPip         = SimPip.size();
      fSim_nPim         = SimPim.size();
      fSim_nNeutron     = SimNeutrons.size();
      fSim_nProton      = SimProtons.size();

      // If multiple sim particles present in event, sort momentum from high to low
      if ( fSim_nNumu > 1 ) std::sort(SimNumus.begin(), SimNumus.end(), MomentumOrderMCParticle);
      if ( fSim_nMu > 1 )   std::sort(SimMuons.begin(), SimMuons.end(), MomentumOrderMCParticle);

      // Store info for leading E sim numu GEANT 4 level
      if ( fSim_nNumu > 0 ) {
        const simb::MCParticle& leadingnumu = *(SimNumus[0]);
        fSim_numu_E = leadingnumu.E();
      }

      // Store info for leading momentum sim muon
      if ( fSim_nMu > 0 ) {
        for (int i = 0; i < fSim_nMu; i++){

        const simb::MCParticle& leadingmu = *(SimMuons[i]);

        // A muon track consists of a set of 4-positions and 4-mommenta.
        const size_t munumberTrajectoryPoints = leadingmu.NumberTrajectoryPoints();

        // For trajectories, as for vectors and arrays, the first point is #0, not #1.
        const int mulast = munumberTrajectoryPoints - 1;
        const TLorentzVector& mupositionStart = leadingmu.Position(0);
        const TLorentzVector& mupositionEnd = leadingmu.Position(mulast);
        const TLorentzVector& mumomentumStart = leadingmu.Momentum(0);
        const TLorentzVector& mumomentumEnd = leadingmu.Momentum(mulast);

        fSim_mu_start_4position.push_back(mupositionStart.X());
		fSim_mu_start_4position.push_back(mupositionStart.Y());
		fSim_mu_start_4position.push_back(mupositionStart.Z());
		fSim_mu_start_4position.push_back(mupositionStart.T());
        fSim_mu_end_4position.push_back(mupositionEnd.X());
		fSim_mu_end_4position.push_back(mupositionEnd.Y());
		fSim_mu_end_4position.push_back(mupositionEnd.Z());
		fSim_mu_end_4position.push_back(mupositionEnd.T());
        fSim_mu_start_4mommenta.push_back(mumomentumStart.Px());
		fSim_mu_start_4mommenta.push_back(mumomentumStart.Py());
		fSim_mu_start_4mommenta.push_back(mumomentumStart.Pz());
		fSim_mu_start_4mommenta.push_back(mumomentumStart.E());
        fSim_mu_end_4mommenta.push_back(mumomentumEnd.Px());
		fSim_mu_end_4mommenta.push_back(mumomentumEnd.Py());
		fSim_mu_end_4mommenta.push_back(mumomentumEnd.Pz());
		fSim_mu_end_4mommenta.push_back(mumomentumEnd.E());

        }
        /*Fill position and momentum components
        fSim_mu_start_vx = leadingmu.Vx(0); // unit?
        fSim_mu_start_vy = leadingmu.Vy(0);
        fSim_mu_start_vz = leadingmu.Vz(0);
        fSim_mu_end_vx   = leadingmu.Vx(mulast);
        fSim_mu_end_vy   = leadingmu.Vy(mulast);
        fSim_mu_end_vz   = leadingmu.Vz(mulast);
        fSim_mu_start_px = leadingmu.Px(0);
        fSim_mu_start_py = leadingmu.Py(0);
        fSim_mu_start_pz = leadingmu.Pz(0);
        fSim_mu_start_E  = leadingmu.E(0);
        fSim_mu_end_px   = leadingmu.Px(mulast);
        fSim_mu_end_py   = leadingmu.Py(mulast);
        fSim_mu_end_pz   = leadingmu.Pz(mulast);
        fSim_mu_end_E    = leadingmu.E(mulast);
        */
        // Fill arrays with the 4-values.


        // Calculate length using spherical cooridnate system: assume straight track? time negligible?
        //const double mutrackLength = (mupositionEnd - mupositionStart).Rho();
        //fSim_mu_track_length = mutrackLength;
      
      }// End if muon exists

//Collecting four-vector information for each particle

	if ( fSim_nProton > 0 ) {
    for (int i = 0; i < fSim_nProton; i++){
		  const simb::MCParticle& leadingP = *(SimProtons[i]);

	    const size_t PnumberTrajectoryPoints = leadingP.NumberTrajectoryPoints();

	    const int Plast = PnumberTrajectoryPoints - 1;
	    const TLorentzVector& PpositionStart = leadingP.Position(0);
	    const TLorentzVector& PpositionEnd = leadingP.Position(Plast);
	    const TLorentzVector& PmomentumStart = leadingP.Momentum(0);
	    const TLorentzVector& PmomentumEnd = leadingP.Momentum(Plast);

        fSim_P_start_4position.push_back(PpositionStart.X());
		fSim_P_start_4position.push_back(PpositionStart.Y());
		fSim_P_start_4position.push_back(PpositionStart.Z());
		fSim_P_start_4position.push_back(PpositionStart.T());
        fSim_P_end_4position.push_back(PpositionEnd.X());
		fSim_P_end_4position.push_back(PpositionEnd.Y());
		fSim_P_end_4position.push_back(PpositionEnd.Z());
		fSim_P_end_4position.push_back(PpositionEnd.T());
        fSim_P_start_4mommenta.push_back(PmomentumStart.Px());
		fSim_P_start_4mommenta.push_back(PmomentumStart.Py());
		fSim_P_start_4mommenta.push_back(PmomentumStart.Pz());
		fSim_P_start_4mommenta.push_back(PmomentumStart.E());
        fSim_P_end_4mommenta.push_back(PmomentumEnd.Px());
		fSim_P_end_4mommenta.push_back(PmomentumEnd.Py());
		fSim_P_end_4mommenta.push_back(PmomentumEnd.Pz());
		fSim_P_end_4mommenta.push_back(PmomentumEnd.E());
    }
	      /*fSim_P_start_vx = leadingP.Vx(0); // unit?
        fSim_P_start_vy = leadingP.Vy(0);
        fSim_P_start_vz = leadingP.Vz(0);
        fSim_P_end_vx   = leadingP.Vx(Plast);
        fSim_P_end_vy   = leadingP.Vy(Plast);
        fSim_P_end_vz   = leadingP.Vz(Plast);
        fSim_P_start_px = leadingP.Px(0);
        fSim_P_start_py = leadingP.Py(0);
        fSim_P_start_pz = leadingP.Pz(0);
        fSim_P_start_E  = leadingP.E(0);
        fSim_P_end_px   = leadingP.Px(Plast);
        fSim_P_end_py   = leadingP.Py(Plast);
        fSim_P_end_pz   = leadingP.Pz(Plast);
        fSim_P_end_E    = leadingP.E(Plast);


	const double PtrackLength = (PpositionEnd - PpositionStart).Rho();
		fSim_P_track_length = PtrackLength;
    */
	}

	if ( fSim_nNeutron > 0 ) {
    for (int i = 0; i < fSim_nNeutron; i++){
		const simb::MCParticle& leadingN = *(SimNeutrons[i]);

		const size_t NnumberTrajectoryPoints = leadingN.NumberTrajectoryPoints();

	  const int Nlast = NnumberTrajectoryPoints - 1;
	  const TLorentzVector& NpositionStart = leadingN.Position(0);
	  const TLorentzVector& NpositionEnd = leadingN.Position(Nlast);
	  const TLorentzVector& NmomentumStart = leadingN.Momentum(0);
	  const TLorentzVector& NmomentumEnd = leadingN.Momentum(Nlast);

        fSim_N_start_4position.push_back(NpositionStart.X());
		fSim_N_start_4position.push_back(NpositionStart.Y());
		fSim_N_start_4position.push_back(NpositionStart.Z());
		fSim_N_start_4position.push_back(NpositionStart.T());
        fSim_mu_end_4position.push_back(NpositionEnd.X());
		fSim_N_end_4position.push_back(NpositionEnd.Y());
		fSim_N_end_4position.push_back(NpositionEnd.Z());
		fSim_N_end_4position.push_back(NpositionEnd.T());
        fSim_N_start_4mommenta.push_back(NmomentumStart.Px());
		fSim_N_start_4mommenta.push_back(NmomentumStart.Py());
		fSim_N_start_4mommenta.push_back(NmomentumStart.Pz());
		fSim_N_start_4mommenta.push_back(NmomentumStart.E());
        fSim_N_end_4mommenta.push_back(NmomentumEnd.Px());
		fSim_N_end_4mommenta.push_back(NmomentumEnd.Py());
		fSim_N_end_4mommenta.push_back(NmomentumEnd.Pz());
		fSim_N_end_4mommenta.push_back(NmomentumEnd.E());
    }
	      /*fSim_N_start_vx = leadingN.Vx(0); // unit?
        fSim_N_start_vy = leadingN.Vy(0);
        fSim_N_start_vz = leadingN.Vz(0);
        fSim_N_end_vx   = leadingN.Vx(Nlast);
        fSim_N_end_vy   = leadingN.Vy(Nlast);
        fSim_N_end_vz   = leadingN.Vz(Nlast);
        fSim_N_start_px = leadingN.Px(0);
        fSim_N_start_py = leadingN.Py(0);
        fSim_N_start_pz = leadingN.Pz(0);
        fSim_N_start_E  = leadingN.E(0);
        fSim_N_end_px   = leadingN.Px(Nlast);
        fSim_N_end_py   = leadingN.Py(Nlast);
        fSim_N_end_pz   = leadingN.Pz(Nlast);
        fSim_N_end_E    = leadingN.E(Nlast);

        NpositionStart.GetXYZT(fSim_N_start_4position);
        NpositionEnd.GetXYZT(fSim_N_end_4position);
        NmomentumStart.GetXYZT(fSim_N_start_4mommenta);
        NmomentumEnd.GetXYZT(fSim_N_end_4mommenta);

	const double NtrackLength = (NpositionEnd - NpositionStart).Rho();
		fSim_N_track_length = NtrackLength;
    */
	}

	if ( fSim_nPionNeutral > 0 ) {

    for (int i = 0; i < fSim_nPionNeutral; i++){
		  const simb::MCParticle& leadingPi0 = *(SimNeutralPions[i]);

	    const size_t Pi0numberTrajectoryPoints = leadingPi0.NumberTrajectoryPoints();

  	  const int Pi0last = Pi0numberTrajectoryPoints - 1;
	    const TLorentzVector& Pi0positionStart = leadingPi0.Position(0);
	    const TLorentzVector& Pi0positionEnd = leadingPi0.Position(Pi0last);
	    const TLorentzVector& Pi0momentumStart = leadingPi0.Momentum(0);
	    const TLorentzVector& Pi0momentumEnd = leadingPi0.Momentum(Pi0last);


        fSim_Pi0_start_4position.push_back(Pi0positionStart.X());
		fSim_Pi0_start_4position.push_back(Pi0positionStart.Y());
		fSim_Pi0_start_4position.push_back(Pi0positionStart.Z());
		fSim_Pi0_start_4position.push_back(Pi0positionStart.T());
        fSim_Pi0_end_4position.push_back(Pi0positionEnd.X());
		fSim_Pi0_end_4position.push_back(Pi0positionEnd.Y());
		fSim_Pi0_end_4position.push_back(Pi0positionEnd.Z());
		fSim_Pi0_end_4position.push_back(Pi0positionEnd.T());
        fSim_Pi0_start_4mommenta.push_back(Pi0momentumStart.Px());
		fSim_Pi0_start_4mommenta.push_back(Pi0momentumStart.Py());
		fSim_Pi0_start_4mommenta.push_back(Pi0momentumStart.Pz());
		fSim_Pi0_start_4mommenta.push_back(Pi0momentumStart.E());
        fSim_Pi0_end_4mommenta.push_back(Pi0momentumEnd.Px());
		fSim_Pi0_end_4mommenta.push_back(Pi0momentumEnd.Py());
		fSim_Pi0_end_4mommenta.push_back(Pi0momentumEnd.Pz());
		fSim_Pi0_end_4mommenta.push_back(Pi0momentumEnd.E());
    }
	/*fSim_Pi0_start_vx = leadingPi0.Vx(0); // unit?
        fSim_Pi0_start_vy = leadingPi0.Vy(0);
        fSim_Pi0_start_vz = leadingPi0.Vz(0);
        fSim_Pi0_end_vx   = leadingPi0.Vx(Pi0last);
        fSim_Pi0_end_vy   = leadingPi0.Vy(Pi0last);
        fSim_Pi0_end_vz   = leadingPi0.Vz(Pi0last);
        fSim_Pi0_start_px = leadingPi0.Px(0);
        fSim_Pi0_start_py = leadingPi0.Py(0);
        fSim_Pi0_start_pz = leadingPi0.Pz(0);
        fSim_Pi0_start_E  = leadingPi0.E(0);
        fSim_Pi0_end_px   = leadingPi0.Px(Pi0last);
        fSim_Pi0_end_py   = leadingPi0.Py(Pi0last);
        fSim_Pi0_end_pz   = leadingPi0.Pz(Pi0last);
        fSim_Pi0_end_E    = leadingPi0.E(Pi0last);

        Pi0positionStart.GetXYZT(fSim_Pi0_start_4position);
        Pi0positionEnd.GetXYZT(fSim_Pi0_end_4position);
        Pi0momentumStart.GetXYZT(fSim_Pi0_start_4mommenta);
        Pi0momentumEnd.GetXYZT(fSim_Pi0_end_4mommenta);

	const double Pi0trackLength = (Pi0positionEnd - Pi0positionStart).Rho();
		fSim_Pi0_track_length = Pi0trackLength;
    */
	}

  if (fSim_nPim > 0){
    for (int i = 0; i < fSim_nPim; i++){
     const simb::MCParticle& Pimvecs = *(SimPim[i]);
        
      const size_t pimnumberTrajectoryPoints = Pimvecs.NumberTrajectoryPoints();

	    const int pimlast = pimnumberTrajectoryPoints - 1;
	    const TLorentzVector& pimpositionStart = Pimvecs.Position(0);
	    const TLorentzVector& pimpositionEnd = Pimvecs.Position(pimlast);
	    const TLorentzVector& pimmomentumStart = Pimvecs.Momentum(0);
	    const TLorentzVector& pimmomentumEnd = Pimvecs.Momentum(pimlast);

        fSim_pim_start_4position.push_back(pimpositionStart.X());
		fSim_pim_start_4position.push_back(pimpositionStart.Y());
		fSim_pim_start_4position.push_back(pimpositionStart.Z());
		fSim_pim_start_4position.push_back(pimpositionStart.T());
        fSim_pim_end_4position.push_back(pimpositionEnd.X());
		fSim_pim_end_4position.push_back(pimpositionEnd.Y());
		fSim_pim_end_4position.push_back(pimpositionEnd.Z());
		fSim_pim_end_4position.push_back(pimpositionEnd.Y());
        fSim_pim_start_4mommenta.push_back(pimmomentumStart.Px());
		fSim_pim_start_4mommenta.push_back(pimmomentumStart.Py());
		fSim_pim_start_4mommenta.push_back(pimmomentumStart.Pz());
		fSim_pim_start_4mommenta.push_back(pimmomentumStart.E());
        fSim_pim_end_4mommenta.push_back(pimmomentumEnd.Px());
		fSim_pim_end_4mommenta.push_back(pimmomentumEnd.Py());
		fSim_pim_end_4mommenta.push_back(pimmomentumEnd.Pz());
		fSim_pim_end_4mommenta.push_back(pimmomentumEnd.E());
    }

       /* fSim_pim_start_vx = Pimvecs.Vx(0); // unit?
        fSim_pim_start_vy = Pimvecs.Vy(0);
        fSim_pim_start_vz = Pimvecs.Vz(0);
        fSim_pim_end_vx   = Pimvecs.Vx(pimlast);
        fSim_pim_end_vy   = Pimvecs.Vy(pimlast);
        fSim_pim_end_vz   = Pimvecs.Vz(pimlast);
        fSim_pim_start_px = Pimvecs.Px(0);
        fSim_pim_start_py = Pimvecs.Py(0);
        fSim_pim_start_pz = Pimvecs.Pz(0);
        fSim_pim_start_E  = Pimvecs.E(0);
        fSim_pim_end_px   = Pimvecs.Px(pimlast);
        fSim_pim_end_py   = Pimvecs.Py(pimlast);
        fSim_pim_end_pz   = Pimvecs.Pz(pimlast);
        fSim_pim_end_E    = Pimvecs.E(pimlast);

	const double pimtrackLength = (pimpositionEnd - pimpositionStart).Rho();
		fSim_pim_track_length = pimtrackLength;
    */
  }

  if (fSim_nPip > 0){
    for (int i = 0; i < fSim_nPip; i++){
      const simb::MCParticle& Pipvecs = *(SimPip[i]);
        
      const size_t pipnumberTrajectoryPoints = Pipvecs.NumberTrajectoryPoints();

	    const int piplast = pipnumberTrajectoryPoints - 1;
	    const TLorentzVector& pippositionStart = Pipvecs.Position(0);
	    const TLorentzVector& pippositionEnd = Pipvecs.Position(piplast);
	    const TLorentzVector& pipmomentumStart = Pipvecs.Momentum(0);
	    const TLorentzVector& pipmomentumEnd = Pipvecs.Momentum(piplast);

        fSim_pip_start_4position.push_back(pippositionStart.X());
		fSim_pip_start_4position.push_back(pippositionStart.Y());
		fSim_pip_start_4position.push_back(pippositionStart.Z());
		fSim_pip_start_4position.push_back(pippositionStart.T());
        fSim_pip_end_4position.push_back(pippositionEnd.X());
		fSim_pip_end_4position.push_back(pippositionEnd.Y());
		fSim_pip_end_4position.push_back(pippositionEnd.Z());
		fSim_pip_end_4position.push_back(pippositionEnd.T());
        fSim_pip_start_4mommenta.push_back(pipmomentumStart.Px());
		fSim_pip_start_4mommenta.push_back(pipmomentumStart.Py());
		fSim_pip_start_4mommenta.push_back(pipmomentumStart.Pz());
		fSim_pip_start_4mommenta.push_back(pipmomentumStart.E());
        fSim_pip_end_4mommenta.push_back(pipmomentumEnd.Px());
		fSim_pip_end_4mommenta.push_back(pipmomentumEnd.Py());
		fSim_pip_end_4mommenta.push_back(pipmomentumEnd.Pz());
		fSim_pip_end_4mommenta.push_back(pipmomentumEnd.E());
    }

      /*  fSim_pip_start_vx = Pipvecs.Vx(0); // unit?
        fSim_pip_start_vy = Pipvecs.Vy(0);
        fSim_pip_start_vz = Pipvecs.Vz(0);
        fSim_pip_end_vx   = Pipvecs.Vx(piplast);
        fSim_pip_end_vy   = Pipvecs.Vy(piplast);
        fSim_pip_end_vz   = Pipvecs.Vz(piplast);
        fSim_pip_start_px = Pipvecs.Px(0);
        fSim_pip_start_py = Pipvecs.Py(0);
        fSim_pip_start_pz = Pipvecs.Pz(0);
        fSim_pip_start_E  = Pipvecs.E(0);
        fSim_pip_end_px   = Pipvecs.Px(piplast);
        fSim_pip_end_py   = Pipvecs.Py(piplast);
        fSim_pip_end_pz   = Pipvecs.Pz(piplast);
        fSim_pip_end_E    = Pipvecs.E(piplast);


	const double piptrackLength = (pippositionEnd - pippositionStart).Rho();
		fSim_pip_track_length = piptrackLength;
    */
  }
  //End four-vector collection
  
  //Begin interaction classification and energy calculation
  /*
	    std::string combined_string = ""; 		//Stores the interaction classification code
	    std::string primary_particle = ""; 		//Primary particle of interaction
	    std::string daughter_particle = "";		//Current daughter particle being processed
	    std::string daughter_particles = "";	//String of daughter particles per primary particle
	    unsigned long long combined_int = 0;
	    double daughter_begin_sum = 0;
	    double primary_end_energy = 0;
//Loop through particle list and classify primary particle
	    for(size_t k = 0; k < fSimP_TrackID_vec.size(); k++){
		    	switch(fSimP_PDG_vec[k]){
				        case 22: primary_particle = "0";	//photon
                                                break;
                                        case 11: primary_particle = "1";	//electron
                                                break;
                                        case -11: primary_particle = "2";	//positron
                                                break;
                                        case 13: primary_particle = "3";	//muon 
                                                break;
                                        case -13: primary_particle = "4";	//mu+
                                                break;
                                        case 211: primary_particle = "5";	//pi+
                                                break;4
                                        case -211: primary_particle = "6";	//pi-
                                                break;
                                        case 2212: primary_particle = "7";	//proton
                                                break;
                                        case 2112: primary_particle = "8";	//neutron
                                                break;
                                        case 1000180400: primary_particle = "9"; //Argon nucleus
                                                break;
                                        default: break;
						}
			const simb::MCParticle& primaryVec = *(SimParticles[k]);	//Store primary particle's MCParticle information
			const size_t NPrimaryPoints = primaryVec.NumberTrajectoryPoints();	//Number of trajectory points for the primary particle
			/*const int primary_end = NPrimaryPoints - 1;
			const TLorentzVector& primary_end_4vector = primaryVec.Momentum(primary_end);
			double primary_end_energy = primary_end_4vector.E();
			if(primary_particle == "1" || "2" || "3" || "4" || "7" || "8" || "9"){
				primary_end_energy = primary_end_energy - primaryVec.Mass();
			}
			combined_string += primary_particle;
	if(primary_particle != "8"){						//Exclude neutron interactions for now
		for(size_t j = 0; j < fSimP_Mom_vec.size(); j++){
			if(fSimP_Mom_vec[j] == fSimP_TrackID_vec[k]){		//Loop through particles again to see which particles are tagged with primary as mom
				switch(fSimP_PDG_vec[j]){
					case 22: daughter_particle = "0";
						break;
					case 11: daughter_particle = "1";
						break;
					case -11: daughter_particle = "2";
						break;
					case 13: daughter_particle = "3";
						break;
					case -13: daughter_particle = "4";
						break;
					case 211: daughter_particle = "5";
						break;
					case -211: daughter_particle = "6";
						break;
                               		case 2212: daughter_particle = "7";
                                        	break;
                               		case 2112: daughter_particle = "8";
                          		        break;
					case 1000180400: daughter_particle = "9";
						break;	
					default: break;
				}
			const simb::MCParticle& daughterVec = *(SimParticles[j]);	//Daughter particle information
			for(size_t l = 0; l <= NPrimaryPoints; l++){
				const TLorentzVector& primary_position = primaryVec.Position(l);	//Store particles four-vectors
				const TLorentzVector& primary_momentum = primaryVec.Momentum(l);
				const TLorentzVector& daughter_position_start = daughterVec.Position(0)		//Match final primary position with initial daughter position
				if(primary_position.X() == daughter_position_start.X() && primary_position.Y() == daughter_position_start.Y() && primary_position.Z() == daughter_position_start.Z()){
				primary_end_energy = primary_momentum.E();	//Store Primary energy
				
				daughter_particles += daughter_particle;	//Store daughter
				daughter_particle= "";
			
			const TLorentzVector& daughter_begin_4vector = daughterVec.Momentum(0);
			double daughter_begin_energy = daughter_begin_4vector.E();
			if(daughter_particle == "1" || "2" || "3" || "4" || "7" || "8" || "9"){		//Subtract rest mass if not pion
				daughter_begin_energy = daughter_begin_energy - daughterVec.Mass();
			}
			daughter_begin_sum += daughter_begin_energy;		//sum daughter particle's energy
			daughter_begin_energy = 0;
		}
	    }
	    }
	    }
                                std::sort(daughter_particles.begin(), daughter_particles.end(), [](char a, char b){	//Sort daughter code from low to high mass
                                return std::stoull(std::string(1,a)) < std::stoull(std::string(1, b));
                        });
			combined_string += daughter_particles;
			
		  if(combined_string.length() <= 19 && combined_string.length() > 1){
		    fP_int_class_string.push_back(combined_string);
		    combined_int = std::stoull(combined_string);
		    fP_int_class.push_back(combined_int);
			std::cout << combined_string << std::endl;
		    fSim_primary_end_energy.push_back(primary_end_energy);
		    fSim_daughter_begin_energy.push_back(daughter_begin_sum);
			}
			}
		    combined_int = 0;
		    daughter_particles = "";
		    primary_particle = "";
		    combined_string = "";
		    daughter_begin_sum = 0;
		    primary_end_energy = 0;
	    
	    }
     */
      // Calculate sim hadronic deposit energy
      //

      // Loop over the SimChannel objects in the event to look at the energy deposited by particle's track.
      for ( auto const& channel : (*simChannelHandle) )
      {

        // Get the numeric ID associated with this channel.
        // See methods at https://internal.dunescience.org/doxygen/SimChannel_8h_source.html
        auto const channelNumber = channel.Channel();

        // Each channel has a map inside it that connects a time slice to energy deposits in the detector.
        // The full type of this map is std::map<unsigned short, std::vector<sim::IDE>>; we'll use "auto" here
        auto const& timeSlices = channel.TDCIDEMap();
        for ( auto const& timeSlice : timeSlices ) {

          // For the timeSlices map, the 'first' is a time slice number; The 'second' is a vector of IDE objects.
          auto const& energyDeposits = timeSlice.second;

          for ( auto const& energyDeposit : energyDeposits ) {

            // Method b: First check if it's on collection plane
            std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
            if ( Wires[0].planeID().Plane == 0 ) {

              // All EM shower are treated as secondary interactions, and their particles are not saved in the MC particle list
              // Still do the search, but now only for primary lepton (particleMap trkID is always positive)
              // Also search for EM shower particles from primary lepton, these deposits has trkID that's negative of the primary lepton trkID
              auto search = particleMap.find( abs(energyDeposit.trackID) );


              if ( search != particleMap.end() ) { // found match in map

                const simb::MCParticle& particle = *((*search).second);
                // if the energy deposit is from primary lepton,
                // or its ancestor mother particle is the primary lepton (e.g., from muon decays)
                if ( ( particle.Process() == "primary" && abs(particle.PdgCode()) == 13 ) || IsAncestorMotherPrimaryLep(particle, primarylep_trkID, particleMap) ) {
                  fSim_mu_Edep_b2 += energyDeposit.energy;
                  // now continue to the next energy deposit
                  // continue here to avoid counting into fSim_hadronic_Edep_b2
                  continue;
                } // end lepton deposited energy

                //if ( particle.PdgCode() == 22 && particle.Mother() == pi0_trkID ) {
                  //std::cout << "EDep MeV: "<< energyDeposit.energy << " from gamma: Mother trkid: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
                //}

                // if the energy deposit is from neutron
                // or its ancestor mother particle is the neutron
                if      ( particle.PdgCode() == 2112 || IsAncestorMotherNeutron(particle, neutron_trkID, particleMap) ) { fSim_n_Edep_b2 += energyDeposit.energy; } //std::cout << "fire n! " << std::endl;  // end neutron deposited energy
                // if the energy deposit is from proton
                // or its ancestor mother particle is the proton
                else if ( particle.PdgCode() == 2212 || IsAncestorMotherProton(particle, proton_trkID, particleMap) )   { fSim_p_Edep_b2 += energyDeposit.energy; } //std::cout << "fire p! " << std::endl;  // end proton deposited energy
                // if the energy deposit is from primary pip
                // or its ancestor mother particle is the pip
                else if ( particle.PdgCode() == 211  || IsAncestorMotherPip(particle, pip_trkID, particleMap) )         { fSim_pip_Edep_b2 += energyDeposit.energy; } // std::cout << "fire pip! pdg: "<< particle.PdgCode() << ", pip_trkID: " << pip_trkID << ", IsAncestorMotherPip: " << IsAncestorMotherPip(particle, pip_trkID, particleMap) << std::endl; } // end pip deposited energy
                // if the energy deposit is from primary pim
                // or its ancestor mother particle is the pim (pi0 decays into two photons)
                else if ( particle.PdgCode() == -211 || IsAncestorMotherPim(particle, pim_trkID, particleMap) )         { fSim_pim_Edep_b2 += energyDeposit.energy; } //std::cout << "fire pim! " << std::endl; } // end pim deposited energy
                // if the energy deposit is from primary pi0
                // or its ancestor mother particle is the pi0 (pi0 decays into two photons)
                else if ( particle.PdgCode() == 111  || IsAncestorMotherPi0(particle, pi0_trkID, particleMap) )         { fSim_pi0_Edep_b2 += energyDeposit.energy; } //std::cout << "fire pi0! " << std::endl; } // end pi0 deposited energy
                else if ( particle.PdgCode() == 321  || particle.PdgCode() == -321 || particle.PdgCode() == 311 || particle.PdgCode() == -311 || particle.PdgCode() == 130 || particle.PdgCode() == 310 || particle.PdgCode() == 22 || (particle.PdgCode() >= 100 && particle.PdgCode() <= -9999) || (particle.PdgCode() >= -9999 && particle.PdgCode() <= -100)) //eOther which includes: kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)
                {
                  fSim_Other_Edep_b2 += energyDeposit.energy;
                  //std::cout << "fire Other! " << std::endl;
                }
              } // end found match

              // if it's not, count as hadronic
              // If the energyDeposit made this far, it's counted as hadronic deposits (primary+secondary), do not involve particleMap
              fSim_hadronic_Edep_b2 += energyDeposit.energy;
              fSim_hadronic_hit_x_b.push_back(energyDeposit.x);
              fSim_hadronic_hit_y_b.push_back(energyDeposit.y);
              fSim_hadronic_hit_z_b.push_back(energyDeposit.z);
              fSim_hadronic_hit_Edep_b2.push_back(energyDeposit.energy);

              // Store trackID for debug
              EDepTrackID = energyDeposit.trackID;
              auto exist = EDepMap.find( EDepTrackID );
              // if can't find, store it and add the edep
              if ( exist == EDepMap.end() ) {
                EDep_TrackID_vec.push_back( EDepTrackID ); // negative trackID exists
                EDepMap[EDepTrackID] = energyDeposit.energy;
              } else { // find the same track id
                EDepMap[EDepTrackID] += energyDeposit.energy;
              }

            } // end plane == 0
          } // end energy deposit loop
        } // end For each time slice
      } // end For each SimChannel

      fSim_n_hadronic_Edep_b = fSim_hadronic_hit_x_b.size();

      // Print out EDepMap
      //if ( false ) std::cout << "fGen_numu_E: "<< fGen_numu_E << ", fSim_mu_Edep_b2_debug: " << fSim_mu_Edep_b2_debug<< ", fSim_hadronic_Edep_b2_debug: " << fSim_hadronic_Edep_b2_debug << ", Tot had E track id: " << EDep_TrackID_vec.size() << std::endl;
      if ( false ) {
        for (long unsigned int i=0; i< EDep_TrackID_vec.size(); i++) {
          std::cout << "Evt track id: " << EDep_TrackID_vec.at(i) << std::endl;
        }
        std::map<int, double>::iterator it;
        std::cout << "TrackID" << " | " << "Tot EDep" << std::endl;
        for (it = EDepMap.begin(); it != EDepMap.end(); it++) std::cout << "    " << it->first << " | " << it->second << std::endl;
      }

      // In general, objects in the LArSoft reconstruction chain are linked using the art::Assns class:
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artAssns>
      // The following statement will find the simb::MCTruth associated with the simb::MCParticle
      const art::FindManyP<simb::MCTruth> findManyTruth(particleHandle, event, fSimulationProducerLabel);

      if ( !findManyTruth.isValid() ) {
        std::cout << "findManyTruth simb::MCTruth for simb::MCParticle failed!" << std::endl;
      }

      size_t particle_index = 0; // only look at first particle in particleHandle's vector.
      auto const& truth = findManyTruth.at(particle_index);

      // Make sure there's no problem.
      if ( truth.empty() ) {
        std::cout << "Particle ID=" << particleHandle->at(particle_index).TrackId() << " has no primary!" << std::endl;
      }

      fNtuple->Fill();

    } // MyEnergyAnalysis::analyze()

    // This macro has to be defined for this module to be invoked from a
    // .fcl file; see MyEnergyAnalysis.fcl for more information.
    DEFINE_ART_MODULE(MyEnergyAnalysis)

  } // namespace example
} // namespace lar

// Back to our local namespace.
namespace {

  double DetectorDiagonal(geo::GeometryCore const& geom)
  {
    const double length = geom.DetLength();
    const double width = 2. * geom.DetHalfWidth();
    const double height = 2. * geom.DetHalfHeight();

    return std::sqrt(cet::sum_of_squares(length, width, height));
  }

  bool MomentumOrderMCParticle(const simb::MCParticle* p1, const simb::MCParticle* p2) {
    return ( p1->P(0) > p2->P(0) );
  }

  // If this returns true, then the energy deposit is associated with primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle& p1, int primarylep_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    // Immediate mother is the primary lep
    if ( MothertrkID == primarylep_trkID )  return true;
    // Immediate mother is not primary lep, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPrimaryLep(tmp_mother, primarylep_trkID, particleMap);
    }
  } // end GetAncestorMotherLeptonTrkID

    // If this returns true, then the energy deposit is associated with neutron
  bool IsAncestorMotherNeutron(const simb::MCParticle& p1, std::vector<int> neutron_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the neutrons in the event
    for (long unsigned int i=0; i< neutron_trkID.size(); i++){
      if ( MothertrkID == neutron_trkID.at(i) ) MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true) return true;

    // Immediate mother is not neutron, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherNeutron(tmp_mother, neutron_trkID, particleMap);
    }
  } // end GetAncestorMotherNeutronTrkID

      // If this returns true, then the energy deposit is associated with proton
  bool IsAncestorMotherProton(const simb::MCParticle& p1, std::vector<int> proton_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the protons in the event
    for (long unsigned int i=0; i< proton_trkID.size(); i++){
      if ( MothertrkID == proton_trkID.at(i) ) MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true) return true;

    // Immediate mother is not proton, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherProton(tmp_mother, proton_trkID, particleMap);
    }
  } // end GetAncestorMotherProtonTrkID

    // If this returns true, then the energy deposit is associated with pion+
  bool IsAncestorMotherPip(const simb::MCParticle& p1, std::vector<int> pip_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pips in the event
    for (long unsigned int i=0; i< pip_trkID.size(); i++){
      if ( MothertrkID == pip_trkID.at(i) ) MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true) return true;

    // Immediate mother is not pip, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPip(tmp_mother, pip_trkID, particleMap);
    }
  } // end GetAncestorMotherPipTrkID

  // If this returns true, then the energy deposit is associated with pion+
  bool IsAncestorMotherPim(const simb::MCParticle& p1, std::vector<int> pim_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pims in the event
    for (long unsigned int i=0; i< pim_trkID.size(); i++){
      if ( MothertrkID == pim_trkID.at(i) ) MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true) return true;


    // Immediate mother is not pim, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPim(tmp_mother, pim_trkID, particleMap);
    }
  } // end GetAncestorMotherPimTrkID


  // If this returns true, then the energy deposit is associated with pion-
  bool IsAncestorMotherPi0(const simb::MCParticle& p1, std::vector<int> pi0_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pi0s in the event
    for (long unsigned int i=0; i< pi0_trkID.size(); i++){
      if ( MothertrkID == pi0_trkID.at(i) ) MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true) return true;

    // Immediate mother is not pi0, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPi0(tmp_mother, pi0_trkID, particleMap);
    }
  } // end GetAncestorMotherPi0TrkID

} // local namespace
