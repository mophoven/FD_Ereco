
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
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

// C++ includes
#include <cmath>
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

namespace
{

  // This is a local namespace. Stuff declared here will not be
  // visible beyond this file. We will define functions at the end,
  // but we declare them here so that the module can freely use them.

  struct Vertex{
    float x, y, z, t;
    std::vector<const simb::MCParticle*> daughters;
  };

  // struct primaryVertex{
  //   float x, y, z, t;
  //   const simb::MCParticle* incoming;
  //   std::vector<const simb::MCParticle*> daughters;
  // };

  // Utility function to get the diagonal of the detector
  double DetectorDiagonal(geo::GeometryCore const &geom);

  // Sort MC particles based on its start momentum P(0)
  // bool MomentumOrderMCParticle(const simb::MCParticle*, const simb::MCParticle*);

  // Ancestor Mother is primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle &, int, std::map<int, const simb::MCParticle *>);

  // Ancestor Mother is Neutron
  bool IsAncestorMotherNeutron(const simb::MCParticle &, std::vector<int>, std::map<int, const simb::MCParticle *>);

  // Ancestor Mother is Proton
  bool IsAncestorMotherProton(const simb::MCParticle &, std::vector<int>, std::map<int, const simb::MCParticle *>);

  // Ancestor Mother is Pi+
  bool IsAncestorMotherPip(const simb::MCParticle &, std::vector<int>, std::map<int, const simb::MCParticle *>);

  // Ancestor Mother is Pi-
  bool IsAncestorMotherPim(const simb::MCParticle &, std::vector<int>, std::map<int, const simb::MCParticle *>);

  // Ancestor Mother is pi0
  bool IsAncestorMotherPi0(const simb::MCParticle &, std::vector<int>, std::map<int, const simb::MCParticle *>);

  // void getHadronicInformation(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&, int, double);

  void fillInteractionTree(const simb::MCParticle *, const Vertex &, const std::map<int, const simb::MCParticle *> &, TTree *,
                           float &, float &, float &, float &, float &, float &, float &, float &, int &, std::vector<float> &, std::vector<float> &,
                           std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<int> &);

  std::vector<Vertex> clusterVertices(const std::vector<const simb::MCParticle *> &);

  //void getHadronicInformation(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&, int, double);

  void fillInteractionTree(const simb::MCParticle*, const Vertex&, const std::map<int, const simb::MCParticle*>&, TTree*, 
                            float&, float&, float&, float&, float&, float&, float&, float&, int&, std::string&, std::vector<float>&, std::vector<float>&,
                            std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<int>&, std::vector<std::string>&);

  std::vector<Vertex> clusterVertices(const std::vector<const simb::MCParticle*>&);

  double getPrimaryKE(const simb::MCParticle*, double, double, double);

  void getHadronic02(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&, int&, int&, double&);

  void getHadronic02(const simb::MCParticle *, const std::vector<const simb::MCParticle *> &, int &, double &);

  //std::vector<primaryVertex> clusterPrimaryVertices(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&);
  void getDescendants(int, const std::vector<int> &, const std::vector<int> &, const std::map<int, const simb::MCParticle *> &, std::vector<const simb::MCParticle *> &);

  // std::vector<primaryVertex> clusterPrimaryVertices(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&);

} // local namespace

// An outside package call this module like lar::example::MyEnergyAnalysis

namespace lar
{
  namespace example
  {

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
    class MyEnergyAnalysis : public art::EDAnalyzer
    {
    public:
      // This structure describes the configuration parameters of the module.
      // Any missing or unknown parameters will generate a configuration error.

      struct Config
      {

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
      explicit MyEnergyAnalysis(Parameters const &config);

      // This method is called once, at the start of the job. In this
      // example, it will define the histograms and n-tuples we'll
      // write.
      virtual void beginJob() override;

      // This method is called once, at the start of each run. It's a
      // good place to read databases or files that may have
      // run-dependent information.
      virtual void beginRun(const art::Run &run) override;

      // The analysis routine, called once per event.
      virtual void analyze(const art::Event &event) override;

    private:
      // Step-KE histograms (one per particle type)
 
      // The parameters we will read from the .fcl file.
      art::InputTag fGenieGenModuleLabel;     // The name of the producer that generated particles e.g. GENIE
      art::InputTag fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector

      // The n-tuple to create
      TTree *fNtuple;

      TTree *fInteractionTree; // Tree for interaction information

      float fInX, fInY, fInZ, fInT;
      float fInPx, fInPy, fInPz, fInE;
      int fInPDG;

      std::vector<float> fOutX, fOutY, fOutZ, fOutT;
      std::vector<float> fOutPx, fOutPy, fOutPz, fOutE;
      std::vector<int> fOutPDG;

      TTree* fInteractionTree; // Tree for interaction information

      float fInX, fInY, fInZ, fInT;
      float fInPx, fInPy, fInPz, fInE;
      int fInPDG;
      std::string fInProcess;
      std::vector<std::string> fOutProcess;

      std::vector<float> fOutX, fOutY, fOutZ, fOutT;
      std::vector<float> fOutPx, fOutPy, fOutPz, fOutE;
      std::vector<int> fOutPDG;
      //std::vector<std::char> fOutProcess;

      // Event info
      int fEvent;  // number of the event being processed
      int fRun;    // number of the run being processed
      int fSubRun; // number of the sub-run being processed

      // Add nu information
      double eP, eN, ePip, ePim, ePi0, eOther;    // Energy of particles
      int nLep, nP, nN, nPip, nPim, nPi0, nOther; // number of particles
      double E_vis_true;                          // True vis energy [GeV]

      //
      // Variables related to geneator/simulation
      //
      int fSimPDG; // MCParticle PDG ID
      std::vector<int> fSimP_TrackID_vec;
      std::vector<int> EDep_TrackID_vec;
      std::vector<int> fSimP_PDG_vec;
      std::vector<int> fSimP_Mom_vec;
      std::vector<std::vector<int>> fSimP_Daughter_vec;
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

      int fSimTrackID; // GEANT ID of the particle being processed
      int EDepTrackID;

      int primarylep_trkID;
      std::vector<int> neutron_trkID;
      std::vector<int> proton_trkID;
      std::vector<int> pip_trkID;
      std::vector<int> pim_trkID;
      std::vector<int> pi0_trkID;
      // Next is to code it in vectors
      // std::vector<int> pi0_trkID;

      int fSim_nEle;         // No. of Sim electrons (e+/e-)
      int fSim_nNue;         // No. of Sim electron neutrinos (nue and nuebar)
      int fSim_nMu;          // No. of Sim muons (mu+/mu-)
      int fSim_nNumu;        // No. of Sim muon neutrinos (numu and numubar)
      int fSim_nTau;         // No. of Sim tau leptons (+/-)
      int fSim_nNutau;       // No. of Sim tau neutrinos (nutau and nutaubar)
      int fSim_nPhoton;      // No. of Sim photons
      int fSim_nPionNeutral; // No. of Sim pi+/pi-
      int fSim_nPip;
      int fSim_nPim;               // No. of Sim pi0
      int fSim_nNeutron;           // No. of Sim neutrons
      int fSim_nProton;            // No. of Sim protons
      double fSim_LepE, fSim_HadE; // Energy of Sim lep and had

      int fCCNC_truth;                        // 0=CC 1=NC
      int fMode_truth;                        // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      int fInteractionType;                   // Interaction type
      double fNuvtxx_truth;                   // Genie true neutrino interaction vertex x
      double fNuvtxy_truth;                   // Genie true neutrino interaction vertex y
      double fNuvtxz_truth;                   // Genie true neutrino interaction vertex z
      int fNuPDG;                             // Generator level neutrino PDG code
      int fLepPDG;                            // Generator level outgoing lepton PDG code
      double fLepMomX, fLepMomY, fLepMomZ;    // Generator level outgoing lepton momentum
      double fLepvtx_x, fLepvtx_y, fLepvtx_z; // Generator level outgoing lepton vtx
      double fVis_LepE;                       // Generator level neutrino lepton energy [GeV]
      double fLepMass;                        // Generator level neutrino lepton mass [GeV]
      int fStatusCode;                        // Generator level neutrino lepton statuscode
      double fLepNuAngle;                     // Angle b/w nu and lepton

      double fGen_numu_E; // Energy of generator level neutrino [GeV]
      double fSim_numu_E; // Energy of leading muon (anti) neutrino

      // For now these just store the first Particle created's Energy, since each particle would be far more complex

      // Each Particle Processed's vertexes, momenta, and four vectors

      int fSim_nParticles;

      std::vector<double> fSim_start_4position;
      std::vector<double> fSim_end_4position;
      std::vector<double> fSim_start_4mommenta;
      std::vector<double> fSim_end_4mommenta;

      std::vector<double> fSim_primary_end_energy;    // Primary particle in interaction's final energy
      std::vector<double> fSim_daughter_begin_energy; // Sum of daughter particle's energy per interaction

      double fSim_mu_Edep_b2;    // [MeV]
      double fSim_n_Edep_b2;     // [MeV] Energy Deposit of neutron
      double fSim_p_Edep_b2;     // [MeV] Energy Deposit of proton
      double fSim_pip_Edep_b2;   // [MeV] Energy Deposity of Pion+
      double fSim_pim_Edep_b2;   // [MeV] Energy Deposity of Pion-
      double fSim_pi0_Edep_b2;   // [MeV] Energy Deposity of Pion0
      double fSim_Other_Edep_b2; // [MeV] Energy Deposity of eOther ; includes kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)

      // Two ways (a, b) to access collection plane +
      // Two ways (1, 2) of get E deposit for sim::IDE
      // Method b
      // double fSim_hadronic_Edep_b1;
      double fSim_hadronic_Edep_b2;
      // double fSim_hadronic_Edep_NonCollectionPlane_b2;  // [MeV]
      // double fSim_hadronic_Edep_b2_debug;          // [MeV]
      int fSim_n_hadronic_Edep_b; // Number of hadronic energy deposits
      std::vector<float> fSim_hadronic_hit_x_b;
      std::vector<float> fSim_hadronic_hit_y_b;
      std::vector<float> fSim_hadronic_hit_z_b;
      // std::vector<float> fSim_hadronic_hit_Edep_b1;
      std::vector<float> fSim_hadronic_hit_Edep_b2;

      std::vector<std::string> fP_int_class_string;
      std::vector<unsigned long long> fP_int_class;

      //
      // Other variables that will be shared between different methods.
      //
      geo::GeometryCore const *fGeometryService; // pointer to Geometry provider
      double fElectronsToGeV;                    // conversion factor for no. of ionization electrons to energy deposited in GeV

      // True info for each particle generated
      std::vector<int> fP_PDG;        // PDG code for each particle
      std::vector<int> fP_TrackID;    // TrackID for each particle
      int fP_num;                     // Number of types of particle
      std::vector<int> fP_StatusCode; // Status code for each particle, https://internal.dunescience.org/doxygen/GENIEGen__module_8cc_source.html
      std::vector<float> fP_vtx_x;    // Position: x component for each particle
      std::vector<float> fP_vtx_y;    // Position: y component for each particle
      std::vector<float> fP_vtx_z;    // Position: z component for each particle
      std::vector<float> fP_ptot;     // Total momentum for each particle
      std::vector<float> fP_px;       // Momentum: x component for each particle
      std::vector<float> fP_py;       // Momentum: y component for each particle
      std::vector<float> fP_pz;       // Momentum: z component for each particle
      std::vector<float> fP_E;        // Energy for each particle [GeV]
      std::vector<float> fP_mass;     // Mass for each particle [GeV/c^2]
      std::vector<float> fP_Ek;       // Kinetic Energy for each particle [GeV]
      std::vector<int> fP_mother;     // Find the parent of the produced particle. -1 means this particle has no mother

      // True info for energy
      double fTrue_HadE; // True had E by adding all fP_E (!=lepton)
      double fTrue_LepE; // True Lep E by adding all fP_E (==lepton)
      double fVis_HadE;  // Visible had E

    }; // class MyEnergyAnalysis

    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor

    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const &config)
        : EDAnalyzer(config), fGenieGenModuleLabel(config().GenieGenModuleLabel()), fSimulationProducerLabel(config().SimulationLabel())
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
      fInteractionTree = tfs->make<TTree>("HadronicTree", "Handronic Interaction Information");

      fInteractionTree->Branch("InX", &fInX, "InX/F");
      fInteractionTree->Branch("InY", &fInY, "InY/F");
      fInteractionTree->Branch("InZ", &fInZ, "InZ/F");
      fInteractionTree->Branch("InT", &fInT, "InT/F");
      fInteractionTree->Branch("InPx", &fInPx, "InPx/F");
      fInteractionTree->Branch("InPy", &fInPy, "InPy/F");
      fInteractionTree->Branch("InPz", &fInPz, "InPz/F");
      fInteractionTree->Branch("InE", &fInE, "InE/F");
      fInteractionTree->Branch("InPDG", &fInPDG, "InPDG/I");
      fInteractionTree->Branch("InProcess", &fInProcess, "InProcess/C");
      fInteractionTree->Branch("OutX", &fOutX);
      fInteractionTree->Branch("OutY", &fOutY);
      fInteractionTree->Branch("OutZ", &fOutZ);
      fInteractionTree->Branch("OutT", &fOutT);
      fInteractionTree->Branch("OutPx", &fOutPx);
      fInteractionTree->Branch("OutPy", &fOutPy);
      fInteractionTree->Branch("OutPz", &fOutPz);
      fInteractionTree->Branch("OutE", &fOutE);
      fInteractionTree->Branch("OutPDG", &fOutPDG);
      fInteractionTree->Branch("OutProcess", &fOutProcess);

      fNtuple = tfs->make<TTree>("MyTree", "MyTree");

      fNtuple->Branch("Event", &fEvent, "Event/I");
      fNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
      fNtuple->Branch("Run", &fRun, "Run/I");
      // Add true nu information
      fNtuple->Branch("Vis_LepE", &fVis_LepE, "Vis_LepE/D");
      fNtuple->Branch("LepMass", &fLepMass, "LepMass/D");

      fNtuple->Branch("eP", &eP, "eP/D");
      fNtuple->Branch("eN", &eN, "eN/D");
      fNtuple->Branch("ePip", &ePip, "ePip/D");
      fNtuple->Branch("ePim", &ePim, "ePim/D");
      fNtuple->Branch("ePi0", &ePi0, "ePi0/D");
      fNtuple->Branch("eOther", &eOther, "eOther/D");
      fNtuple->Branch("nLep", &nLep, "nLep/I");
      fNtuple->Branch("nP", &nP, "nP/I");
      fNtuple->Branch("nN", &nN, "nN/I");
      fNtuple->Branch("nPip", &nPip, "nPip/I");
      fNtuple->Branch("nPim", &nPim, "nPim/I");
      fNtuple->Branch("nPi0", &nPi0, "nPi0/I");
      fNtuple->Branch("nOther", &nOther, "nOther/D");
      fNtuple->Branch("E_vis_true", &E_vis_true, "E_vis_true/D");

      // GEN neutrino E
      fNtuple->Branch("Gen_numu_E", &fGen_numu_E, "Gen_numu_E/D");
      fNtuple->Branch("CCNC_truth", &fCCNC_truth, "CCNC_truth/I");
      fNtuple->Branch("Mode_truth", &fMode_truth, "Mode_truth/I");
      fNtuple->Branch("InteractionType", &fInteractionType, "InteractionType/I");
      fNtuple->Branch("Nuvtxx_truth", &fNuvtxx_truth, "Nuvtxx_truth/D");
      fNtuple->Branch("Nuvtxy_truth", &fNuvtxy_truth, "Nuvtxy_truth/D");
      fNtuple->Branch("Nuvtxz_truth", &fNuvtxz_truth, "Nuvtxz_truth/D");
      // Generator level PDG code
      fNtuple->Branch("LepPDG", &fLepPDG, "LepPDG/I");
      fNtuple->Branch("neuPDG", &fNuPDG, "neuPDG/I");
      fNtuple->Branch("LepNuAngle", &fLepNuAngle, "LepNuAngle/D");
      fNtuple->Branch("LepMomX", &fLepMomX, "LepMomX/D");
      fNtuple->Branch("LepMomY", &fLepMomY, "LepMomY/D");
      fNtuple->Branch("LepMomZ", &fLepMomZ, "LepMomZ/D");
      fNtuple->Branch("Lepvtx_x", &fLepvtx_x, "Lepvtx_x/D");
      fNtuple->Branch("Lepvtx_y", &fLepvtx_y, "Lepvtx_y/D");
      fNtuple->Branch("Lepvtx_z", &fLepvtx_z, "Lepvtx_z/D");
      fNtuple->Branch("StatusCode", &fStatusCode, "StatusCode/I");

      // Simulation branches Sim*
      fNtuple->Branch("SimP_TrackID_vec", &fSimP_TrackID_vec);
      fNtuple->Branch("SimP_PDG_vec", &fSimP_PDG_vec);
      fNtuple->Branch("SimP_Mom_vec", &fSimP_Mom_vec);
      fNtuple->Branch("SimP_Daughter_vec", &fSimP_Daughter_vec);
      fNtuple->Branch("SimP_SC_vec", &fSimP_SC_vec);
      fNtuple->Branch("SimP_vtx_x_vec", &fSimP_vtx_x_vec);
      fNtuple->Branch("SimP_vtx_y_vec", &fSimP_vtx_y_vec);
      fNtuple->Branch("SimP_vtx_z_vec", &fSimP_vtx_z_vec);
      fNtuple->Branch("SimP_ptot_vec", &fSimP_ptot_vec);
      fNtuple->Branch("SimP_px_vec", &fSimP_px_vec);
      fNtuple->Branch("SimP_py_vec", &fSimP_py_vec);
      fNtuple->Branch("SimP_pz_vec", &fSimP_pz_vec);
      fNtuple->Branch("SimP_E_vec", &fSimP_E_vec);
      fNtuple->Branch("SimP_M_vec", &fSimP_M_vec);
      fNtuple->Branch("SimP_Ek_vec", &fSimP_Ek_vec);

      fNtuple->Branch("Sim_nEle", &fSim_nEle, "Sim_nEle/I");
      fNtuple->Branch("Sim_nNue", &fSim_nNue, "Sim_nNue/I");
      fNtuple->Branch("Sim_nMu", &fSim_nMu, "Sim_nMu/I");
      fNtuple->Branch("Sim_nNumu", &fSim_nNumu, "Sim_nNumu/I");
      fNtuple->Branch("Sim_nTau", &fSim_nTau, "Sim_nTau/I");
      fNtuple->Branch("Sim_nNutau", &fSim_nNutau, "Sim_nNutau/I");
      fNtuple->Branch("Sim_nPhoton", &fSim_nPhoton, "Sim_nPhoton/I");
      fNtuple->Branch("Sim_nPionNeutral", &fSim_nPionNeutral, "Sim_nPionNeutral/I");
      fNtuple->Branch("Sim_nPip", &fSim_nPip, "Sim_nPip/I");
      fNtuple->Branch("Sim_nPim", &fSim_nPim, "Sim_nPim/I");
      fNtuple->Branch("Sim_nNeutron", &fSim_nNeutron, "Sim_nNeutron/I");
      fNtuple->Branch("Sim_nProton", &fSim_nProton, "Sim_nProton/I");
      fNtuple->Branch("Sim_LepE", &fSim_LepE, "Sim_LepE/D");
      fNtuple->Branch("Sim_HadE", &fSim_HadE, "Sim_HadE/D");

      // GEANT level neutrino E
      fNtuple->Branch("Sim_numu_E", &fSim_numu_E, "Sim_numu_E/D");

      fNtuple->Branch("Sim_nParticles", &fSim_nParticles);

      fNtuple->Branch("Sim_start_4position", &fSim_start_4position);
      fNtuple->Branch("Sim_end_4position", &fSim_end_4position);
      fNtuple->Branch("Sim_start_4mommenta", &fSim_start_4mommenta);
      fNtuple->Branch("Sim_end_4mommenta", &fSim_end_4mommenta);

      fNtuple->Branch("Sim_primary_end_energy", &fSim_primary_end_energy);
      fNtuple->Branch("Sim_daughter_begin_energy", &fSim_daughter_begin_energy);

      fNtuple->Branch("Sim_mu_Edep_b2", &fSim_mu_Edep_b2, "Sim_mu_Edep_b2/D");
      fNtuple->Branch("Sim_n_Edep_b2", &fSim_n_Edep_b2, "Sim_n_Edep_b2/D");
      fNtuple->Branch("Sim_p_Edep_b2", &fSim_p_Edep_b2, "Sim_p_Edep_b2/D");
      fNtuple->Branch("Sim_pip_Edep_b2", &fSim_pip_Edep_b2, "Sim_pip_Edep_b2/D");
      fNtuple->Branch("Sim_pim_Edep_b2", &fSim_pim_Edep_b2, "Sim_pim_Edep_b2/D");
      fNtuple->Branch("Sim_pi0_Edep_b2", &fSim_pi0_Edep_b2, "Sim_pi0_Edep_b2/D");
      fNtuple->Branch("Sim_Other_Edep_b2", &fSim_Other_Edep_b2, "Sim_Other_Edep_b2/D");

      fNtuple->Branch("Sim_hadronic_Edep_b2", &fSim_hadronic_Edep_b2, "Sim_hadronic_Edep_b2/D");
      fNtuple->Branch("Sim_n_hadronic_Edep_b", &fSim_n_hadronic_Edep_b, "Sim_n_hadronic_Edep_b/I");
      fNtuple->Branch("Sim_hadronic_hit_x_b", &fSim_hadronic_hit_x_b);
      fNtuple->Branch("Sim_hadronic_hit_y_b", &fSim_hadronic_hit_y_b);
      fNtuple->Branch("Sim_hadronic_hit_z_b", &fSim_hadronic_hit_z_b);
      fNtuple->Branch("Sim_hadronic_hit_Edep_b2", &fSim_hadronic_hit_Edep_b2);

      // True info for each particle
      fNtuple->Branch("P_num", &fP_num, "P_num/I");
      fNtuple->Branch("P_mother", &fP_mother);
      fNtuple->Branch("P_TrackID", &fP_TrackID);
      fNtuple->Branch("P_PDG", &fP_PDG);
      fNtuple->Branch("P_StatusCode", &fP_StatusCode);
      fNtuple->Branch("P_vtx_x", &fP_vtx_x);
      fNtuple->Branch("P_vtx_y", &fP_vtx_y);
      fNtuple->Branch("P_vtx_z", &fP_vtx_z);
      fNtuple->Branch("P_ptot", &fP_ptot);
      fNtuple->Branch("P_px", &fP_px);
      fNtuple->Branch("P_py", &fP_py);
      fNtuple->Branch("P_pz", &fP_pz);
      fNtuple->Branch("P_E", &fP_E);
      fNtuple->Branch("P_mass", &fP_mass);
      fNtuple->Branch("P_Ek", &fP_Ek);

      // Reconstruction branches
      fNtuple->Branch("True_HadE", &fTrue_HadE, "True_HadE/D");
      fNtuple->Branch("True_LepE", &fTrue_LepE, "True_LepE/D");
      fNtuple->Branch("Vis_HadE", &fVis_HadE, "Vis_HadE/D");

      fNtuple->Branch("P_int_class_string", &fP_int_class_string);
      fNtuple->Branch("P_int_class", &fP_int_class);
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::beginRun(const art::Run & /*run*/)
    {
      // Conversion factor for no. of ionization electrons to energy deposited in GeV
      // The ultimate source of this conversion factor is
      // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
      art::ServiceHandle<sim::LArG4Parameters const> larParameters;
      fElectronsToGeV = 1. / larParameters->GeVToElectrons();
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::analyze(const art::Event &event)
    {
      // Fetching basic event information.
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      // Initialize
      fGen_numu_E = 0.;
      fCCNC_truth = -9999.;
      fMode_truth = -9999.;
      fInteractionType = -9999.;
      fNuvtxx_truth = -9999.;
      fNuvtxy_truth = -9999.;
      fNuvtxz_truth = -9999.;
      fSim_numu_E = 0.;

      fSim_LepE = 0.;
      fSim_HadE = 0.;

      // Initialize track ID
      primarylep_trkID = -1;
      neutron_trkID.clear();
      proton_trkID.clear();
      pip_trkID.clear();
      pim_trkID.clear();
      pi0_trkID.clear();

      // Initialize true info
      fLepNuAngle = -9999.;
      fLepMomX = -9999.;
      fLepMomY = -9999.;
      fLepMomZ = -9999.;
      fLepvtx_x = -9999.;
      fLepvtx_y = -9999.;
      fLepvtx_z = -9999.;
      fVis_LepE = -9999.;
      fLepMass = -9999.;

      fP_num = 0;
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
      fSimP_Daughter_vec.clear();
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

      fSim_mu_Edep_b2 = 0.;
      fSim_n_Edep_b2 = 0.;
      fSim_p_Edep_b2 = 0.;
      fSim_pip_Edep_b2 = 0.;
      fSim_pim_Edep_b2 = 0.;
      fSim_pi0_Edep_b2 = 0.;
      fSim_Other_Edep_b2 = 0.;
      fSim_hadronic_Edep_b2 = 0.;

      fSim_nParticles = 0;

      fSim_start_4position.clear();
      fSim_end_4position.clear();
      fSim_start_4mommenta.clear();
      fSim_end_4mommenta.clear();

      fSim_primary_end_energy.clear();
      fSim_daughter_begin_energy.clear();

      fP_int_class_string.clear();
      fP_int_class.clear();

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
      if (event.getByLabel(fGenieGenModuleLabel, mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);

      // There could be more than one MCTruth, e.g., you might have multiple neutrino interactions per spill,
      // in which case you’d run GENIE multiple times and have one MCTruth per interaction.
      // Or you might want one MCTruth information for the GENIE event and another that overlays cosmic simulation or data onto the same event
      if (mclist.size())
      {
        fGen_numu_E = mclist[0]->GetNeutrino().Nu().E();               // true neutrino energy
        fCCNC_truth = mclist[0]->GetNeutrino().CCNC();                 // CC or NC interaction
        fMode_truth = mclist[0]->GetNeutrino().Mode();                 // Interaction mode (QE/1-pi/DIS...)
        fInteractionType = mclist[0]->GetNeutrino().InteractionType(); // Interaction type
        fNuvtxx_truth = mclist[0]->GetNeutrino().Nu().Vx();            // Genie true neutrino interaction vertex x
        fNuvtxy_truth = mclist[0]->GetNeutrino().Nu().Vy();            // Genie true neutrino interaction vertex y
        fNuvtxz_truth = mclist[0]->GetNeutrino().Nu().Vz();            // Genie true neutrino interaction vertex z
        fNuPDG = mclist[0]->GetNeutrino().Nu().PdgCode();              // Generator level neutrino PDG code
        fLepPDG = mclist[0]->GetNeutrino().Lepton().PdgCode();         // Generator level lepton PDG code
        fLepMomX = mclist[0]->GetNeutrino().Lepton().Momentum().X();   // Generator level lepton momentum x
        fLepMomY = mclist[0]->GetNeutrino().Lepton().Momentum().Y();   // Generator level lepton momentum y
        fLepMomZ = mclist[0]->GetNeutrino().Lepton().Momentum().Z();   // Generator level lepton momentum z
        fLepvtx_x = mclist[0]->GetNeutrino().Lepton().Vx();            // Generator level lepton vtx x
        fLepvtx_y = mclist[0]->GetNeutrino().Lepton().Vy();            // Generator level lepton vtx y
        fLepvtx_z = mclist[0]->GetNeutrino().Lepton().Vz();            // Generator level lepton vtx z
        fLepMass = mclist[0]->GetNeutrino().Lepton().Mass();
        fVis_LepE = mclist[0]->GetNeutrino().Lepton().Momentum().T() - fLepMass;                                                  // Generator level neutrino lepton kinetic energy
        fStatusCode = mclist[0]->GetNeutrino().Lepton().StatusCode();                                                             // Generator level neutrino lepton statuscode
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
      if (fCCNC_truth == 0)
      {
        for (int p = 0; p < mclist[0]->NParticles(); p++)
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
          if (fP_StatusCode.at(p) == 1) // Stable Final State
          {

            // Calculate true Lep E
            if (abs(fP_PDG.at(p)) == 13)
            {
              fTrue_LepE += fP_E.at(p);
            }

            if (abs(fP_PDG.at(p)) <= 999 && abs(fP_PDG.at(p)) >= 100) // kPdgMeson
            {
              fTrue_HadE += fP_E.at(p);
            }
            else if (fP_PDG.at(p) == 2212 || fP_PDG.at(p) == 2112) // kPdgProton or kPdgNeutron
            {
              fTrue_HadE += fP_Ek.at(p);
            }
            else if (fP_PDG.at(p) <= 9999 && fP_PDG.at(p) >= 1000) // kPdgBaryon except proton and neutron
            {
              fTrue_HadE += fP_Ek.at(p) + (fP_mass.at(p) - proton_mass);
            }
            else if (fP_PDG.at(p) >= -9999 && fP_PDG.at(p) <= -1000) // kPdgAntiBaryon except proton and neutron, antihyperon
            {
              fTrue_HadE += fP_Ek.at(p) + 2 * fP_mass.at(p) + (fP_mass.at(p) - proton_mass);
            }
            else if (fP_PDG.at(p) == 22) // kPdgGamma
            {
              fTrue_HadE += fP_E.at(p);
            }

            if (abs(fP_PDG.at(p)) == 13) // kPdgMuon
            {
              nLep++;
            }
            if (fP_PDG.at(p) == 2212) // kPdgProton
            {
              eP += fP_Ek.at(p);
              nP++;
            }
            else if (fP_PDG.at(p) == 2112) // kPdgNeutron
            {
              eN += fP_Ek.at(p);
              nN++;
            }
            else if (fP_PDG.at(p) == 211) // kPdgPiP
            {
              ePip += fP_Ek.at(p);
              nPip++;
            }
            else if (fP_PDG.at(p) == -211) // kPdgPiM
            {
              ePim += fP_Ek.at(p);
              nPim++;
            }
            else if (fP_PDG.at(p) == 111) // kPdgPi0
            {
              ePi0 += fP_Ek.at(p);
              nPi0++;
            }
            else if (fP_PDG.at(p) == 321 || fP_PDG.at(p) == -321 || fP_PDG.at(p) == 311 || fP_PDG.at(p) == -311 || fP_PDG.at(p) == 130 || fP_PDG.at(p) == 310 || fP_PDG.at(p) == 22 || (fP_PDG.at(p) >= 100 && fP_PDG.at(p) <= 9999) || (fP_PDG.at(p) >= -9999 && fP_PDG.at(p) <= -100)) // kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)
            {
              eOther += fP_Ek.at(p);
              nOther++;
            }
          } // end kIStHadronInTheNucleus
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
      std::map<int, const simb::MCParticle *> particleMap;
      // Create a map of energy deposits to its track ID
      std::map<int, double> EDepMap;

      //
      // Process Sim MCparticles info
      //

      art::Handle<std::vector<simb::MCParticle>> particleHandle; // GEANT 4 level truth

      // Then fill the vector with all the objects
      if (!event.getByLabel(fSimulationProducerLabel, particleHandle))
      {
        // If no MCParticles in an event, throw an exception to force this module to stop.
        throw cet::exception("MyEnergyAnalysis") << " No simb::MCParticle objects in this event - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Store specific particles
      std::vector<const simb::MCParticle *> SimParticles;
      std::vector<const simb::MCParticle *> SimElectrons;
      std::vector<const simb::MCParticle *> SimNues;
      std::vector<const simb::MCParticle *> SimMuons;
      std::vector<const simb::MCParticle *> SimNumus;
      std::vector<const simb::MCParticle *> SimTaus;
      std::vector<const simb::MCParticle *> SimNutaus;
      std::vector<const simb::MCParticle *> SimPhotons;
      std::vector<const simb::MCParticle *> SimNeutralPions;
      std::vector<const simb::MCParticle *> SimPip;
      std::vector<const simb::MCParticle *> SimPim;
      std::vector<const simb::MCParticle *> SimNeutrons;
      std::vector<const simb::MCParticle *> SimProtons;

      // Loop over the list of particles in the event
      // GENIE: primary process; GEANT4: primary+secondary
      for (auto const &particle : (*particleHandle))
      {

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
        fSimP_Ek_vec.push_back(particle.E() - particle.Mass());

        // Take note of primary lepton track id, to be used later
        if (particle.Process() == "primary" && abs(fSimPDG) == 13)
        {
          primarylep_trkID = fSimTrackID;
          if (false)
            std::cout << "primarylep_trkID: " << primarylep_trkID << std::endl; // the primary lep should always have trk id = 1
        }

        // Take note of neutron trackID
        if (fSimPDG == 2112)
        {
          neutron_trkID.push_back(fSimTrackID);
          // std::cout << "neutron_trkID: " << neutron_trkID << std::endl;
        }

        // Take note of proton trackID
        if (fSimPDG == 2212)
        {
          proton_trkID.push_back(fSimTrackID);
          // std::cout << "proton_trkID: " << proton_trkID << std::endl;
        }

        // Take note of pip track ID
        if (fSimPDG == 211)
        {
          pip_trkID.push_back(fSimTrackID);
          // std::cout << "pip_trkID: " << pip_trkID << std::endl;
        }

        // Take note of primary pim track ID
        if (fSimPDG == -211)
        {
          pim_trkID.push_back(fSimTrackID);
          // std::cout << "pim_trkID: " << pim_trkID << std::endl;
        }

        // Take note of primary pi0 track ID
        if (fSimPDG == 111)
        {
          pi0_trkID.push_back(fSimTrackID);
          // Could add an counter here to see how many pi0 in the event. If no pi0s, when calculate energy deposit later you don't need to check pi0 at all
          if (pi0_trkID.size() == 0)
          {
            fSim_pi0_Edep_b2 = 0;
          }

          // std::cout << "pi0_trkID: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
        }

        // if ( fSimPDG == 22 && particle.Mother() == pi0_trkID) {

        // std::cout << "gamma: Mother trkid: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
        //}

        // Calculate sim_lepE and sim_hadE
        if (particle.StatusCode() == 1)
        {
          // Sim_LepE
          if (abs(fSimPDG) == 13)
            fSim_LepE += particle.E();
          // Sim_HadE
          if (abs(fSimPDG) <= 999 && abs(fSimPDG) >= 100) // kPdgMeson
          {
            fSim_HadE += particle.E();
          }
          else if (fSimPDG == 2212 || fSimPDG == 2112) // kPdgProton or kPdgNeutron
          {
            fSim_HadE += particle.E() - particle.Mass();
          }
          else if (fSimPDG <= 9999 && fSimPDG >= 1000) // kPdgBaryon except proton and neutron
          {
            fSim_HadE += particle.E() - particle.Mass() + (particle.Mass() - proton_mass);
          }
          else if (fSimPDG >= -9999 && fSimPDG <= -1000) // kPdgAntiBaryon except proton and neutron, antihyperon
          {
            fSim_HadE += particle.E() - particle.Mass() + 2 * particle.Mass() + (particle.Mass() - proton_mass);
          }
          else if (fSimPDG == 22) // kPdgGamma
          {
            fSim_HadE += particle.E();
          }
        }
        SimParticles.push_back(&particle);
        // if ( particle.Process() == "primary" ) {
        if (abs(fSimPDG) == 11)
          SimElectrons.push_back(&particle);
        if (abs(fSimPDG) == 12)
          SimNues.push_back(&particle);
        if (abs(fSimPDG) == 13)
          SimMuons.push_back(&particle);
        if (abs(fSimPDG) == 14)
          SimNumus.push_back(&particle);
        if (abs(fSimPDG) == 15)
          SimTaus.push_back(&particle);
        if (abs(fSimPDG) == 16)
          SimNutaus.push_back(&particle);
        if (abs(fSimPDG) == 22)
          SimPhotons.push_back(&particle);
        if (abs(fSimPDG) == 111)
          SimNeutralPions.push_back(&particle);
        if (fSimPDG == 211)
          SimPip.push_back(&particle);
        if (fSimPDG == -211)
          SimPim.push_back(&particle);
        if (abs(fSimPDG) == 2112)
          SimNeutrons.push_back(&particle);
        if (abs(fSimPDG) == 2212)
          SimProtons.push_back(&particle);
        //}

      } // end loop over all particles in the event.

      fSim_nEle = SimElectrons.size();
      fSim_nNue = SimNues.size();
      fSim_nMu = SimMuons.size();
      fSim_nNumu = SimNumus.size();
      fSim_nTau = SimTaus.size();
      fSim_nNutau = SimNutaus.size();
      fSim_nPhoton = SimPhotons.size();
      fSim_nPionNeutral = SimNeutralPions.size();
      fSim_nPip = SimPip.size();
      fSim_nPim = SimPim.size();
      fSim_nNeutron = SimNeutrons.size();
      fSim_nProton = SimProtons.size();
      fSim_nParticles = SimParticles.size();

      // Store info for leading E sim numu GEANT 4 level
      for (int i = 0; i < fSim_nParticles; i++)
      {
        

        // const int last = Ntrajpoints - 1;
        // const TLorentzVector& positionStart = particleVec.Position(0);
        // const TLorentzVector& positionEnd = particleVec.Position(last);
        // const TLorentzVector& momentumStart = particleVec.Momentum(0);
        // const TLorentzVector& momentumEnd = particleVec.Momentum(last);
        //  New stuff
        // double fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;
        // auto const &geom = *fGeometryService;
        // fXmin = -geom.DetLength();
        // fXmax = geom.DetLength();
        // fYmin = -geom.DetHalfWidth()*2;
        // fYmax = geom.DetHalfWidth()*2;
        // fZmin = -geom.DetHalfHeight()*2;
        // fZmax = geom.DetHalfHeight()*2;
        // std::cout << fXmax << " ," << fYmin << "," << fYmax << "," << fZmin << "," << fZmax << std::endl;

        // 2) Loop over each particle
        // for (int l=0; l<fSim_nParticles; l++) {

        /*fSim_start_4position.push_back(positionStart.X());
    fSim_start_4position.push_back(positionStart.Y());
    fSim_start_4position.push_back(positionStart.Z());
    fSim_start_4position.push_back(positionStart.T());
        fSim_end_4position.push_back(positionEnd.X());
    fSim_end_4position.push_back(positionEnd.Y());
    fSim_end_4position.push_back(positionEnd.Z());
    fSim_end_4position.push_back(positionEnd.T());
        fSim_start_4mommenta.push_back(momentumStart.Px());
    fSim_start_4mommenta.push_back(momentumStart.Py());
    fSim_start_4mommenta.push_back(momentumStart.Pz());
    fSim_start_4mommenta.push_back(momentumStart.E());
        fSim_end_4mommenta.push_back(momentumEnd.Px());
    fSim_end_4mommenta.push_back(momentumEnd.Py());
    fSim_end_4mommenta.push_back(momentumEnd.Pz());
    fSim_end_4mommenta.push_back(momentumEnd.E());
    }*/
        // loop over every trajectory point, compare to geometry,
        size_t Ntraj = particleVec.NumberTrajectoryPoints();
        art::ServiceHandle<geo::Geometry const> geom;
        bool hasEntered = false;
        for (size_t ipt = 0; ipt < Ntraj; ++ipt)
        {
          // std::cout<<Ntraj<<std::endl;

          const TLorentzVector &pos = particleVec.Position(ipt);
          double localX = pos.X();
          double localY = pos.Y();
          double localZ = pos.Z();
          double X_MIN = -400.0, X_MAX = 400.0;
          double Y_MIN = -600.0, Y_MAX = 600.0;
          double Z_MIN = 0.0, Z_MAX = 1300.0;

          // std::cout << pos.X() << " ," << pos.Y() << "," << pos.Z() << std::endl;
          // std::cout << localX << " ," << localY << "," << localZ << std::endl;
          // std::cout << std::abs(centerX) << " ," << std::abs(centerY) << "," << std::abs(centerZ) << std::endl;
          bool inside =
              localX >= X_MIN && localX <= X_MAX &&
              localY >= Y_MIN && localY <= Y_MAX &&
              localZ >= Z_MIN && localZ <= Z_MAX;
          // std::abs(localX) <= tpc.HalfWidth() * 2 && std::abs(localY) <= tpc.HalfHeight() * 2 && std::abs(localZ) <= tpc.HalfLength() * 2;

          if (!hasEntered)
          {
            if (inside)
            {
              hasEntered = true;
              auto const &mom = particleVec.Momentum(ipt);
              double stepKE = (mom.E() - particleVec.Mass()); // in GeV
              std::cout << "Particle TRKID " << particleVec.TrackId() << ", PDG: " << particleVec.PdgCode()
                        << ", ENTERED at pt " << ipt << ", Position (" << pos.X() << "," << pos.Y() << "," << pos.Z() << "), Step Energy: " << stepKE << " GeV" << std ::endl;
              for (size_t i = 0; i + 1 < Ntraj; ++i)
              {
                auto const &a = particleVec.Position(i);
                auto const &b = particleVec.Position(i + 1);
                double dx = b.X() - a.X(), dy = b.Y() - a.Y(), dz = b.Z() - a.Z();
                double ds = std::sqrt(dx * dx + dy * dy + dz * dz);
                std::cout << "seg " << i << "->" << (i + 1) << "  ds=" << ds << " cm\n";
              }
              // fill out Energy(stepKE) histograms for protons, neutrons, electrons, muons, pions
              // go back to my branch
            }
            else
            {
              std::cout << " ! Particle TRKID " << particleVec.TrackId() << ", PDG: " << particleVec.PdgCode()
                        << ", Not ENTERED yet at pt " << ipt << ", Position (" << pos.X() << "," << pos.Y() << "," << pos.Z() << ") " << std ::endl;
            }
          }
          else
          {
            if (!inside)
            {
              // compute KE as before
              auto const &mom = particleVec.Momentum(ipt);
              double KE = mom.E() - particleVec.Mass();
              std::cout << "Particle " << particleVec.TrackId()
                        << " EXITED at pt " << ipt
                        << " with KE=" << KE << " GeV\n";
              break;
              // std::cout << pos.X() << " ," << pos.Y() << "," << pos.Z() << std::endl;
              // std::cout << localX << " ," << localY << "," << localZ << std::endl;
              // std::cout << "Particle: " << particleVec.TrackId() << ", PDG: " << particleVec.PdgCode() << ", Trajectory point: " << ipt << " Ntraj:" << Ntraj << std::endl;
            }
          }
        }
      }
      // End four-vector collection

      // Collecting all Daughters of Each primary

      std::vector<std::vector<const simb::MCParticle *>> DaughterpartVec;
      std::vector<const simb::MCParticle *> primary_vec;

      for (size_t i = 0; i < fSimP_TrackID_vec.size(); i++)
      {
        int currentMom = fSimP_Mom_vec[i];
        std::vector<const simb::MCParticle *> CurrentDaughters;
        CurrentDaughters.clear();
        const simb::MCParticle *currentpart = SimParticles[i];
        getDescendants(fSimP_TrackID_vec[i], fSimP_Mom_vec, fSimP_TrackID_vec, particleMap, CurrentDaughters);
        std::vector<Vertex> interactionVertices = clusterVertices(CurrentDaughters);
        // std::cout << "Number of Interaction Vertices for particle: " << fSimP_TrackID_vec[i] << " is: " << interactionVertices.size() << std::endl;
        for (const Vertex &vtx : interactionVertices)
        {
          fillInteractionTree(currentpart, vtx, particleMap, fInteractionTree, fInX, fInY, fInZ, fInT, fInPx, fInPy, fInPz, fInE, fInPDG, fOutX, fOutY, fOutZ, fOutT, fOutPx, fOutPy, fOutPz, fOutE, fOutPDG);
        }
        if (currentMom == 0)
        {
          int primary = fSimP_TrackID_vec[i];
          getDescendants(primary, fSimP_Mom_vec, fSimP_TrackID_vec, particleMap, CurrentDaughters);
          DaughterpartVec.push_back(CurrentDaughters);
          primary_vec.push_back(SimParticles[i]);
          int NHad = 0;
          double BindingE = 0.0;
          getHadronic02(SimParticles[i], SimParticles, NHad, BindingE);
          // std::cout << "Number Had interactions per primary: " << NHad << ", BindingE: " << BindingE << std::endl;
        }
      }

      // for(size_t n = 0; n < DaughterpartVec.size(); n++){
      //   int NHad = 0;
      //   int HadE = 0;
      //  getHadronicInformation(primary_vec[n], DaughterpartVec[n], NHad, HadE);
      //}

      // Calculate sim hadronic deposit energy
      //

      // Loop over the SimChannel objects in the event to look at the energy deposited by particle's track.
      for (auto const &channel : (*simChannelHandle))
      {

        // Get the numeric ID associated with this channel.
        // See methods at https://internal.dunescience.org/doxygen/SimChannel_8h_source.html
        auto const channelNumber = channel.Channel();

        // Each channel has a map inside it that connects a time slice to energy deposits in the detector.
        // The full type of this map is std::map<unsigned short, std::vector<sim::IDE>>; we'll use "auto" here
        auto const &timeSlices = channel.TDCIDEMap();
        for (auto const &timeSlice : timeSlices)
        {

          // For the timeSlices map, the 'first' is a time slice number; The 'second' is a vector of IDE objects.
          auto const &energyDeposits = timeSlice.second;

          for (auto const &energyDeposit : energyDeposits)
          {

            // Method b: First check if it's on collection plane
            std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
            if (Wires[0].planeID().Plane == 0)
            {

              // All EM shower are treated as secondary interactions, and their particles are not saved in the MC particle list
              // Still do the search, but now only for primary lepton (particleMap trkID is always positive)
              // Also search for EM shower particles from primary lepton, these deposits has trkID that's negative of the primary lepton trkID
              auto search = particleMap.find(abs(energyDeposit.trackID));

              // std::cout << "Time Slice Number: " << timeSlice.first << "Energy Deposit TrackID: " << energyDeposit.trackID << "Energy Deposit Energy: "<< energyDeposit.energy << std::endl;

              if (search != particleMap.end())
              { // found match in map

                const simb::MCParticle &particle = *((*search).second);

                // std::cout << particle.PdgCode() << std::endl;

                // if the energy deposit is from primary lepton,
                // or its ancestor mother particle is the primary lepton (e.g., from muon decays)
                if ((particle.Process() == "primary" && abs(particle.PdgCode()) == 13) || IsAncestorMotherPrimaryLep(particle, primarylep_trkID, particleMap))
                {
                  fSim_mu_Edep_b2 += energyDeposit.energy;
                  // now continue to the next energy deposit
                  // continue here to avoid counting into fSim_hadronic_Edep_b2
                  continue;
                } // end lepton deposited energy

                // if ( particle.PdgCode() == 22 && particle.Mother() == pi0_trkID ) {
                // std::cout << "EDep MeV: "<< energyDeposit.energy << " from gamma: Mother trkid: " << pi0_trkID << ", E: " << particle.E() << ", mass: " << particle.Mass() << std::endl;
                //}

                // if the energy deposit is from neutron
                // or its ancestor mother particle is the neutron
                if (particle.PdgCode() == 2112 || IsAncestorMotherNeutron(particle, neutron_trkID, particleMap))
                {
                  fSim_n_Edep_b2 += energyDeposit.energy;
                } // std::cout << "fire n! " << std::endl;  // end neutron deposited energy
                // if the energy deposit is from proton
                // or its ancestor mother particle is the proton
                else if (particle.PdgCode() == 2212 || IsAncestorMotherProton(particle, proton_trkID, particleMap))
                {
                  fSim_p_Edep_b2 += energyDeposit.energy;
                } // std::cout << "fire p! " << std::endl;  // end proton deposited energy
                // if the energy deposit is from primary pip
                // or its ancestor mother particle is the pip
                else if (particle.PdgCode() == 211 || IsAncestorMotherPip(particle, pip_trkID, particleMap))
                {
                  fSim_pip_Edep_b2 += energyDeposit.energy;
                } // std::cout << "fire pip! pdg: "<< particle.PdgCode() << ", pip_trkID: " << pip_trkID << ", IsAncestorMotherPip: " << IsAncestorMotherPip(particle, pip_trkID, particleMap) << std::endl; } // end pip deposited energy
                // if the energy deposit is from primary pim
                // or its ancestor mother particle is the pim (pi0 decays into two photons)
                else if (particle.PdgCode() == -211 || IsAncestorMotherPim(particle, pim_trkID, particleMap))
                {
                  fSim_pim_Edep_b2 += energyDeposit.energy;
                } // std::cout << "fire pim! " << std::endl; } // end pim deposited energy
                // if the energy deposit is from primary pi0
                // or its ancestor mother particle is the pi0 (pi0 decays into two photons)
                else if (particle.PdgCode() == 111 || IsAncestorMotherPi0(particle, pi0_trkID, particleMap))
                {
                  fSim_pi0_Edep_b2 += energyDeposit.energy;
                } // std::cout << "fire pi0! " << std::endl; } // end pi0 deposited energy
                else if (particle.PdgCode() == 321 || particle.PdgCode() == -321 || particle.PdgCode() == 311 || particle.PdgCode() == -311 || particle.PdgCode() == 130 || particle.PdgCode() == 310 || particle.PdgCode() == 22 || (particle.PdgCode() >= 100 && particle.PdgCode() <= -9999) || (particle.PdgCode() >= -9999 && particle.PdgCode() <= -100)) // eOther which includes: kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)
                {
                  fSim_Other_Edep_b2 += energyDeposit.energy;
                  // std::cout << "fire Other! " << std::endl;
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
              auto exist = EDepMap.find(EDepTrackID);
              // if can't find, store it and add the edep
              if (exist == EDepMap.end())
              {
                EDep_TrackID_vec.push_back(EDepTrackID); // negative trackID exists
                EDepMap[EDepTrackID] = energyDeposit.energy;
              }
              else
              { // find the same track id
                EDepMap[EDepTrackID] += energyDeposit.energy;
              }

            } // end plane == 0
          } // end energy deposit loop
        } // end For each time slice
      } // end For each SimChannel

      fSim_n_hadronic_Edep_b = fSim_hadronic_hit_x_b.size();

      // Print out EDepMap
      // if ( false ) std::cout << "fGen_numu_E: "<< fGen_numu_E << ", fSim_mu_Edep_b2_debug: " << fSim_mu_Edep_b2_debug<< ", fSim_hadronic_Edep_b2_debug: " << fSim_hadronic_Edep_b2_debug << ", Tot had E track id: " << EDep_TrackID_vec.size() << std::endl;
      if (false)
      {
        for (long unsigned int i = 0; i < EDep_TrackID_vec.size(); i++)
        {
          std::cout << "Evt track id: " << EDep_TrackID_vec.at(i) << std::endl;
        }
        std::map<int, double>::iterator it;
        std::cout << "TrackID" << " | " << "Tot EDep" << std::endl;
        for (it = EDepMap.begin(); it != EDepMap.end(); it++)
          std::cout << "    " << it->first << " | " << it->second << std::endl;
      }

      // In general, objects in the LArSoft reconstruction chain are linked using the art::Assns class:
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artAssns>
      // The following statement will find the simb::MCTruth associated with the simb::MCParticle
      const art::FindManyP<simb::MCTruth> findManyTruth(particleHandle, event, fSimulationProducerLabel);

      if (!findManyTruth.isValid())
      {
        std::cout << "findManyTruth simb::MCTruth for simb::MCParticle failed!" << std::endl;
      }

      size_t particle_index = 0; // only look at first particle in particleHandle's vector.
      auto const &truth = findManyTruth.at(particle_index);

      // Make sure there's no problem.
      if (truth.empty())
      {
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
namespace
{

  double DetectorDiagonal(geo::GeometryCore const &geom)
  {
    const double length = geom.DetLength();
    const double width = 2. * geom.DetHalfWidth();
    const double height = 2. * geom.DetHalfHeight();

    return std::sqrt(cet::sum_of_squares(length, width, height));
  }

  // bool MomentumOrderMCParticle(const simb::MCParticle* p1, const simb::MCParticle* p2) {
  //   return ( p1->P(0) > p2->P(0) );
  // }

  // If this returns true, then the energy deposit is associated with primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle &p1, int primarylep_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    // Immediate mother is the primary lep
    if (MothertrkID == primarylep_trkID)
      return true;
    // Immediate mother is not primary lep, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPrimaryLep(tmp_mother, primarylep_trkID, particleMap);
    }
  } // end GetAncestorMotherLeptonTrkID

  // If this returns true, then the energy deposit is associated with neutron
  bool IsAncestorMotherNeutron(const simb::MCParticle &p1, std::vector<int> neutron_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the neutrons in the event
    for (long unsigned int i = 0; i < neutron_trkID.size(); i++)
    {
      if (MothertrkID == neutron_trkID.at(i))
        MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true)
      return true;

    // Immediate mother is not neutron, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherNeutron(tmp_mother, neutron_trkID, particleMap);
    }
  } // end GetAncestorMotherNeutronTrkID

  // If this returns true, then the energy deposit is associated with proton
  bool IsAncestorMotherProton(const simb::MCParticle &p1, std::vector<int> proton_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the protons in the event
    for (long unsigned int i = 0; i < proton_trkID.size(); i++)
    {
      if (MothertrkID == proton_trkID.at(i))
        MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true)
      return true;

    // Immediate mother is not proton, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherProton(tmp_mother, proton_trkID, particleMap);
    }
  } // end GetAncestorMotherProtonTrkID

  // If this returns true, then the energy deposit is associated with pion+
  bool IsAncestorMotherPip(const simb::MCParticle &p1, std::vector<int> pip_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pips in the event
    for (long unsigned int i = 0; i < pip_trkID.size(); i++)
    {
      if (MothertrkID == pip_trkID.at(i))
        MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true)
      return true;

    // Immediate mother is not pip, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPip(tmp_mother, pip_trkID, particleMap);
    }
  } // end GetAncestorMotherPipTrkID

  // If this returns true, then the energy deposit is associated with pion+
  bool IsAncestorMotherPim(const simb::MCParticle &p1, std::vector<int> pim_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pims in the event
    for (long unsigned int i = 0; i < pim_trkID.size(); i++)
    {
      if (MothertrkID == pim_trkID.at(i))
        MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true)
      return true;

    // Immediate mother is not pim, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPim(tmp_mother, pim_trkID, particleMap);
    }
  } // end GetAncestorMotherPimTrkID

  // If this returns true, then the energy deposit is associated with pion-
  bool IsAncestorMotherPi0(const simb::MCParticle &p1, std::vector<int> pi0_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    // Immediate mother is one of the pi0s in the event
    for (long unsigned int i = 0; i < pi0_trkID.size(); i++)
    {
      if (MothertrkID == pi0_trkID.at(i))
        MatchMultipleTrkID = true;
    }

    if (MatchMultipleTrkID == true)
      return true;

    // Immediate mother is not pi0, but other primary particles from genie
    else if (MothertrkID == 0)
      return false;
    // Keep looking upstream, find it in particleMap
    else
    {
      auto tmp_search = particleMap.find(MothertrkID); // this search must be found, can't be null
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPi0(tmp_mother, pi0_trkID, particleMap);
    }

  } // end GetAncestorMotherPi0TrkID

  // void getHadronicInformation(const simb::MCParticle* primary, const std::vector<const simb::MCParticle*>& daughters, int& NHad, double& BindingE){
  //   int pNTP = primary->NumberTrajectoryPoints();
  //   int pLast = pNTP - 1;
  //   for(size_t k = 0; k < daughters.size(); k++){
  //     //int dNTP = daughters[k]->NumberTrajectoryPoints();
  //     //int dLast = dNTP - 1;
  //     const TLorentzVector& daughterstart = daughters[k]->Position(0);
  //     const TLorentzVector& Edaughterstart = daughters[k]->Momentum(0);
  //     std::vector<float> X;
  //     std::vector<float> Y;
  //     std::vector<float> Z;
  //     std::vector<float> T;
  //     for(int l = 0; l < pLast; l++){
  //       const TLorentzVector& pripos = primary->Position(l);
  //       float epsilon = 0.01;
  //       double Ein = 0;
  //       double Eout = 0;
  //       if(abs(pripos.X() - daughterstart.X()) < epsilon && abs(pripos.Y() - daughterstart.Y()) < epsilon && abs(pripos.Z() - daughterstart.Z()) < epsilon){
  //         if(abs(primary->PdgCode()) == 211){
  //           Ein = primary->E(l);
  //         }
  //         else{
  //           Ein = primary->E(l) - primary->Mass();
  //         }
  //         std::cout << "PDG of primary: " << primary->PdgCode() << std::endl;
  //         std::cout << "PDG of daughter: " << daughters[k]->PdgCode() << std::endl;
  //         if(abs(daughters[k]->PdgCode()) == 211){
  //          Eout = Edaughterstart.E();
  //         }
  //         else{
  //           Eout = Edaughterstart.E() - daughters[k]->Mass();
  //         }
  //         std::cout << "Ein: " << Ein << std::endl;
  //         std::cout << "Eout: " << Eout << std::endl;
  //         double currentBindingE = Ein - Eout;
  //         std::cout << "Current Binding Energy: " << currentBindingE << std::endl;
  //         BindingE += currentBindingE;
  //         std::cout << "BindingE (inside func) : " << BindingE << std::endl;
  //         X.push_back(daughterstart.X());
  //         Y.push_back(daughterstart.Y());
  //         Z.push_back(daughterstart.Z());
  //         T.push_back(daughterstart.T());
  //         int Xsize = X.size();
  //         int Xlast = Xsize - 1;
  //         int iterator = 0;
  //         for(size_t m = 0; m < X.size(); m++){
  //           if(abs(X[Xlast] - X[m]) < epsilon && abs(Y[Xlast] - Y[m]) < epsilon && abs(Z[Xlast] - Z[m]) < epsilon && abs(T[Xlast] - T[m]) < epsilon) iterator = iterator +1;
  //           std::cout << "iterator: " << iterator << std::endl;
  //         }
  //         if(iterator == 1) NHad = NHad +1;
  //       }
  //     }
  //     std::cout << "NHad (inside): " << NHad << std::endl;
  //   }
  // }

  void fillInteractionTree(const simb::MCParticle* incoming,
    const Vertex& vertex,
    const std::map<int, const simb::MCParticle*>& particleMap,
    TTree* fInteractionTree,
    float& fInX, float& fInY, float& fInZ, float& fInT,
    float& fInPx, float& fInPy, float& fInPz, float& fInE, int& fInPDG,
    std::string& fInProcess,
    std::vector<float>& fOutX, std::vector<float>& fOutY,
    std::vector<float>& fOutZ, std::vector<float>& fOutT,
    std::vector<float>& fOutPx, std::vector<float>& fOutPy,
    std::vector<float>& fOutPz, std::vector<float>& fOutE,
    std::vector<int>& fOutPDG, std::vector<std::string>& fOutProcess) {

  // Clear outgoing particle containers
  fOutX.clear(); fOutY.clear(); fOutZ.clear(); fOutT.clear();
  fOutPx.clear(); fOutPy.clear(); fOutPz.clear(); fOutE.clear();
  fOutPDG.clear(); fOutProcess.clear();

  // Basic incoming particle info
  fInX = vertex.x; 
  fInY = vertex.y; 
  fInZ = vertex.z;
  fInT = vertex.t;
  fInPDG = incoming->PdgCode();
  fInProcess = incoming->EndProcess();

  int incomingID = incoming->TrackId();
  double minDist = 1e10;
  TLorentzVector bestMom;
  bool dies = false;
  TLorentzVector nextpos;
  TLorentzVector nextmom;

  for (unsigned int i = 0; i <= incoming->NumberTrajectoryPoints(); i++) {
    TLorentzVector pos = incoming->Position(i);
    double dist = std::hypot(pos.X() - vertex.x, pos.Y() - vertex.y, pos.Z() - vertex.z);
    if (dist < minDist) {
      minDist = dist;
      bestMom = incoming->Momentum(i);
      if (i == incoming->NumberTrajectoryPoints()) {
      dies = true;
      }
      if (i != incoming-> NumberTrajectoryPoints()){
        dies = false;
      nextpos = incoming->Position(i+1);
      nextmom = incoming->Momentum(i+1);
      }
    }
  }
  


  // --- Group daughters by production time ---
  const double timeEpsilon = 1e-3; // ns, small tolerance for clustering
  std::map<double, std::vector<const simb::MCParticle*>> timeGroups;

  for (const simb::MCParticle* daughter : vertex.daughters) {
    if (daughter->TrackId() == incomingID) continue;        // skip self
    if (daughter->Mother() != incoming->TrackId()) continue; // skip indirect descendants

    double t = daughter->Position(0).T();
    bool added = false;

    // Look for an existing time bin within tolerance
    for (auto& kv : timeGroups) {
      if (std::fabs(kv.first - t) < timeEpsilon) {
        kv.second.push_back(daughter);
        added = true;
        break;
      }
    }

    // If none exists, create a new time bin
    if (!added) {
      timeGroups[t].push_back(daughter);
    }
  }

  fInPx = bestMom.Px();
  fInPy = bestMom.Py();
  fInPz = bestMom.Pz();
  fInE  = bestMom.E();

  if(!dies){
    std::cout << "Dies false" << std::endl;
  }


  // --- Fill one TTree entry per time group ---
  for (const auto& kv : timeGroups) {
    // Clear outgoing vectors for this time cluster
    fOutX.clear(); fOutY.clear(); fOutZ.clear(); fOutT.clear();
    fOutPx.clear(); fOutPy.clear(); fOutPz.clear(); fOutE.clear();
    fOutPDG.clear(); fOutProcess.clear();

    for (const simb::MCParticle* daughter : kv.second) {
      const TLorentzVector& pos = daughter->Position(0);
      const TLorentzVector& mom = daughter->Momentum(0);

      fOutX.push_back(pos.X());
      fOutY.push_back(pos.Y());
      fOutZ.push_back(pos.Z());
      fOutT.push_back(pos.T());

      fOutPx.push_back(mom.Px());
      fOutPy.push_back(mom.Py());
      fOutPz.push_back(mom.Pz());
      fOutE.push_back(mom.E());
      fOutPDG.push_back(daughter->PdgCode());
      fOutProcess.push_back(daughter->EndProcess());
    }
    std::cout << "Incoming particle process: " << fInProcess << std::endl;
    if(!fOutT.empty()){
      if(!dies){
        std::cout << "Incoming scattered, adding incoming particle to outgoing list to preserve energy/momentum conservation" << std::endl;
        fOutX.push_back(nextpos.X());
        fOutY.push_back(nextpos.Y());
        fOutZ.push_back(nextpos.Z());
        fOutT.push_back(nextpos.T());
        fOutPx.push_back(nextmom.Px());
        fOutPy.push_back(nextmom.Py());
        fOutPz.push_back(nextmom.Pz());
        fOutE.push_back(nextmom.E());
        fOutPDG.push_back(incoming->PdgCode());
        fOutProcess.push_back("nucleonScat");
      }

      // if(dies && fInProcess == "Decay" && std::abs(fOutT.back() - fInT) < timeEpsilon){
      //   std::cout << "Incoming particle decayed at rest, adding artificial daughter to preserve energy/momentum conservation" << std::endl;
      //  fOutX.push_back(fOutX.back());
      //  fOutY.push_back(fOutY.back());
      //  fOutZ.push_back(fOutZ.back());
      //  fOutT.push_back(fOutT.back());

      //  fOutPx.push_back(0.0);
      //  fOutPy.push_back(0.0);
      //  fOutPz.push_back(0.0);
      //  fOutE.push_back(incoming->Mass());
      //  fOutPDG.push_back(incoming->PdgCode());
      //  fOutProcess.push_back("artificialAtRest");
      // }
      if(dies && fInProcess == "Decay" && std::abs(fInT - fOutT.back()) > 3){
        std::cout << "Incoming particle decayed at rest with time delay, setting incoming momentum to 0 to preserve energy/momentum conservation" << std::endl;
        fInPx = 0.0;
        fInPy = 0.0;
        fInPz = 0.0;
        fInE = incoming->Mass();
        fInProcess = "artificialAtRest";
      }
    }
    // Only fill if we have outgoing particles for this time group
    if (!fOutX.empty()) {
      fInteractionTree->Fill();
    }
  }
}


  // std::vector<primaryVertex> clusterPrimaryVertices(const simb::MCParticle* incoming, const std::vector<const simb::MCParticle*>& daughters){
  //   float epsilon = 0.01;
  //   float tepsilon = 1e-3;
  //   std::vector<primaryVertex> vtxs;

  //   for (const simb::MCParticle* d : daughters){
  //     const TLorentzVector& pos = d->Position(0);
  //     float x = pos.X(), y = pos.Y(), z = pos.Z(), t = pos.T();
  //     bool found = false;

  //     for (primaryVertex& v : vtxs) {
  //       if (std::abs(v.x - x) < epsilon && std::abs(v.y - y) < epsilon && std::abs(v.z - z) < epsilon && std::abs(v.t - t) < tepsilon) {
  //         v.daughters.push_back(d);
  //         found = true;
  //         break;
  //       }
  //     }

  //     if (!found) {
  //       primaryVertex vtx = {x, y, z, t, incoming, {d}};
  //       vtxs.push_back(vtx);
  //     }
  //   }
  //     return vtxs;
  // }

  std::vector<Vertex> clusterVertices(const std::vector<const simb::MCParticle*>& daughters){
    std::vector<Vertex> vertices;

    float epsilon = 0.01; 
    //float tepsilon = 1e-3;
    
    for (const simb::MCParticle* d : daughters){
      const TLorentzVector& pos = d->Position(0);
      float x = pos.X(), y = pos.Y(), z = pos.Z(), t = pos.T();
      bool found = false;
      for (Vertex& v : vertices) {
        if (std::abs(v.x - x) < epsilon && std::abs(v.y - y) < epsilon && std::abs(v.z - z) < epsilon) {
          v.daughters.push_back(d);
          found = true;
          break;
        }
      }
      if (!found){
        Vertex vert = {x, y, z, t, {d}};
        vertices.push_back(vert);
      }
  }
  return vertices;
}

double getPrimaryKE(const simb::MCParticle* primary, double x, double y, double z){
  double minDist = 1e10;
  int closestDist = 0;

  for(unsigned int n = 0; n < primary->NumberTrajectoryPoints(); ++n){
    const TLorentzVector& position = primary->Position(n);
    double dist = std::sqrt(std::pow(position.X() - x, 2) + std::pow(position.Y() - y, 2) + std::pow(position.Z() - z, 2));
    if (dist < minDist) {
      minDist = dist;
      closestDist = n;
    }
  }
  const TLorentzVector& ClosestMom = primary->Momentum(closestDist);
  return ClosestMom.E() - primary->Mass();
}

void getHadronic02(const simb::MCParticle* particle, const std::vector<const simb::MCParticle*>& allPart, int& NHad, int& Nintlow, double& totalBindingE){
  std::vector<const simb::MCParticle*> daughters;
  TLorentzVector currentPos = particle->Position(0);

  for(const simb::MCParticle* p : allPart){
    if(p->Mother() == particle->TrackId()){
      daughters.push_back(p);
    }
  }
  
  if(!daughters.empty()){
    std::vector<Vertex> vertices = clusterVertices(daughters);
    float BindingE = 0.0;

    for(const auto& vertex : vertices){
      double Ein = getPrimaryKE(particle, vertex.x, vertex.y, vertex.z);
      double Eout = 0.0;

      for(const simb::MCParticle* daughter : vertex.daughters){
        if(daughter->PdgCode() == 211){ // Check if daughter is a pion
          Eout += daughter->Momentum(0).E();
        } else {
          Eout += daughter->Momentum(0).E() - daughter->Mass(); // For other particles, subtract mass
        }
      }
      BindingE = Ein - Eout;
      if(BindingE > 0.001){
        totalBindingE += BindingE;
        NHad++;
      }
      if(BindingE < 0.001 && BindingE > -0.001){
        Nintlow++;
      }
    }
  }
  for(const simb::MCParticle* daughter : daughters){
    getHadronic02(daughter, allPart, NHad, Nintlow, totalBindingE);
  }
}

  void getDescendants(int motherID, const std::vector<int> &momVec, const std::vector<int> &TrkIDvec, const std::map<int, const simb::MCParticle *> &particleMap, std::vector<const simb::MCParticle *> &primaryDaughters)
  {
    for (size_t j = 0; j < TrkIDvec.size(); j++)
    {
      if (momVec[j] == motherID)
      {
        int daughterID = TrkIDvec[j];

        auto it = particleMap.find(daughterID);
        if (it != particleMap.end())
        {
          primaryDaughters.push_back(it->second);

          getDescendants(daughterID, momVec, TrkIDvec, particleMap, primaryDaughters);
        }
      }
    }
  }

} // local namespace
