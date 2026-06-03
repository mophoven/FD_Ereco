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

  //void getHadronicInformation(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&, int, double);

  double getMassFromPDG(int);

  // Particle role classification for energy budget accounting
  enum class ParticleRole {
    VisibleCharged,   // p, pi+-, mu, e — deposit in TPC
    NeutralEscaping,  // n, gamma — leave no track
    Pi0,              // pi0 — decays to 2 gammas, handle separately
    NuclearRecoil,    // PDG >= 1e9 — nucleus absorbs energy
    Neutrino,         // leaves entirely
    Unknown           // unrecognized PDG
  };

  // Classify a PDG code into a ParticleRole for energy budget accounting
  ParticleRole classifyParticle(int pdg);

  // Interpolate the incoming particle's 4-momentum at a vertex position
  TLorentzVector getMomentumAtVertex(const simb::MCParticle* p,
                                      float vx, float vy, float vz,
                                      bool& dies);

  // Fill one TTree entry per interaction vertex with full 4-momentum energy budget
  void fillInteractionTree(
      const simb::MCParticle*                    incoming,
      const Vertex&                              vertex,
      const std::map<int, const simb::MCParticle*>& particleMap,
      TTree*                                     tree,
      float& fVtx_x,  float& fVtx_y,  float& fVtx_z,  float& fVtx_t,
      int&   fIn_PDG, int&   fIn_TrackID,
      float& fIn_Px,  float& fIn_Py,  float& fIn_Pz,  float& fIn_E,
      float& fIn_Mass,
      std::string& fIn_Process,
      bool&  fIn_Died,
      float& fE_visible,
      float& fE_neutral,
      float& fE_pi0,
      float& fE_recoil,
      float& fE_neutrino,
      float& fE_binding,
      float& fQValue,
      float& fmichDifference,
      std::vector<int>&   fOut_PDG,
      std::vector<float>& fOut_Px,
      std::vector<float>& fOut_Py,
      std::vector<float>& fOut_Pz,
      std::vector<float>& fOut_E,
      std::vector<float>& fOut_Mass,
      std::vector<int>&   fOut_Role,
      std::vector<int>&   fOut_TrackID,
      int& fEvent);

  std::vector<Vertex> clusterVertices(const std::vector<const simb::MCParticle*>&);
  double getPrimaryKE(const simb::MCParticle*, double, double, double);
  void getHadronic02(const simb::MCParticle*, const std::vector<const simb::MCParticle*>&, int&, double&);
  void getDescendants(int, const std::vector<int> &, const std::vector<int> &, const std::map<int, const simb::MCParticle *> &, std::vector<const simb::MCParticle *> &);
  void getAncestors(const simb::MCParticle *, std::vector<int> &, const std::map<int, const simb::MCParticle *> &);
  void ReportFirstExitRootOnly(const simb::MCParticle &, const std::map<int, const simb::MCParticle *> &, double &);

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

        fhicl::Atom<art::InputTag> SimChannelLabel{
            Name("SimChannelLabel"),
            Comment("tag of the SimChannel data product, e.g. tpcrawdecoder:simpleSC")};

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
      art::InputTag fSimChannelLabel;         // The name of the producer that made SimChannels e.g. tpcrawdecoder:simpleSC

      // --- Interaction tree (one entry per hadronic vertex) ---
      TTree* fInteractionTree;

      // Vertex location
      float fVtx_x, fVtx_y, fVtx_z, fVtx_t;

      // Incoming particle at vertex
      int         fIn_PDG;
      int         fIn_TrackID;
      float       fIn_Px, fIn_Py, fIn_Pz, fIn_E;
      float       fIn_Mass;
      std::string fIn_Process;
      bool        fIn_Died;

      // Per-category energy sums at this vertex (GeV)
      float fE_visible;    // KE of charged visible particles (p, pi+-, mu, e)
      float fE_neutral;    // E  of neutral escaping particles (n, gamma)
      float fE_pi0;        // E  of pi0s (partially visible via EM showers)
      float fE_recoil;     // KE of nuclear recoil fragments
      float fE_neutrino;   // E  of neutrinos
      float fE_binding;    // residual: E_in - all of the above - Q-value
      float fQValue;       // nuclear Q-value: sum(M_nuclear_out) - M_nuclear_in

      // Outgoing particle list (one entry per daughter)
      std::vector<int>   fOut_PDG;
      std::vector<float> fOut_Px, fOut_Py, fOut_Pz, fOut_E;
      std::vector<float> fOut_Mass;
      std::vector<int>   fOut_Role;     // ParticleRole cast to int
      std::vector<int>   fOut_TrackID;

      float fmichDifference; // Michel electron energy difference (mu -> e)

      // The n-tuple to create
      TTree *fNtuple;

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
      std::vector<simb::MCTrajectory> fSimP_Traj_vec; // saving trajectory points (from V1)

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
      double fSim_nuclei_Edep_b2; // [MeV] Energy Deposit of nuclei recoil

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

      double ftotalExited; // Total KE of primary particles exiting the detector (from V2)

    }; // class MyEnergyAnalysis

    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor

    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const &config)
        : EDAnalyzer(config),
          fGenieGenModuleLabel(config().GenieGenModuleLabel()),
          fSimulationProducerLabel(config().SimulationLabel()),
          fSimChannelLabel(config().SimChannelLabel())
    {
      // Get a pointer to the geometry service provider.
      fGeometryService = &*art::ServiceHandle<geo::Geometry>();

      consumes<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<sim::SimChannel>>(fSimChannelLabel);
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
      fInteractionTree = tfs->make<TTree>("HadronicTree", "Hadronic Interaction Vertices");

      // Vertex location
      fInteractionTree->Branch("Vtx_x",    &fVtx_x,    "Vtx_x/F");
      fInteractionTree->Branch("Vtx_y",    &fVtx_y,    "Vtx_y/F");
      fInteractionTree->Branch("Vtx_z",    &fVtx_z,    "Vtx_z/F");
      fInteractionTree->Branch("Vtx_t",    &fVtx_t,    "Vtx_t/F");

      // Incoming particle at vertex
      fInteractionTree->Branch("In_PDG",     &fIn_PDG,     "In_PDG/I");
      fInteractionTree->Branch("In_TrackID", &fIn_TrackID, "In_TrackID/I");
      fInteractionTree->Branch("In_Px",      &fIn_Px,      "In_Px/F");
      fInteractionTree->Branch("In_Py",      &fIn_Py,      "In_Py/F");
      fInteractionTree->Branch("In_Pz",      &fIn_Pz,      "In_Pz/F");
      fInteractionTree->Branch("In_E",       &fIn_E,       "In_E/F");
      fInteractionTree->Branch("In_Mass",    &fIn_Mass,    "In_Mass/F");
      fInteractionTree->Branch("In_Process", &fIn_Process);
      fInteractionTree->Branch("In_Died",    &fIn_Died,    "In_Died/O");

      // Per-category energy sums
      fInteractionTree->Branch("E_visible",  &fE_visible,  "E_visible/F");
      fInteractionTree->Branch("E_neutral",  &fE_neutral,  "E_neutral/F");
      fInteractionTree->Branch("E_pi0",      &fE_pi0,      "E_pi0/F");
      fInteractionTree->Branch("E_recoil",   &fE_recoil,   "E_recoil/F");
      fInteractionTree->Branch("E_neutrino", &fE_neutrino, "E_neutrino/F");
      fInteractionTree->Branch("E_binding",  &fE_binding,  "E_binding/F");
      fInteractionTree->Branch("QValue",     &fQValue,     "QValue/F");

      // Outgoing particle list
      fInteractionTree->Branch("Out_PDG",     &fOut_PDG);
      fInteractionTree->Branch("Out_Px",      &fOut_Px);
      fInteractionTree->Branch("Out_Py",      &fOut_Py);
      fInteractionTree->Branch("Out_Pz",      &fOut_Pz);
      fInteractionTree->Branch("Out_E",       &fOut_E);
      fInteractionTree->Branch("Out_Mass",    &fOut_Mass);
      fInteractionTree->Branch("Out_Role",    &fOut_Role);
      fInteractionTree->Branch("Out_TrackID", &fOut_TrackID);

      // Event context
      fInteractionTree->Branch("michDifference", &fmichDifference, "michDifference/F");
      fInteractionTree->Branch("Event",          &fEvent,          "Event/I");


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
      fNtuple->Branch("SimP_Traj_vec", &fSimP_Traj_vec); // from V1
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
      fNtuple->Branch("Sim_nuclei_Edep_b2", &fSim_nuclei_Edep_b2, "Sim_nuclei_Edep_b2/D");

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
      fNtuple->Branch("totalExited", &ftotalExited); // from V2
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

      ftotalExited = 0.; // from V2

      // Reset interaction tree per-event variables
      fVtx_x = fVtx_y = fVtx_z = fVtx_t = 0.f;
      fIn_PDG = fIn_TrackID = 0;
      fIn_Px = fIn_Py = fIn_Pz = fIn_E = fIn_Mass = 0.f;
      fIn_Process = "";
      fIn_Died = false;
      fE_visible = fE_neutral = fE_pi0 = fE_recoil = fE_neutrino = 0.f;
      fE_binding = fQValue = 0.f;
      fmichDifference = 0.f;
      fOut_PDG.clear();
      fOut_Px.clear(); fOut_Py.clear(); fOut_Pz.clear(); fOut_E.clear();
      fOut_Mass.clear(); fOut_Role.clear(); fOut_TrackID.clear();

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
      fSimP_Traj_vec.clear(); // from V1
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
      // in which case you'd run GENIE multiple times and have one MCTruth per interaction.
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
      auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);

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
        fSimP_Traj_vec.push_back(particle.Trajectory()); // from V1
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
        }

        // Take note of proton trackID
        if (fSimPDG == 2212)
        {
          proton_trkID.push_back(fSimTrackID);
        }

        // Take note of pip track ID
        if (fSimPDG == 211)
        {
          pip_trkID.push_back(fSimTrackID);
        }

        // Take note of primary pim track ID
        if (fSimPDG == -211)
        {
          pim_trkID.push_back(fSimTrackID);
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
        }

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

      // Compute total KE escaping the detector for root-level primary particles (from V2)
      double KE = 0;
      for (int i = 0; i < fSim_nParticles; i++)
      {
        const simb::MCParticle &particleVec = *(SimParticles[i]);
        if (particleVec.Process() == "primary")
        {
          ReportFirstExitRootOnly(particleVec, particleMap, KE);
          ftotalExited += KE;
        }
      }

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
        for (const Vertex &vtx : interactionVertices)
        {
          fillInteractionTree(
              currentpart, vtx, particleMap, fInteractionTree,
              fVtx_x, fVtx_y, fVtx_z, fVtx_t,
              fIn_PDG, fIn_TrackID,
              fIn_Px, fIn_Py, fIn_Pz, fIn_E, fIn_Mass,
              fIn_Process, fIn_Died,
              fE_visible, fE_neutral, fE_pi0, fE_recoil, fE_neutrino,
              fE_binding, fQValue, fmichDifference,
              fOut_PDG, fOut_Px, fOut_Py, fOut_Pz, fOut_E,
              fOut_Mass, fOut_Role, fOut_TrackID,
              fEvent);
        }
        if (currentMom == 0)
        {
          int primary = fSimP_TrackID_vec[i];
          getDescendants(primary, fSimP_Mom_vec, fSimP_TrackID_vec, particleMap, CurrentDaughters);
          DaughterpartVec.push_back(CurrentDaughters);
          primary_vec.push_back(SimParticles[i]);
          int NHad = 0;
          double BindingE;
          getHadronic02(SimParticles[i], SimParticles, NHad, BindingE);
        }
      }

      // Calculate sim hadronic deposit energy
      //

      // Loop over the SimChannel objects in the event to look at the energy deposited by particle's track.
      for (auto const &channel : (*simChannelHandle))
      {

        auto const channelNumber = channel.Channel();

        auto const &timeSlices = channel.TDCIDEMap();
        for (auto const &timeSlice : timeSlices)
        {

          auto const &energyDeposits = timeSlice.second;

          for (auto const &energyDeposit : energyDeposits)
          {

            // Method b: First check if it's on collection plane
            std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
            if (Wires[0].planeID().Plane == 0)
            {

              auto search = particleMap.find(abs(energyDeposit.trackID));

              if (search != particleMap.end())
              { // found match in map

                const simb::MCParticle &particle = *((*search).second);

                if ((particle.Process() == "primary" && abs(particle.PdgCode()) == 13) || IsAncestorMotherPrimaryLep(particle, primarylep_trkID, particleMap))
                {
                  fSim_mu_Edep_b2 += energyDeposit.energy;
                  continue;
                } // end lepton deposited energy

                if (particle.PdgCode() == 2112 || IsAncestorMotherNeutron(particle, neutron_trkID, particleMap))
                {
                  fSim_n_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() == 2212 || IsAncestorMotherProton(particle, proton_trkID, particleMap))
                {
                  fSim_p_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() == 211 || IsAncestorMotherPip(particle, pip_trkID, particleMap))
                {
                  fSim_pip_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() == -211 || IsAncestorMotherPim(particle, pim_trkID, particleMap))
                {
                  fSim_pim_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() == 111 || IsAncestorMotherPi0(particle, pi0_trkID, particleMap))
                {
                  fSim_pi0_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() == 321 || particle.PdgCode() == -321 || particle.PdgCode() == 311 || particle.PdgCode() == -311 || particle.PdgCode() == 130 || particle.PdgCode() == 310 || particle.PdgCode() == 22 || (particle.PdgCode() >= 100 && particle.PdgCode() <= -9999) || (particle.PdgCode() >= -9999 && particle.PdgCode() <= -100)) // eOther
                {
                  fSim_Other_Edep_b2 += energyDeposit.energy;
                }
                else if (particle.PdgCode() >= 1000000000 && particle.PdgCode() <= 9999999999) // nucleus
                {
                  fSim_nuclei_Edep_b2 += energyDeposit.energy;
                }
              } // end found match

              fSim_hadronic_Edep_b2 += energyDeposit.energy;
              fSim_hadronic_hit_x_b.push_back(energyDeposit.x);
              fSim_hadronic_hit_y_b.push_back(energyDeposit.y);
              fSim_hadronic_hit_z_b.push_back(energyDeposit.z);
              fSim_hadronic_hit_Edep_b2.push_back(energyDeposit.energy);

              EDepTrackID = energyDeposit.trackID;
              auto exist = EDepMap.find(EDepTrackID);
              if (exist == EDepMap.end())
              {
                EDep_TrackID_vec.push_back(EDepTrackID);
                EDepMap[EDepTrackID] = energyDeposit.energy;
              }
              else
              {
                EDepMap[EDepTrackID] += energyDeposit.energy;
              }

            } // end plane == 0
          } // end energy deposit loop
        } // end For each time slice
      } // end For each SimChannel

      fSim_n_hadronic_Edep_b = fSim_hadronic_hit_x_b.size();

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

      const art::FindManyP<simb::MCTruth> findManyTruth(particleHandle, event, fSimulationProducerLabel);

      if (!findManyTruth.isValid())
      {
        std::cout << "findManyTruth simb::MCTruth for simb::MCParticle failed!" << std::endl;
      }

      size_t particle_index = 0;
      auto const &truth = findManyTruth.at(particle_index);

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

  // If this returns true, then the energy deposit is associated with primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle &p1, int primarylep_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    if (MothertrkID == primarylep_trkID)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPrimaryLep(tmp_mother, primarylep_trkID, particleMap);
    }
  }

  bool IsAncestorMotherNeutron(const simb::MCParticle &p1, std::vector<int> neutron_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    for (long unsigned int i = 0; i < neutron_trkID.size(); i++)
    {
      if (MothertrkID == neutron_trkID.at(i))
        MatchMultipleTrkID = true;
    }
    if (MatchMultipleTrkID == true)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherNeutron(tmp_mother, neutron_trkID, particleMap);
    }
  }

  bool IsAncestorMotherProton(const simb::MCParticle &p1, std::vector<int> proton_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    for (long unsigned int i = 0; i < proton_trkID.size(); i++)
    {
      if (MothertrkID == proton_trkID.at(i))
        MatchMultipleTrkID = true;
    }
    if (MatchMultipleTrkID == true)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherProton(tmp_mother, proton_trkID, particleMap);
    }
  }

  bool IsAncestorMotherPip(const simb::MCParticle &p1, std::vector<int> pip_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    for (long unsigned int i = 0; i < pip_trkID.size(); i++)
    {
      if (MothertrkID == pip_trkID.at(i))
        MatchMultipleTrkID = true;
    }
    if (MatchMultipleTrkID == true)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPip(tmp_mother, pip_trkID, particleMap);
    }
  }

  bool IsAncestorMotherPim(const simb::MCParticle &p1, std::vector<int> pim_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    for (long unsigned int i = 0; i < pim_trkID.size(); i++)
    {
      if (MothertrkID == pim_trkID.at(i))
        MatchMultipleTrkID = true;
    }
    if (MatchMultipleTrkID == true)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPim(tmp_mother, pim_trkID, particleMap);
    }
  }

  bool IsAncestorMotherPi0(const simb::MCParticle &p1, std::vector<int> pi0_trkID, std::map<int, const simb::MCParticle *> particleMap)
  {
    int MothertrkID = p1.Mother();
    bool MatchMultipleTrkID = false;
    for (long unsigned int i = 0; i < pi0_trkID.size(); i++)
    {
      if (MothertrkID == pi0_trkID.at(i))
        MatchMultipleTrkID = true;
    }
    if (MatchMultipleTrkID == true)
      return true;
    else if (MothertrkID == 0)
      return false;
    else
    {
      auto tmp_search = particleMap.find(MothertrkID);
      const simb::MCParticle &tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPi0(tmp_mother, pi0_trkID, particleMap);
    }
  }

  // getMassFromPDG: uses V2's corrected strange baryon masses; isotope masses
  // converted from amu to GeV (factor 0.9314941), rounded to 6 sig figs.
  // Nucleon masses from https://www.chemlin.org/chemical-elements/isotopes.php
  double getMassFromPDG(int pdg)
  {
    double mass;
    switch (pdg)
    {
      case 11:   mass = 0.000511;  break; // electron
      case -11:  mass = 0.000511;  break; // positron
      case 12:   mass = 0.0;       break; // electron neutrino
      case -12:  mass = 0.0;       break;
      case 13:   mass = 0.105658;  break; // muon
      case -13:  mass = 0.105658;  break;
      case 14:   mass = 0.0;       break; // muon neutrino
      case -14:  mass = 0.0;       break;
      case 22:   mass = 0.0;       break; // photon
      case 111:  mass = 0.134977;  break; // pi0
      case 211:  mass = 0.139570;  break; // pi+
      case -211: mass = 0.139570;  break; // pi-
      case 221:  mass = 0.547862;  break; // eta
      case -221: mass = 0.547862;  break;
      case 331:  mass = 0.95778;   break; // eta'
      case -331: mass = 0.95778;   break;
      case 321:  mass = 0.493677;  break; // K+
      case -321: mass = 0.493677;  break; // K-
      case 130:  mass = 0.497677;  break; // K_L
      case 310:  mass = 0.497677;  break; // K_S
      case 311:  mass = 0.497677;  break; // K0
      case -311: mass = 0.493677;  break;
      case 2112:  mass = 0.939565; break; // neutron
      case -2112: mass = 0.939565; break; // antineutron
      case 2212:  mass = 0.938272; break; // proton
      case -2212: mass = 0.938272; break; // antiproton
      // Strange baryons: physical masses (from V2)
      case 3122:  mass = 1.11568;  break; // Lambda
      case -3122: mass = 1.11568;  break; // anti-Lambda
      case 3212:  mass = 1.31486;  break; // Sigma0
      case -3212: mass = 0.939565; break; // anti-Sigma0 (V2 value)
      case 3222:  mass = 1.18937;  break; // Sigma+
      case -3222: mass = 1.18937;  break;
      case 3112:  mass = 1.19745;  break; // Sigma-
      case -3112: mass = 1.19745;  break;
      // Isotopes (masses from V2 where they differ, otherwise V1)
      case 1000010020: mass = 1.87561;  break; // deuterium
      case 1000010030: mass = 2.80892;  break; // tritium
      case 1000020030: mass = 2.80839;  break; // He-3
      case 1000020040: mass = 3.72738;  break; // He-4
      case 1000040080: mass = 7.45486;  break; // Be-8
      case 1000040090: mass = 8.39276;  break; // Be-9
      case 1000050100: mass = 9.32444;  break; // B-10
      case 1000050110: mass = 10.2526;  break; // B-11
      case 1000050120: mass = 11.1888;  break; // B-12
      case 1000060110: mass = 10.2540;  break; // C-11
      case 1000060120: mass = 11.1749;  break; // C-12
      case 1000060130: mass = 12.1095;  break; // C-13
      case 1000060140: mass = 13.0409;  break; // C-14
      case 1000070130: mass = 12.1112;  break; // N-13
      case 1000070140: mass = 13.0402;  break; // N-14
      case 1000070150: mass = 13.9690;  break; // N-15
      case 1000070160: mass = 14.9060;  break; // N-16
      case 1000080150: mass = 13.9712;  break; // O-15
      case 1000080160: mass = 14.8951;  break; // O-16
      case 1000080180: mass = 16.7620;  break; // O-18
      case 1000090180: mass = 16.7632;  break; // F-18
      case 1000090190: mass = 17.6923;  break; // F-19
      case 1000100200: mass = 18.6178;  break; // Ne-20
      case 1000100210: mass = 19.5506;  break; // Ne-21
      case 1000100220: mass = 20.4798;  break; // Ne-22
      case 1000110220: mass = 20.4821;  break; // Na-22
      case 1000110230: mass = 21.4092;  break; // Na-23
      case 1000110240: mass = 22.3418;  break; // Na-24
      case 1000110250: mass = 23.2724;  break; // Na-25
      case 1000120220: mass = 20.4864;  break; // Mg-22
      case 1000120230: mass = 21.4128;  break; // Mg-23
      case 1000120240: mass = 22.3358;  break; // Mg-24
      case 1000120250: mass = 23.2680;  break; // Mg-25
      case 1000120260: mass = 24.1965;  break; // Mg-26
      case 1000120270: mass = 25.1297;  break; // Mg-27
      case 1000130260: mass = 24.2000;  break; // Al-26
      case 1000130270: mass = 25.1265;  break; // Al-27
      case 1000130280: mass = 26.0584;  break; // Al-28
      case 1000130290: mass = 26.9885;  break; // Al-29
      case 1000130300: mass = 27.9223;  break; // Al-30
      case 1000130310: mass = 28.8548;  break; // Al-31
      case 1000130320: mass = 29.7901;  break; // Al-32
      case 1000140270: mass = 25.1308;  break; // Si-27
      case 1000140280: mass = 26.0532;  break; // Si-28
      case 1000140290: mass = 26.9843;  break; // Si-29
      case 1000140300: mass = 27.9133;  break; // Si-30
      case 1000140310: mass = 28.8462;  break; // Si-31
      case 1000140320: mass = 29.7766;  break; // Si-32
      case 1000140330: mass = 30.7117;  break; // Si-33
      case 1000150300: mass = 27.9170;  break; // P-30
      case 1000150310: mass = 28.8442;  break; // P-31
      case 1000150320: mass = 29.7759;  break; // P-32
      case 1000150330: mass = 30.7053;  break; // P-33
      case 1000150340: mass = 31.6386;  break; // P-34
      case 1000150350: mass = 32.5698;  break; // P-35
      case 1000150360: mass = 33.5059;  break; // P-36
      case 1000150370: mass = 34.4387;  break; // P-37
      case 1000150380: mass = 35.3745;  break; // P-38
      case 1000160320: mass = 29.7736;  break; // S-32
      case 1000160330: mass = 30.7046;  break; // S-33
      case 1000160340: mass = 31.6327;  break; // S-34
      case 1000160350: mass = 32.5653;  break; // S-35
      case 1000160360: mass = 33.4950;  break; // S-36
      case 1000160370: mass = 34.4302;  break; // S-37
      case 1000160380: mass = 35.3618;  break; // S-38
      case 1000170340: mass = 31.6377;  break; // Cl-34
      case 1000170350: mass = 32.5646;  break; // Cl-35
      case 1000170360: mass = 33.4956;  break; // Cl-36
      case 1000170370: mass = 34.4252;  break; // Cl-37
      case 1000170380: mass = 35.3583;  break; // Cl-38
      case 1000170390: mass = 36.2898;  break; // Cl-39
      case 1000170400: mass = 37.2236;  break; // Cl-40
      case 1000180360: mass = 33.4944;  break; // Ar-36
      case 1000180370: mass = 34.4252;  break; // Ar-37
      case 1000180380: mass = 35.3529;  break; // Ar-38
      case 1000180390: mass = 36.2859;  break; // Ar-39
      case 1000180400: mass = 37.2156;  break; // Ar-40
      case 1000180350: mass = 32.5638;  break; // Ar-35
      case 1000180410: mass = 38.1491;  break; // Ar-41
      case 1000190380: mass = 35.3583;  break; // K-38
      case 1000190390: mass = 36.2848;  break; // K-39
      case 1000190400: mass = 37.2166;  break; // K-40
      case 1000190410: mass = 38.1634;  break; // K-41 (V2 corrected value)
      case 1000200400: mass = 37.2147;  break; // Ca-40
      case 1000220480: mass = 44.6520;  break; // Ti-48
      case 1000220490: mass = 45.5835;  break; // Ti-49
      case 1000230490: mass = 45.5835;  break; // V-49
      case 1000230500: mass = 46.518;   break; // V-50
      case 1000230510: mass = 47.4423;  break; // V-51
      case 1000240490: mass = 45.5857;  break; // Cr-49
      case 1000240500: mass = 46.5122;  break; // Cr-50
      case 1000240510: mass = 47.4425;  break; // Cr-51
      case 1000240520: mass = 48.3701;  break; // Cr-52
      case 1000240530: mass = 49.3017;  break; // Cr-53
      case 1000240540: mass = 50.2315;  break; // Cr-54
      case 1000250530: mass = 49.3018;  break; // Mn-53
      case 1000250540: mass = 50.2324;  break; // Mn-54
      case 1000250550: mass = 51.1617;  break; // Mn-55
      case 1000250560: mass = 52.0940;  break; // Mn-56
      case 1000260520: mass = 48.3761;  break; // Fe-52
      case 1000260530: mass = 49.3050;  break; // Fe-53
      case 1000260540: mass = 50.2312;  break; // Fe-54
      case 1000260550: mass = 51.1615;  break; // Fe-55
      case 1000260560: mass = 52.0898;  break; // Fe-56
      case 1000260570: mass = 53.0217;  break; // Fe-57
      case 1000260580: mass = 53.9513;  break; // Fe-58
      case 1000270570: mass = 53.0221;  break; // Co-57
      case 1000270600: mass = 55.8142;  break; // Co-60
      case 1000280580: mass = 53.9522;  break; // Ni-58
      case 1000280590: mass = 54.8827;  break; // Ni-59
      case 1000280600: mass = 55.8109;  break; // Ni-60
      case 1000280610: mass = 56.7427;  break; // Ni-61
      default:
        std::cerr << "Error in PDG code " << pdg << std::endl;
        return -1.0;
    }
    return mass;
  }

  // --------------------------------------------------------------------------
  // classifyParticle
  // --------------------------------------------------------------------------
  ParticleRole classifyParticle(int pdg)
  {
    // Check for nuclear fragments FIRST, before taking std::abs —
    // nuclear PDG codes are always positive and can be close to INT_MAX,
    // so std::abs on them can overflow a 32-bit int.
    if (pdg >= 1000000000)
      return ParticleRole::NuclearRecoil;

    int absPDG = std::abs(pdg);
    if (absPDG == 12 || absPDG == 14 || absPDG == 16)
      return ParticleRole::Neutrino;
    if (pdg == 111)
      return ParticleRole::Pi0;
    if (pdg == 2112 || pdg == -2112 || pdg == 22 || pdg == 130)
      return ParticleRole::NeutralEscaping;
    if (absPDG == 11 || absPDG == 13 ||
        absPDG == 211 || absPDG == 321 || absPDG == 2212 ||
        absPDG == 3122 || absPDG == 3112 || absPDG == 3222 || absPDG == 3212 ||
        pdg == 310)
      return ParticleRole::VisibleCharged;
    return ParticleRole::Unknown;
  }

  // --------------------------------------------------------------------------
  // getMomentumAtVertex
  //
  // Finds the closest trajectory segment to (vx,vy,vz) and linearly
  // interpolates the 4-momentum along it.
  // Sets 'dies' = true if the vertex coincides with the last trajectory point.
  // --------------------------------------------------------------------------
  TLorentzVector getMomentumAtVertex(const simb::MCParticle* p,
                                      float vx, float vy, float vz,
                                      bool& dies)
  {
    dies = false;
    const unsigned int Ntraj = p->NumberTrajectoryPoints();
    if (Ntraj == 0) return TLorentzVector();
    if (Ntraj == 1) { dies = true; return p->Momentum(0); }

    double minDist = 1e10;
    int    bestIdx = 0;
    for (unsigned int i = 0; i < Ntraj; ++i) {
      auto pos    = p->Position(i);
      double dist = std::hypot(pos.X()-vx, pos.Y()-vy, pos.Z()-vz);
      if (dist < minDist) { minDist = dist; bestIdx = static_cast<int>(i); }
    }

    if (bestIdx == static_cast<int>(Ntraj)-1) {
      dies = true;
      return p->Momentum(bestIdx);
    }

    auto pos0 = p->Position(bestIdx);
    auto pos1 = p->Position(bestIdx+1);
    auto mom0 = p->Momentum(bestIdx);
    auto mom1 = p->Momentum(bestIdx+1);

    double segLen = std::hypot(pos1.X()-pos0.X(), pos1.Y()-pos0.Y(), pos1.Z()-pos0.Z());
    if (segLen < 1e-6) return mom0;

    double dx = vx-pos0.X(), dy = vy-pos0.Y(), dz = vz-pos0.Z();
    double sdx = pos1.X()-pos0.X(), sdy = pos1.Y()-pos0.Y(), sdz = pos1.Z()-pos0.Z();
    double frac = (dx*sdx + dy*sdy + dz*sdz) / (segLen*segLen);
    frac = std::max(0.0, std::min(1.0, frac));

    TLorentzVector interp;
    interp.SetPxPyPzE(
        mom0.Px() + frac*(mom1.Px()-mom0.Px()),
        mom0.Py() + frac*(mom1.Py()-mom0.Py()),
        mom0.Pz() + frac*(mom1.Pz()-mom0.Pz()),
        mom0.E()  + frac*(mom1.E() -mom0.E()));
    return interp;
  }

  // --------------------------------------------------------------------------
  // fillInteractionTree
  //
  // Energy budget per vertex using KE convention throughout:
  //
  //   E_binding = KE_in - sum(KE_out_all) - Q_std
  //
  //   Q_std = sum(M_out_ALL) - M_incoming - M_target
  //
  //   M_target is the implicit target nucleus, reconstructed from
  //   baryon number and charge conservation over ALL outgoing particles
  //   (nuclear fragments + free nucleons).
  //
  // KE convention:
  //   Baryons, mesons, nuclei: contribution = E - M  (KE)
  //   Photons, neutrinos:      contribution = E       (M=0 so KE=E)
  //   Pi0:                     contribution = E - M   (mass goes into Q_std)
  //
  // --------------------------------------------------------------------------
  void fillInteractionTree(
      const simb::MCParticle*                       incoming,
      const Vertex&                                 vertex,
      const std::map<int, const simb::MCParticle*>& particleMap,
      TTree*                                        tree,
      float& fVtx_x,  float& fVtx_y,  float& fVtx_z,  float& fVtx_t,
      int&   fIn_PDG, int&   fIn_TrackID,
      float& fIn_Px,  float& fIn_Py,  float& fIn_Pz,  float& fIn_E,
      float& fIn_Mass,
      std::string& fIn_Process,
      bool&  fIn_Died,
      float& fE_visible,
      float& fE_neutral,
      float& fE_pi0,
      float& fE_recoil,
      float& fE_neutrino,
      float& fE_binding,
      float& fQValue,
      float& fmichDifference,
      std::vector<int>&   fOut_PDG,
      std::vector<float>& fOut_Px,
      std::vector<float>& fOut_Py,
      std::vector<float>& fOut_Pz,
      std::vector<float>& fOut_E,
      std::vector<float>& fOut_Mass,
      std::vector<int>&   fOut_Role,
      std::vector<int>&   fOut_TrackID,
      int& fEvent)
  {
    // ---- Helper lambdas ----
    // Nuclear PDG convention: 10LZZZAAAI  →  Z=(pdg/10000)%1000, A=(pdg/10)%1000
    auto getNuclearZ = [](int pdg) -> int { return (pdg / 10000) % 1000; };
    auto getNuclearA = [](int pdg) -> int { return (pdg / 10)    % 1000; };

    // Baryon number of incoming hadron (0 for mesons/photons/leptons)
    auto getBaryonNumber = [](int pdg) -> int {
      int a = std::abs(pdg);
      if (a == 2212 || a == 2112)                             return 1;
      if (a == 3122 || a == 3112 || a == 3212 || a == 3222)  return 1;
      return 0;
    };

    // Electric charge of incoming hadron
    auto getChargeOf = [](int pdg) -> int {
      switch (pdg) {
        case  2212: return  1;   case -2212: return -1;
        case  2112: return  0;   case -2112: return  0;
        case  3122: return  0;
        case  3222: return  1;   case  3112: return -1;
        case   211: return  1;   case  -211: return -1;
        case   321: return  1;   case  -321: return -1;
        default:    return  0;
      }
    };

    // ---- Incoming particle mass ----
    double inMass = getMassFromPDG(incoming->PdgCode());
    if (inMass < 0.0) {
      inMass = incoming->Mass();
      if (inMass < 0.0) {
        std::cerr << "[fillInteractionTree] Warning: unknown incoming PDG "
                  << incoming->PdgCode() << " — skipping vertex\n";
        return;
      }
    }

    // ---- Interpolate incoming 4-momentum at vertex ----
    bool dies = false;
    TLorentzVector inMom = getMomentumAtVertex(
        incoming, vertex.x, vertex.y, vertex.z, dies);

    double inKE = inMom.E() - inMass;
    if (inKE < -1e-4 && !dies) {
      std::cerr << "[fillInteractionTree] Warning: negative KE ("
                << inKE << " GeV) for PDG " << incoming->PdgCode()
                << " — skipping vertex\n";
      return;
    }

    // ---- Fill incoming particle info ----
    fVtx_x = vertex.x;  fVtx_y = vertex.y;
    fVtx_z = vertex.z;  fVtx_t = vertex.t;
    fIn_PDG     = incoming->PdgCode();
    fIn_TrackID = incoming->TrackId();
    fIn_Px      = static_cast<float>(inMom.Px());
    fIn_Py      = static_cast<float>(inMom.Py());
    fIn_Pz      = static_cast<float>(inMom.Pz());
    fIn_E       = static_cast<float>(inMom.E());
    fIn_Mass    = static_cast<float>(inMass);
    fIn_Process = incoming->EndProcess();
    fIn_Died    = dies;

    // ---- Incoming kinetic energy (KE convention throughout) ----
    // KE = E - M for all massive particles; for photons/neutrinos M=0 so KE=E.
    // Using KE here is consistent with sumOut below — rest masses are
    // accounted for entirely through Q_std.
    const double E_in = inMom.E() - inMass;

    // ---- Clear output vectors ----
    fOut_PDG.clear();
    fOut_Px.clear(); fOut_Py.clear(); fOut_Pz.clear(); fOut_E.clear();
    fOut_Mass.clear(); fOut_Role.clear(); fOut_TrackID.clear();

    // ---- Energy and mass accumulators ----
    double sumVisible    = 0.0;  // KE of charged visible particles
    double sumNeutral    = 0.0;  // KE of neutral escaping (n, gamma)
    double sumPi0        = 0.0;  // KE of pi0 (mass → Q_std)
    double sumRecoil     = 0.0;  // KE of nuclear recoil fragments
    double sumNeutrino   = 0.0;  // total E of neutrinos (M≈0)
    double sumMassAllOut = 0.0;  // sum of ALL outgoing particle masses for Q_std
    int    sumA_all_out  = 0;    // total baryon number of all outgoing particles
    int    sumZ_all_out  = 0;    // total charge of all outgoing particles

    // ---- Loop over daughters ----
    for (const simb::MCParticle* daughter : vertex.daughters)
    {
      if (daughter->TrackId() == incoming->TrackId()) continue;
      if (daughter->Mother()  != incoming->TrackId()) continue;

      double dMass = getMassFromPDG(daughter->PdgCode());
      if (dMass < 0.0) {
        dMass = daughter->Mass();
        if (dMass < 0.0) {
          std::cerr << "[fillInteractionTree] Warning: unknown daughter PDG "
                    << daughter->PdgCode() << " — skipping\n";
          continue;
        }
      }

      const TLorentzVector& dMom = daughter->Momentum(0);
      ParticleRole role = classifyParticle(daughter->PdgCode());

      // KE for all massive particles; E for massless
      double dKE = dMom.E() - dMass;

      switch (role)
      {
        case ParticleRole::VisibleCharged:
          sumVisible += dKE;
          break;
        case ParticleRole::NeutralEscaping:
          // Use KE here — neutron rest mass is NOT deposited and is
          // already accounted for by adding dMass to sumMassAllOut below
          sumNeutral += dKE;
          break;
        case ParticleRole::Pi0:
          // pi0 rest mass goes into Q_std via sumMassAllOut
          sumPi0 += dKE;
          break;
        case ParticleRole::NuclearRecoil:
          sumRecoil += dKE;
          sumA_all_out += getNuclearA(daughter->PdgCode());
          sumZ_all_out += getNuclearZ(daughter->PdgCode());
          break;
        case ParticleRole::Neutrino:
          sumNeutrino += dMom.E();  // M_nu ≈ 0 so E ≈ KE
          break;
        case ParticleRole::Unknown:
          std::cerr << "[fillInteractionTree] Warning: unclassified PDG "
                    << daughter->PdgCode() << " — treating as visible\n";
          sumVisible += std::max(0.0, dKE);
          break;
      }

      // Accumulate ALL outgoing masses for Q_std
      sumMassAllOut += dMass;

      // Track free baryon number and charge for target reconstruction
      // (nuclear fragments already handled above via getNuclearA/Z)
      int dpdg = daughter->PdgCode();
      if      (dpdg ==  2212) { sumA_all_out += 1; sumZ_all_out += 1; }
      else if (dpdg ==  2112) { sumA_all_out += 1; }
      else if (dpdg == -2212) { sumA_all_out -= 1; sumZ_all_out -= 1; }
      else if (dpdg == -2112) { sumA_all_out -= 1; }
      else if (dpdg ==  3122) { sumA_all_out += 1; }                   // Lambda
      else if (dpdg == -3122) { sumA_all_out -= 1; }
      else if (dpdg ==  3222) { sumA_all_out += 1; sumZ_all_out += 1; }// Sigma+
      else if (dpdg ==  3112) { sumA_all_out += 1; sumZ_all_out -= 1; }// Sigma-
      else if (dpdg ==  3212) { sumA_all_out += 1; }                   // Sigma0

      // Fill output vectors
      fOut_PDG.push_back(daughter->PdgCode());
      fOut_Px .push_back(static_cast<float>(dMom.Px()));
      fOut_Py .push_back(static_cast<float>(dMom.Py()));
      fOut_Pz .push_back(static_cast<float>(dMom.Pz()));
      fOut_E  .push_back(static_cast<float>(dMom.E()));
      fOut_Mass.push_back(static_cast<float>(dMass));
      fOut_Role.push_back(static_cast<int>(role));
      fOut_TrackID.push_back(daughter->TrackId());

      // Michel electron detection
      if (std::abs(incoming->PdgCode()) == 13 && std::abs(daughter->PdgCode()) == 11)
        fmichDifference = static_cast<float>(inMom.E() - dMom.E());
    } // end daughter loop

    // ---- Early exit: no valid direct daughters ----
    // Triggered when the vertex is populated only by grandchildren of the
    // incoming particle (mother filter removed everything). Without direct
    // daughters the Q_std formula is undefined and would spuriously give
    // E_binding ≈ M_incoming (≈ 1 GeV for neutrons/protons).
    if (fOut_PDG.empty()) return;

    const bool inIsNucleus = (fIn_PDG >= 1000000000);

    // ---- Early exit: target reconstruction undefined ----
    // A_target < 0 means baryon budget doesn't close — unphysical vertex.
    if (!inIsNucleus) {
      int A_check = sumA_all_out - getBaryonNumber(fIn_PDG);
      if (A_check < 0) return;
    }
    //
    // This is the standard nuclear physics mass-energy balance:
    //   positive Q_std → products are lighter → exothermic → energy released
    //   negative Q_std → products are heavier → endothermic → energy required
    //
    // The target nucleus is reconstructed from baryon/charge conservation
    // over ALL outgoing particles (nuclear fragments + free nucleons).
    //
    double qValue = sumMassAllOut;  // start with sum(M_out_all)

    if (inIsNucleus)
    {
      // Incoming IS a nucleus — its mass is the only incoming mass
      qValue -= inMass;
      // No implicit target: the incoming nucleus IS the target
    }
    else
    {
      // Incoming is a hadron — reconstruct the implicit target nucleus
      // using corrected A/Z that includes all outgoing free nucleons
      int A_target = sumA_all_out - getBaryonNumber(fIn_PDG);
      int Z_target = sumZ_all_out - getChargeOf(fIn_PDG);

      if (A_target > 0 && Z_target >= 0)
      {
        int    targetPDG  = 1000000000 + Z_target * 10000 + A_target * 10;
        double targetMass = getMassFromPDG(targetPDG);

        if (targetMass > 0.0)
        {
          qValue -= targetMass;
        }
        else
        {
          // Isotope not in table — approximate as A * average nucleon mass
          qValue -= A_target * 0.9389;
          std::cerr << "[fillInteractionTree] Note: target PDG " << targetPDG
                    << " (Z=" << Z_target << " A=" << A_target
                    << ") not in mass table, using nucleon approximation\n";
        }
      }
      // Subtract incoming hadron rest mass (part of initial state)
      qValue -= inMass;
    }

    // ---- Binding energy ----
    // E_binding = KE_in - sum(KE_out) - Q_std
    // Should be small and positive (few to ~100 MeV) for physical interactions.
    // Negative values indicate trajectory interpolation error or missing particles.
    double sumOut   = sumVisible + sumNeutral + sumPi0 + sumRecoil + sumNeutrino;
    double eBinding = E_in - sumOut - qValue;

    if (eBinding < -0.2)
      std::cout << "[fillInteractionTree] Note: negative E_binding = "
                << eBinding << " GeV  PDG=" << fIn_PDG
                << "  E_in=" << E_in << "  sumOut=" << sumOut
                << "  Q=" << qValue << "\n";
    if (eBinding > 2.0)
      std::cout << "[fillInteractionTree] Note: large E_binding = "
                << eBinding << " GeV  PDG=" << fIn_PDG << "\n";

    fE_visible  = static_cast<float>(sumVisible);
    fE_neutral  = static_cast<float>(sumNeutral);
    fE_pi0      = static_cast<float>(sumPi0);
    fE_recoil   = static_cast<float>(sumRecoil);
    fE_neutrino = static_cast<float>(sumNeutrino);
    fE_binding  = static_cast<float>(eBinding);
    fQValue     = static_cast<float>(qValue);

    if (!fOut_PDG.empty())
      tree->Fill();
  }
  
  std::vector<Vertex> clusterVertices(const std::vector<const simb::MCParticle*>& daughters)
  {
    std::vector<Vertex> vertices;
    float epsilon = 0.01;

    for (const simb::MCParticle* d : daughters)
    {
      const TLorentzVector& pos = d->Position(0);
      float x = pos.X(), y = pos.Y(), z = pos.Z(), t = pos.T();
      bool found = false;
      for (Vertex& v : vertices)
      {
        if (std::abs(v.x - x) < epsilon && std::abs(v.y - y) < epsilon && std::abs(v.z - z) < epsilon)
        {
          v.daughters.push_back(d);
          found = true;
          break;
        }
      }
      if (!found)
      {
        Vertex vert = {x, y, z, t, {d}};
        vertices.push_back(vert);
      }
    }
    return vertices;
  }

  double getPrimaryKE(const simb::MCParticle* primary, double x, double y, double z)
  {
    double minDist = 1e10;
    int closestDist = 0;

    for (unsigned int n = 0; n < primary->NumberTrajectoryPoints(); ++n)
    {
      const TLorentzVector& position = primary->Position(n);
      double dist = std::sqrt(std::pow(position.X() - x, 2) + std::pow(position.Y() - y, 2) + std::pow(position.Z() - z, 2));
      if (dist < minDist)
      {
        minDist = dist;
        closestDist = n;
      }
    }
    const TLorentzVector& ClosestMom = primary->Momentum(closestDist);
    return ClosestMom.E() - primary->Mass();
  }

  void getHadronic02(const simb::MCParticle* particle, const std::vector<const simb::MCParticle*>& allPart, int& NHad, double& totalBindingE)
  {
    std::vector<const simb::MCParticle*> daughters;
    TLorentzVector currentPos = particle->Position(0);

    for (const simb::MCParticle* p : allPart)
    {
      if (p->Mother() == particle->TrackId())
      {
        daughters.push_back(p);
      }
    }

    if (!daughters.empty())
    {
      std::vector<Vertex> vertices = clusterVertices(daughters);
      float BindingE = 0.0;

      for (const auto& vertex : vertices)
      {
        double Ein  = getPrimaryKE(particle, vertex.x, vertex.y, vertex.z);
        double Eout = 0.0;

        for (const simb::MCParticle* daughter : vertex.daughters)
        {
          if (daughter->PdgCode() == 211)
            Eout += daughter->Momentum(0).E();
          else
            Eout += daughter->Momentum(0).E() - daughter->Mass();
        }
        BindingE = Ein - Eout;
        if (BindingE > 0.001)
        {
          totalBindingE += BindingE;
          NHad++;
        }
      }
    }
    for (const simb::MCParticle* daughter : daughters)
    {
      getHadronic02(daughter, allPart, NHad, totalBindingE);
    }
  }

  void getDescendants(int motherID, const std::vector<int> &momVec, const std::vector<int> &TrkIDvec,
                      const std::map<int, const simb::MCParticle *> &particleMap,
                      std::vector<const simb::MCParticle *> &primaryDaughters)
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

  // Collect all ancestor track IDs of a particle by walking up the mother chain (from V2)
  void getAncestors(const simb::MCParticle *currentpart,
                    std::vector<int> &Mothers,
                    const std::map<int, const simb::MCParticle *> &particleMap)
  {
    if (!currentpart)
      return;

    int momId = currentpart->Mother();
    if (momId <= 0)
      return;
    Mothers.push_back(momId);
    auto it = particleMap.find(momId);
    if (it == particleMap.end())
      return;
    if (it->second == currentpart)
      return;
    getAncestors(it->second, Mothers, particleMap);
  }

  // For root-level (motherless) primary particles, find the first detector-exit
  // trajectory point and return the KE at that point (from V2)
  void ReportFirstExitRootOnly(const simb::MCParticle &part,
                               const std::map<int, const simb::MCParticle *> &particleMap,
                               double &KE)
  {
    const double X_MIN = -400.0, X_MAX = 400.0;
    const double Y_MIN = -600.0, Y_MAX = 600.0;
    const double Z_MIN =    0.0, Z_MAX = 1300.0;

    auto inside = [&](TLorentzVector const &p)
    {
      return (p.X() >= X_MIN && p.X() <= X_MAX) &&
             (p.Y() >= Y_MIN && p.Y() <= Y_MAX) &&
             (p.Z() >= Z_MIN && p.Z() <= Z_MAX);
    };

    // Only act on true root particles (no ancestors in the map)
    std::vector<int> moms;
    getAncestors(&part, moms, particleMap);
    if (!moms.empty())
      return;

    const size_t Ntraj = part.NumberTrajectoryPoints();
    if (Ntraj == 0)
      return;

    bool hasEntered = false;
    for (size_t ipt = 0; ipt < Ntraj; ++ipt)
    {
      const TLorentzVector &pos = part.Position(ipt);
      bool in = inside(pos);
      if (!hasEntered)
      {
        if (in) hasEntered = true;
      }
      else
      {
        if (!in)
        {
          const TLorentzVector &p4 = part.Momentum(ipt);
          KE = p4.E() - part.Mass();
          if (KE < 0) KE = 0;
          return;
        }
      }
    }
  }

} // local namespace
