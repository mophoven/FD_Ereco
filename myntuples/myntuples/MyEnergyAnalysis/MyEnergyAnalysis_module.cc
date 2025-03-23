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
// local namespace

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
      
      if ( !event.getByLabel(fSimulationProducerLabel, particleHandle) ) {
        // If no MCParticles in an event, throw an exception to force this module to stop.
        throw cet::exception("MyEnergyAnalysis") << " No simb::MCParticle objects in this event - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

     std::vector<const simb::MCParticle*> SimParticles;
      // Loop over the list of particles in the event
      // GENIE: primary process; GEANT4: primary+secondary
      for ( auto const& particle : (*particleHandle) ) {

       SimParticles.push_back(&particle);

      } // end loop over all particles in the event
      std::vector<int> daughters;
      for(size_t ip1 = 0, ip1 < SimParticles.size(), ip1++){
        for(size_t ip2 = 0, ip2 < SimParticles.size(), ip2++){
          if(SimParticles(ip2).Mother() == SimParticles(ip1).TrackId() && )
          daughters.push_back(ip2)
        }
      }
     
      
      fNtuple->Fill();

    } // MyEnergyAnalysis::analyze()

    // This macro has to be defined for this module to be invoked from a
    // .fcl file; see MyEnergyAnalysis.fcl for more information.
    DEFINE_ART_MODULE(MyEnergyAnalysis)

  } // namespace example
} // namespace lar

