/**
 * @file   MyEnergyAnalysis_module.cc
 * @brief  Neutrino energy reconstruction / missing-energy accounting for DUNE FD MC.
 *
 * Reworked, "collapsed" version. One pass over the Geant4 MCParticle list builds
 * a genealogy, then three things are produced per event:
 *
 *   (1) EventTree   - one row per event: the closing energy ledger
 *                       E_nu(gen) = E_dep + E_binding + E_escape_neutral
 *                                 + E_escape_charged + E_neutrino_inside + residual
 *                     At truth level (no visibility threshold) the residual is a
 *                     QA number and should close to a few MeV.
 *
 *   (2) VertexTree  - one row per *interesting* interaction vertex (nuclear
 *                     inelastic, capture, elastic, decay). Carries the channel key
 *                     and the per-vertex loss components so the per-interaction
 *                     average-loss table can be built offline by grouping on the key.
 *
 *   (3) EscapeTree  - one row per track that leaves the active volume, charged or
 *                     neutral, with the KE it carried out (KE at the LAST
 *                     inside->outside crossing).
 *
 * Design decisions worth knowing (see inline comments for detail):
 *   - Binding energy is the ground-state nuclear mass change, computed two ways
 *     (mass route, primary; conservation route, cross-check). Disagreement is a flag.
 *   - "Binding" from fillInteractionTree in the old code was identically zero by
 *     construction (it equalled E_in - sumOut - Q_std == 0). That is removed.
 *   - Deposited energy comes from sim::SimEnergyDeposit (true ionization in the
 *     active volume), NOT from a SimChannel/collection-plane loop.
 *   - Escape is a per-track geometric quantity (last exit KE), so neutron energy is
 *     never summed across the cascade (that double-counted in the old code).
 *   - Photons are EM/visible; only their geometric exit counts as escape. Neutrinos
 *     are handled by the neutrino term, and only if created INSIDE the active volume
 *     (a muon that escapes and decays outside already had its energy counted as
 *     charged escape; its decay neutrinos must not be double counted).
 *
 * @author (rework) — building on W. Shi's original; original adapted from
 *         the LArSoft AnalysisExample.
 */

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// art
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TLorentzVector.h"
#include "TTree.h"

// C++
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <functional>

namespace
{

  // ------------------------------------------------------------------
  //  Process categorization
  // ------------------------------------------------------------------
  // The Geant4 *creation process* of a daughter is the channel that made it.
  // We key interactions on that string (n -> daughters created by
  // "neutronInelastic", etc.) rather than guessing from product kinematics.

  // Processes we treat as genuine interaction vertices worth recording.
  const std::set<std::string> &interestingProcesses()
  {
    static const std::set<std::string> s = {
        "neutronInelastic", "protonInelastic",
        "pi+Inelastic", "pi-Inelastic", "pi0Inelastic",
        "kaon+Inelastic", "kaon-Inelastic", "kaon0LInelastic", "kaon0SInelastic",
        "dInelastic", "tInelastic", "He3Inelastic", "alphaInelastic",
        "hadElastic",
        "nCapture", "muMinusCaptureAtRest", "hBertiniCaptureAtRest",
        "Decay", "muMinusDecay",
        "photonNuclear", "muonNuclear"};
    return s;
  }

  bool isDecayProcess(const std::string &p)
  {
    return p == "Decay" || p == "muMinusDecay";
  }

  // Stopping particle absorbed by a nucleus (mu-, stopped hadrons). The PROJECTILE
  // is consumed and (for mu-) most of its rest mass leaves as a neutrino. This is
  // NOT a binding-energy situation — treat like a decay/conversion.
  bool isCaptureAtRest(const std::string &p)
  {
    return p.find("CaptureAtRest") != std::string::npos;
  }

  // A genuine nuclear-target interaction where binding/Q makes sense
  // (in-flight inelastic on a nucleus, or neutron radiative capture nCapture).
  // Excludes decays, capture-at-rest, and elastic (handled separately).
  bool isNuclearTargetProcess(const std::string &p)
  {
    if (isDecayProcess(p))
      return false;
    if (isCaptureAtRest(p))
      return false;
    if (p == "hadElastic")
      return false;
    return interestingProcesses().count(p) > 0;
  }

  // Rest-mass convention for escape: a particle CREATED in the interaction
  // (lepton, meson, photon, hyperon) carries rest mass that came from the
  // neutrino energy, so its full energy E is lost when it exits. A pre-existing
  // nucleon (p, n) or nucleus carries rest mass that was already in the target
  // nucleus (accounted via binding), so only its KE is lost on exit.
  bool isCreatedParticle(int pdg)
  {
    if (pdg == 2212 || pdg == 2112)
      return false; // free nucleons
    if (pdg >= 1000000000)
      return false; // nuclei / fragments (isNucleus)
    return true;    // leptons, mesons, photons, hyperons
  }

  // ------------------------------------------------------------------
  //  Nuclear / hadron bookkeeping helpers
  // ------------------------------------------------------------------
  inline bool isNucleus(int pdg) { return pdg >= 1000000000; }
  inline int nuclearZ(int pdg) { return (pdg / 10000) % 1000; }
  inline int nuclearA(int pdg) { return (pdg / 10) % 1000; }
  inline int makeNucleusPDG(int Z, int A) { return 1000000000 + Z * 10000 + A * 10; }
  inline bool isNeutrino(int pdg)
  {
    int a = std::abs(pdg);
    return a == 12 || a == 14 || a == 16;
  }

  // Baryon number for the incoming hadron (nuclei handled separately).
  int baryonNumber(int pdg)
  {
    int a = std::abs(pdg);
    int s = (pdg < 0) ? -1 : 1;
    if (a == 2212 || a == 2112)
      return s; // p, n
    if (a == 3122 || a == 3112 || a == 3212 || a == 3222)
      return s; // Lambda, Sigma
    return 0;
  }

  // Electric charge of a (non-nuclear) particle.
  int chargeOf(int pdg)
  {
    switch (pdg)
    {
    case 2212:
      return 1;
    case -2212:
      return -1;
    case 2112:
      return 0;
    case -2112:
      return 0;
    case 211:
      return 1;
    case -211:
      return -1;
    case 321:
      return 1;
    case -321:
      return -1;
    case 3222:
      return 1;
    case 3112:
      return -1;
    case 3122:
      return 0;
    case 3212:
      return 0;
    case 13:
      return -1;
    case -13:
      return 1;
    case 11:
      return -1;
    case -11:
      return 1;
    default:
      return 0;
    }
  }

  // ------------------------------------------------------------------
  //  Mass table (GeV). Nuclear masses are full atomic-mass-unit conversions.
  //  Used for the *mass-route* binding energy. Falls back to MCParticle::Mass()
  //  for anything not listed.  (-3212 anti-Sigma0 fixed vs. the old copy/paste.)
  // ------------------------------------------------------------------
  double getMassFromPDG(int pdg)
  {
    switch (pdg)
    {
    case 11:
    case -11:
      return 0.000511;
    case 12:
    case -12:
      return 0.0;
    case 13:
    case -13:
      return 0.105658;
    case 14:
    case -14:
      return 0.0;
    case 16:
    case -16:
      return 0.0;
    case 22:
      return 0.0;
    case 111:
      return 0.134977;
    case 211:
    case -211:
      return 0.139570;
    case 221:
      return 0.547862;
    case 321:
    case -321:
      return 0.493677;
    case 130:
    case 310:
    case 311:
    case -311:
      return 0.497611;
    case 2112:
    case -2112:
      return 0.939565;
    case 2212:
    case -2212:
      return 0.938272;
    case 3122:
    case -3122:
      return 1.115683;
    case 3212:
      return 1.192642;
    case -3212:
      return 1.192642; // anti-Sigma0 (was neutron mass by mistake)
    case 3222:
    case -3222:
      return 1.189370;
    case 3112:
    case -3112:
      return 1.197449;
    // Nuclear masses (GeV, ground-state, electrons removed), derived uniformly
    // from AME2020 mass excesses via  M = [ A*u + Delta - Z*m_e ],
    //   u = 931.4941024 MeV,  m_e = 0.5109989 MeV.  Sorted by Z then A.
    // (Inverse for checking against the NNDC chart: Delta = 1000*M - A*u + Z*m_e.)
    case 1000010020:
      return 1.87561; // d
    case 1000010030:
      return 2.80892; // t
    case 1000020030:
      return 2.80839; // He-3
    case 1000020040:
      return 3.72738; // He-4
    case 1000030070:
      return 6.53383; // Li-7
    case 1000040070:
      return 6.53418; // Be-7
    case 1000040080:
      return 7.45485; // Be-8
    case 1000040090:
      return 8.39275; // Be-9
    case 1000040100:
      return 9.32550; // Be-10
    case 1000050100:
      return 9.32444; // B-10
    case 1000050110:
      return 10.25255; // B-11
    case 1000060110:
      return 10.25402; // C-11
    case 1000060120:
      return 11.17486; // C-12
    case 1000060130:
      return 12.10948; // C-13
    case 1000060140:
      return 13.04087; // C-14
    case 1000070130:
      return 12.11119; // N-13
    case 1000070140:
      return 13.04020; // N-14
    case 1000070150:
      return 13.96894; // N-15
    case 1000070160:
      return 14.90601; // N-16
    case 1000080140:
      return 13.04484; // O-14
    case 1000080150:
      return 13.97118; // O-15
    case 1000080160:
      return 14.89508; // O-16
    case 1000080170:
      return 15.83050; // O-17
    case 1000080180:
      return 16.76202; // O-18
    case 1000090190:
      return 17.69230; // F-19
    case 1000100200:
      return 18.61773; // Ne-20
    case 1000100210:
      return 19.55053; // Ne-21
    case 1000100220:
      return 20.47974; // Ne-22
    case 1000110220:
      return 20.48207; // Na-22
    case 1000110230:
      return 21.40921; // Na-23
    case 1000110240:
      return 22.34182; // Na-24
    case 1000110250:
      return 23.27237; // Na-25
    case 1000120240:
      return 22.33579; // Mg-24
    case 1000120250:
      return 23.26803; // Mg-25
    case 1000120260:
      return 24.19650; // Mg-26
    case 1000120270:
      return 25.12962; // Mg-27
    case 1000120280:
      return 26.06068; // Mg-28
    case 1000130250:
      return 23.27179; // Al-25
    case 1000130260:
      return 24.19999; // Al-26
    case 1000130270:
      return 25.12650; // Al-27
    case 1000130280:
      return 26.05834; // Al-28
    case 1000130300:
      return 27.92231; // Al-30
    case 1000140270:
      return 25.13080; // Si-27
    case 1000140280:
      return 26.05319; // Si-28
    case 1000140290:
      return 26.98428; // Si-29
    case 1000140300:
      return 27.91324; // Si-30
    case 1000140310:
      return 28.84621; // Si-31
    case 1000140320:
      return 29.77658; // Si-32
    case 1000150300:
      return 27.91696; // P-30
    case 1000150310:
      return 28.84421; // P-31
    case 1000150320:
      return 29.77584; // P-32
    case 1000150330:
      return 30.70530; // P-33
    case 1000150340:
      return 31.63858; // P-34
    case 1000160320:
      return 29.77362; // S-32
    case 1000160330:
      return 30.70454; // S-33
    case 1000160340:
      return 31.63269; // S-34
    case 1000160350:
      return 32.56527; // S-35
    case 1000160360:
      return 33.49495; // S-36
    case 1000160370:
      return 34.43021; // S-37
    case 1000160380:
      return 35.36174; // S-38
    case 1000170350:
      return 32.56459; // Cl-35
    case 1000170360:
      return 33.49558; // Cl-36
    case 1000170370:
      return 34.42483; // Cl-37
    case 1000170380:
      return 35.35829; // Cl-38
    case 1000170390:
      return 36.28978; // Cl-39
    case 1000170400:
      return 37.22352; // Cl-40
    case 1000180350:
      return 32.57005; // Ar-35
    case 1000180360:
      return 33.49436; // Ar-36
    case 1000180370:
      return 34.42514; // Ar-37
    case 1000180380:
      return 35.35286; // Ar-38
    case 1000180390:
      return 36.28583; // Ar-39
    case 1000180400:
      return 37.21553; // Ar-40
    case 1000180410:
      return 38.14899; // Ar-41
    case 1000190390:
      return 36.28475; // K-39
    case 1000190400:
      return 37.21652; // K-40
    case 1000190410:
      return 38.14599; // K-41
    case 1000200400:
      return 37.21470; // Ca-40
    case 1000200420:
      return 39.07399; // Ca-42
    case 1000200440:
      return 40.93405; // Ca-44
    case 1000210450:
      return 41.86544; // Sc-45
    case 1000220440:
      return 40.93695; // Ti-44
    case 1000220480:
      return 44.65198; // Ti-48
    case 1000230480:
      return 44.65549; // V-48
    case 1000230500:
      return 46.51373; // V-50
    case 1000230510:
      return 47.44224; // V-51
    case 1000240490:
      return 45.58562; // Cr-49
    case 1000240500:
      return 46.51218; // Cr-50
    case 1000240520:
      return 48.37001; // Cr-52
    case 1000240530:
      return 49.30164; // Cr-53
    case 1000240540:
      return 50.23149; // Cr-54
    case 1000250530:
      return 49.30172; // Mn-53
    case 1000250540:
      return 50.23235; // Mn-54
    case 1000250550:
      return 51.16169; // Mn-55
    case 1000250560:
      return 52.09398; // Mn-56
    case 1000260530:
      return 49.30496; // Fe-53
    case 1000260540:
      return 50.23114; // Fe-54
    case 1000260550:
      return 51.16141; // Fe-55
    case 1000260560:
      return 52.08978; // Fe-56
    case 1000260570:
      return 53.02170; // Fe-57
    case 1000260580:
      return 53.95122; // Fe-58
    case 1000280580:
      return 53.95212; // Ni-58
    case 1000280590:
      return 54.88269; // Ni-59
    default:
      return -1.0; // not in table -> mass route abandoned
    }
  }

  std::map<int, long> &missingMassRegistry()
  {
    static std::map<int, long> m;
    return m;
  }

  // Robust particle mass (GeV). Source priority:
  //   1. getMassFromPDG table (exact mass excesses — needed for accurate binding)
  //   2. the MCParticle's own mass (Geant4's value), if it looks sane
  //   3. for nuclei only, A * u as a crude bound so a miss can never inject a
  //      multi-GeV rest mass into the ledger (binding for that vertex is then only
  //      approximate, but bounded — the guard on the binding sum catches the rest).
  // mcMass is the depositing particle's MCParticle::Mass(); pass <0 if unavailable
  // (e.g. for a reconstructed target nucleus, which then relies on the table only).
  double robustMass(int pdg, double mcMass)
  {
    double m = getMassFromPDG(pdg);
    if (m >= 0.0)
      return m;
    ++missingMassRegistry()[pdg];
    if (pdg >= 1000000000)
    {                                                    // nucleus
      const double expected = nuclearA(pdg) * 0.9314941; // u in GeV ~ A nucleons
      if (mcMass > 0.5 * expected && mcMass < 1.5 * expected)
        return mcMass;
      return expected; // crude, bounded
    }
    if (mcMass > 0.0)
      return mcMass; // trust Geant4 for non-nuclei
    return 0.0;      // unknown -> treat as massless
  }

  // KE of an MCParticle at trajectory point i (GeV).
  inline double pointKE(const simb::MCParticle &p, unsigned int i)
  {
    return p.Momentum(i).E() - p.Mass();
  }

  double incomingKEAtVertex(const simb::MCParticle &p,
                            double vx, double vy, double vz)
  {
    const unsigned int N = p.NumberTrajectoryPoints();
    if (N == 0)
      return 0.0;

    unsigned int iv = 0;
    double best = 1e30;
    for (unsigned int i = 0; i < N; ++i)
    {
      const TLorentzVector &q = p.Position(i);
      double d = std::hypot(q.X() - vx, q.Y() - vy, q.Z() - vz);
      if (d < best)
      {
        best = d;
        iv = i;
      }
    }
    double ke = pointKE(p, iv);
    if (iv > 0)
      ke = std::max(ke, pointKE(p, iv - 1));
    return ke;
  }

  // A clustered interaction vertex: a position where one parent created daughters
  // via one process.
  struct Vertex
  {
    double x, y, z, t;
    std::string process; // daughter creation process
    std::vector<const simb::MCParticle *> daughters;
  };

} // anonymous namespace

namespace lar
{
  namespace example
  {

    class MyEnergyAnalysis : public art::EDAnalyzer
    {
    public:
      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> GenieGenModuleLabel{
            Name("GenieGenModuleLabel"), Comment("generator (GENIE) label"),
            art::InputTag("generator")};

        fhicl::Atom<art::InputTag> SimulationLabel{
            Name("SimulationLabel"), Comment("Geant4 MCParticle label"),
            art::InputTag("largeant")};

        fhicl::Atom<art::InputTag> SimEnergyDepositLabel{
            Name("SimEnergyDepositLabel"),
            Comment("sim::SimEnergyDeposit label — VERIFY with eventdump"),
            art::InputTag("largeant:TPCActive")};

        // Active-volume bounds (cm). Defaults are the user-stated DUNE values;
        // beginJob prints geometry dimensions so you can reconcile them.
        fhicl::Atom<double> ActiveXmin{Name("ActiveXmin"), Comment("cm"), -600.0};
        fhicl::Atom<double> ActiveXmax{Name("ActiveXmax"), Comment("cm"), 600.0};
        fhicl::Atom<double> ActiveYmin{Name("ActiveYmin"), Comment("cm"), -375.0};
        fhicl::Atom<double> ActiveYmax{Name("ActiveYmax"), Comment("cm"), 375.0};
        fhicl::Atom<double> ActiveZmin{Name("ActiveZmin"), Comment("cm"), 0.0};
        fhicl::Atom<double> ActiveZmax{Name("ActiveZmax"), Comment("cm"), 1400.0};

        fhicl::Atom<bool> SelectCC{
            Name("SelectCC"),
            Comment("if true, analyze only charged-current events (skip NC entirely)"),
            true};

        fhicl::Atom<double> FiducialInset{
            Name("FiducialInset"),
            Comment("cm to shrink the active volume by for the vertex fiducial cut (0 = use active bounds)"),
            0.0};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MyEnergyAnalysis(Parameters const &config);
      void beginJob() override;
      void endJob() override;
      void analyze(const art::Event &event) override;

    private:
      bool inside(double x, double y, double z) const
      {
        return x > fXmin && x < fXmax && y > fYmin && y < fYmax && z > fZmin && z < fZmax;
      }
      bool inside(const TLorentzVector &p) const { return inside(p.X(), p.Y(), p.Z()); }

      // Walk to the GENIE-primary ancestor and return its PDG (memoized).
      int rootPrimaryPDG(int trackID);

      // Category index for the deposited-energy breakdown (by primary ancestor).
      int depositCategory(int primaryPDG) const;

      // Labels & geometry
      art::InputTag fGenLabel, fSimLabel, fEdepLabel;
      double fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;
      bool fSelectCC = true;
      long fNskippedNC = 0;
      double fFidInset = 0.0;
      long fNskippedFid = 0;
      geo::GeometryCore const *fGeom = nullptr;

      // Job-level extent of all SimEnergyDeposits (empirical active-volume probe).
      double fDepXmin = 1e30, fDepXmax = -1e30;
      double fDepYmin = 1e30, fDepYmax = -1e30;
      double fDepZmin = 1e30, fDepZmax = -1e30;
      long fDepCount = 0;

      // Per-event scratch shared across helpers
      std::map<int, const simb::MCParticle *> fPmap; // trackID -> particle
      std::map<int, int> fRootCache;

      // ---- Trees ----
      TTree *fEventTree = nullptr;
      TTree *fVertexTree = nullptr;
      TTree *fEscapeTree = nullptr;

      // EventTree branches
      int fEvent, fRun, fSubRun;
      double fGen_nu_E;
      int fGen_nu_PDG, fCCNC, fMode, fInteractionType;
      double fNuVtxX, fNuVtxY, fNuVtxZ;
      int fLep_PDG;
      double fLep_E;

      double fE_dep_total;
      double fE_dep_mu, fE_dep_p, fE_dep_n, fE_dep_pi, fE_dep_em, fE_dep_nuc, fE_dep_other;
      // True category kinetic-energy denominators.
      double fE_true_mu, fE_true_p, fE_true_n, fE_true_pi, fE_true_em, fE_true_nuc, fE_true_other;

      double fE_binding_total;
      double fE_escape_neutral, fE_escape_charged;
      double fE_neutrino_inside;
      double fE_residual;
      double fEscapeFrac; // (E_escape_charged + E_escape_neutral) / E_nu — containment flag
      int fN_michel_inside;
      double fMichel_nu_E;

      // VertexTree branches
      int fV_event;
      double fV_x, fV_y, fV_z, fV_t;
      int fV_in_pdg, fV_in_trk;
      double fV_in_KE;
      std::string fV_process;
      bool fV_inside;
      int fV_target_pdg, fV_residual_pdg;
      int fV_nOut, fV_nNeutron, fV_nProton, fV_nGamma, fV_nNeutrino, fV_nNucleus;
      double fV_Ebind_mass, fV_Ebind_cons, fV_Ebind_diff;
      double fV_Ebind_meson, fV_Ebind_nuclear;
      double fV_Enu_vtx;
      bool fV_isMichel;
      std::string fV_channel;
      unsigned long long fV_channel_hash;
      std::vector<int> fV_out_pdg;
      std::vector<double> fV_out_KE;

      // EscapeTree branches
      int fX_event, fX_trk, fX_pdg;
      double fX_birthKE, fX_exitKE, fX_exitX, fX_exitY, fX_exitZ;
      std::string fX_endProcess;
      std::string fX_birthChannel; // channel of the vertex that created this track
      bool fX_charged;
    };

    // ----------------------------------------------------------------------------
    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const &c)
        : EDAnalyzer(c), fGenLabel(c().GenieGenModuleLabel()), fSimLabel(c().SimulationLabel()), fEdepLabel(c().SimEnergyDepositLabel()), fXmin(c().ActiveXmin()), fXmax(c().ActiveXmax()), fYmin(c().ActiveYmin()), fYmax(c().ActiveYmax()), fZmin(c().ActiveZmin()), fZmax(c().ActiveZmax()), fSelectCC(c().SelectCC()), fFidInset(c().FiducialInset())
    {
      fGeom = &*art::ServiceHandle<geo::Geometry>();
      consumes<std::vector<simb::MCTruth>>(fGenLabel);
      consumes<std::vector<simb::MCParticle>>(fSimLabel);
      consumes<std::vector<sim::SimEnergyDeposit>>(fEdepLabel);
    }

    // ----------------------------------------------------------------------------
    void MyEnergyAnalysis::beginJob()
    {
      mf::LogInfo("MyEnergyAnalysis")
          << "Active-volume bounds in use (cm): "
          << "x[" << fXmin << "," << fXmax << "] "
          << "y[" << fYmin << "," << fYmax << "] "
          << "z[" << fZmin << "," << fZmax << "]\n"
          << "Geometry reference: DetLength=" << fGeom->DetLength()
          << " 2*DetHalfWidth=" << 2 * fGeom->DetHalfWidth()
          << " 2*DetHalfHeight=" << 2 * fGeom->DetHalfHeight() << " cm\n"
          << "SimEnergyDeposit label: " << fEdepLabel.encode()
          << "  (verify this exists in your file)";

      art::ServiceHandle<art::TFileService const> tfs;

      // ---- EventTree ----
      fEventTree = tfs->make<TTree>("EventTree", "per-event energy ledger");
      fEventTree->Branch("Event", &fEvent);
      fEventTree->Branch("Run", &fRun);
      fEventTree->Branch("SubRun", &fSubRun);
      fEventTree->Branch("Gen_nu_E", &fGen_nu_E);
      fEventTree->Branch("Gen_nu_PDG", &fGen_nu_PDG);
      fEventTree->Branch("CCNC", &fCCNC);
      fEventTree->Branch("Mode", &fMode);
      fEventTree->Branch("InteractionType", &fInteractionType);
      fEventTree->Branch("NuVtxX", &fNuVtxX);
      fEventTree->Branch("NuVtxY", &fNuVtxY);
      fEventTree->Branch("NuVtxZ", &fNuVtxZ);
      fEventTree->Branch("Lep_PDG", &fLep_PDG);
      fEventTree->Branch("Lep_E", &fLep_E);
      fEventTree->Branch("E_dep_total", &fE_dep_total);
      fEventTree->Branch("E_dep_mu", &fE_dep_mu);
      fEventTree->Branch("E_dep_p", &fE_dep_p);
      fEventTree->Branch("E_dep_n", &fE_dep_n);
      fEventTree->Branch("E_dep_pi", &fE_dep_pi);
      fEventTree->Branch("E_dep_em", &fE_dep_em);
      fEventTree->Branch("E_dep_nuc", &fE_dep_nuc);
      fEventTree->Branch("E_dep_other", &fE_dep_other);

      fEventTree->Branch("E_true_mu", &fE_true_mu);
      fEventTree->Branch("E_true_p", &fE_true_p);
      fEventTree->Branch("E_true_n", &fE_true_n);
      fEventTree->Branch("E_true_pi", &fE_true_pi);
      fEventTree->Branch("E_true_em", &fE_true_em);
      fEventTree->Branch("E_true_nuc", &fE_true_nuc);
      fEventTree->Branch("E_true_other", &fE_true_other);

      fEventTree->Branch("E_binding_total", &fE_binding_total);
      fEventTree->Branch("E_escape_neutral", &fE_escape_neutral);
      fEventTree->Branch("E_escape_charged", &fE_escape_charged);
      fEventTree->Branch("E_neutrino_inside", &fE_neutrino_inside);
      fEventTree->Branch("E_residual", &fE_residual);
      fEventTree->Branch("EscapeFrac", &fEscapeFrac);
      fEventTree->Branch("N_michel_inside", &fN_michel_inside);
      fEventTree->Branch("Michel_nu_E", &fMichel_nu_E);

      // ---- VertexTree ----
      fVertexTree = tfs->make<TTree>("VertexTree", "per-interaction-vertex channels");
      fVertexTree->Branch("Event", &fV_event);
      fVertexTree->Branch("Vtx_x", &fV_x);
      fVertexTree->Branch("Vtx_y", &fV_y);
      fVertexTree->Branch("Vtx_z", &fV_z);
      fVertexTree->Branch("Vtx_t", &fV_t);
      fVertexTree->Branch("In_PDG", &fV_in_pdg);
      fVertexTree->Branch("In_TrackID", &fV_in_trk);
      fVertexTree->Branch("In_KE", &fV_in_KE);
      fVertexTree->Branch("Process", &fV_process);
      fVertexTree->Branch("Inside", &fV_inside);
      fVertexTree->Branch("Target_PDG", &fV_target_pdg);
      fVertexTree->Branch("Residual_PDG", &fV_residual_pdg);
      fVertexTree->Branch("nOut", &fV_nOut);
      fVertexTree->Branch("nNeutron", &fV_nNeutron);
      fVertexTree->Branch("nProton", &fV_nProton);
      fVertexTree->Branch("nGamma", &fV_nGamma);
      fVertexTree->Branch("nNeutrino", &fV_nNeutrino);
      fVertexTree->Branch("nNucleus", &fV_nNucleus);
      fVertexTree->Branch("E_binding_mass", &fV_Ebind_mass);
      fVertexTree->Branch("E_binding_cons", &fV_Ebind_cons);
      fVertexTree->Branch("E_binding_diff", &fV_Ebind_diff);
      fVertexTree->Branch("E_binding_meson", &fV_Ebind_meson);     // created-meson rest mass removed
      fVertexTree->Branch("E_binding_nuclear", &fV_Ebind_nuclear); // nuclear remainder (added to ledger)
      fVertexTree->Branch("E_neutrino_vtx", &fV_Enu_vtx);
      fVertexTree->Branch("IsMichel", &fV_isMichel);
      fVertexTree->Branch("Channel", &fV_channel);
      fVertexTree->Branch("ChannelHash", &fV_channel_hash);
      fVertexTree->Branch("Out_PDG", &fV_out_pdg);
      fVertexTree->Branch("Out_KE", &fV_out_KE);

      // ---- EscapeTree ----
      fEscapeTree = tfs->make<TTree>("EscapeTree", "tracks leaving the active volume");
      fEscapeTree->Branch("Event", &fX_event);
      fEscapeTree->Branch("TrackID", &fX_trk);
      fEscapeTree->Branch("PDG", &fX_pdg);
      fEscapeTree->Branch("BirthKE", &fX_birthKE);
      fEscapeTree->Branch("ExitKE", &fX_exitKE);
      fEscapeTree->Branch("ExitX", &fX_exitX);
      fEscapeTree->Branch("ExitY", &fX_exitY);
      fEscapeTree->Branch("ExitZ", &fX_exitZ);
      fEscapeTree->Branch("EndProcess", &fX_endProcess);
      fEscapeTree->Branch("BirthChannel", &fX_birthChannel);
      fEscapeTree->Branch("Charged", &fX_charged);
    }

    // ----------------------------------------------------------------------------
    void MyEnergyAnalysis::endJob()
    {

      if (fDepCount > 0)
      {
        mf::LogInfo("MyEnergyAnalysis")
            << "SimEnergyDeposit spatial extent over the job (" << fDepCount << " deposits):\n"
            << "    x [" << fDepXmin << ", " << fDepXmax << "] cm\n"
            << "    y [" << fDepYmin << ", " << fDepYmax << "] cm\n"
            << "    z [" << fDepZmin << ", " << fDepZmax << "] cm\n"
            << "  Escape bounds configured: "
            << "x[" << fXmin << "," << fXmax << "] "
            << "y[" << fYmin << "," << fYmax << "] "
            << "z[" << fZmin << "," << fZmax << "]\n"
            << "  -> set Active{X,Y,Z}{min,max} to match the deposit extent above.";
      }
      else
      {
        mf::LogWarning("MyEnergyAnalysis")
            << "No SimEnergyDeposits seen all job — check SimEnergyDepositLabel.";
      }

      if (fSelectCC)
        mf::LogInfo("MyEnergyAnalysis")
            << "CC-only selection: skipped " << fNskippedNC << " non-CC events.";
      mf::LogInfo("MyEnergyAnalysis")
          << "Fiducial cut (inset " << fFidInset << " cm): skipped " << fNskippedFid
          << " events with vertex outside the active volume.";

      auto const &miss = missingMassRegistry();
      if (!miss.empty())
      {
        std::ostringstream os;
        os << "PDG codes not in getMassFromPDG (mass came from MCParticle/approx):";
        for (auto const &kv : miss)
        {
          int pdg = kv.first;
          os << "\n    " << pdg << "  (x" << kv.second << ")";
          if (pdg >= 1000000000)
            os << "  nucleus Z=" << ((pdg / 10000) % 1000) << " A=" << ((pdg / 10) % 1000);
        }
        mf::LogWarning("MyEnergyAnalysis") << os.str();
      }
    }

    // ----------------------------------------------------------------------------
    int MyEnergyAnalysis::rootPrimaryPDG(int trackID)
    {
      auto cached = fRootCache.find(trackID);
      if (cached != fRootCache.end())
        return cached->second;

      int cur = trackID, result = 0;
      std::vector<int> visited;
      while (true)
      {
        auto it = fPmap.find(cur);
        if (it == fPmap.end())
        {
          result = 0;
          break;
        }
        visited.push_back(cur);
        int mom = it->second->Mother();
        if (mom == 0 || fPmap.find(mom) == fPmap.end())
        {
          result = it->second->PdgCode();
          break;
        }
        cur = mom;
        if (visited.size() > 10000)
        {
          result = it->second->PdgCode();
          break;
        } // loop guard
      }
      for (int v : visited)
        fRootCache[v] = result;
      return result;
    }

    int MyEnergyAnalysis::depositCategory(int pdg) const
    {
      int a = std::abs(pdg);
      if (a == 13)
        return 0; // muon
      if (pdg == 2212)
        return 1; // proton
      if (pdg == 2112)
        return 2; // neutron
      if (a == 211)
        return 3; // charged pion
      if (a == 11 || pdg == 22 || pdg == 111)
        return 4; // EM
      if (isNucleus(pdg))
        return 5; // nuclear
      return 6;   // other
    }

    // ----------------------------------------------------------------------------
    void MyEnergyAnalysis::analyze(const art::Event &event)
    {
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      // reset ledger
      fGen_nu_E = 0;
      fGen_nu_PDG = 0;
      fCCNC = -1;
      fMode = -1;
      fInteractionType = -1;
      fNuVtxX = fNuVtxY = fNuVtxZ = -9999;
      fLep_PDG = 0;
      fLep_E = 0;
      fE_dep_total = 0;
      fE_dep_mu = fE_dep_p = fE_dep_n = fE_dep_pi = 0;
      fE_dep_em = fE_dep_nuc = fE_dep_other = 0;
      fE_true_mu = fE_true_p = fE_true_n = fE_true_pi = 0;
      fE_true_em = fE_true_nuc = fE_true_other = 0;
      fE_binding_total = 0;
      fE_escape_neutral = 0;
      fE_escape_charged = 0;
      fE_neutrino_inside = 0;
      fE_residual = 0;
      fN_michel_inside = 0;
      fMichel_nu_E = 0;
      fPmap.clear();
      fRootCache.clear();

      // trackID -> channel string of the vertex that created it. Populated in stage 4,
      // read in stage 5 so each escaping track can be attributed to its production
      // interaction. Empty for tracks not born at an interesting vertex.
      std::map<int, std::string> birthChannel;

      // ---- (1) generator truth ----
      art::Handle<std::vector<simb::MCTruth>> mcth;
      if (event.getByLabel(fGenLabel, mcth) && !mcth->empty())
      {
        // Dereference the handle directly — no art::Ptr / fill_ptr_vector needed.
        const auto &nu = mcth->at(0).GetNeutrino();
        fGen_nu_E = nu.Nu().E();
        fGen_nu_PDG = nu.Nu().PdgCode();
        fCCNC = nu.CCNC();
        fMode = nu.Mode();
        fInteractionType = nu.InteractionType();
        fNuVtxX = nu.Nu().Vx();
        fNuVtxY = nu.Nu().Vy();
        fNuVtxZ = nu.Nu().Vz();
        fLep_PDG = nu.Lepton().PdgCode();
        fLep_E = nu.Lepton().E();
        // In an NC interaction the "outgoing lepton" is a neutrino that Geant4 does NOT
        // propagate, so it never appears in the MCParticle list scanned in stage 6.
        // Add it here from the generator truth (it is created at the interaction vertex,
        // which is inside the detector). CC leptons are handled via deposits/escape.
        if (fCCNC == 1 && isNeutrino(fLep_PDG) && inside(fNuVtxX, fNuVtxY, fNuVtxZ))
          fE_neutrino_inside += fLep_E;
      }

      // CC-only selection: skip NC (and truth-less) events entirely — no trees filled.
      if (fSelectCC && fCCNC != 0)
      {
        ++fNskippedNC;
        return;
      }

      // Fiducial cut: require the neutrino interaction vertex inside the active volume
      // (optionally inset by FiducialInset cm). This removes "rock"/external events
      // where the neutrino interacted in the surrounding cryostat/concrete/rock — those
      // deposit almost nothing in the active argon, so their residual is ~E_nu and they
      // are not detector events. (You'll see non-argon targets like O-16/Si-28 in their
      // vertex lists.) Skipped events fill no trees.
      {
        const double d = fFidInset;
        bool nuInFid = (fNuVtxX > fXmin + d && fNuVtxX < fXmax - d &&
                        fNuVtxY > fYmin + d && fNuVtxY < fYmax - d &&
                        fNuVtxZ > fZmin + d && fNuVtxZ < fZmax - d);
        if (!nuInFid)
        {
          ++fNskippedFid;
          return;
        }
      }

      // ---- (2) Geant4 particles: build maps ----
      art::Handle<std::vector<simb::MCParticle>> ph;
      if (!event.getByLabel(fSimLabel, ph))
        throw cet::exception("MyEnergyAnalysis") << "No MCParticles (" << fSimLabel.encode() << ")";

      for (auto const &p : *ph)
        fPmap[p.TrackId()] = &p;

      // mother -> children index
      std::map<int, std::vector<int>> children;
      for (auto const &p : *ph)
        children[p.Mother()].push_back(p.TrackId());
for (auto const &p : *ph)
{
  if (p.Mother() != 0)
    continue;

  int pdg = p.PdgCode();

  // Do not count incoming/primary neutrinos as visible particle-category energy.
  if (isNeutrino(pdg))
    continue;

  double trueKE = p.Momentum(0).E() - p.Mass();
  if (trueKE < 0)
    trueKE = 0.0;

  switch (depositCategory(pdg))
  {
  case 0:
    fE_true_mu += trueKE;
    break;
  case 1:
    fE_true_p += trueKE;
    break;
  case 2:
    fE_true_n += trueKE;
    break;
  case 3:
    fE_true_pi += trueKE;
    break;
  case 4:
    fE_true_em += trueKE;
    break;
  case 5:
    fE_true_nuc += trueKE;
    break;
  default:
    fE_true_other += trueKE;
    break;
  }
}

      // ---- (3) deposited (visible) energy from SimEnergyDeposit ----
      // SimEnergyDeposit::Energy() is in MeV; convert to GeV. Deposits exist only in
      // sensitive volumes, so this is already restricted to the active volume.
      art::Handle<std::vector<sim::SimEnergyDeposit>> edh;
      if (event.getByLabel(fEdepLabel, edh))
      {
        for (auto const &d : *edh)
        {
          double e = d.Energy() * 1e-3; // MeV -> GeV
          fE_dep_total += e;
          // empirical active-volume probe: track the spatial extent of deposits
          double dx = d.MidPointX(), dy = d.MidPointY(), dz = d.MidPointZ();
          fDepXmin = std::min(fDepXmin, dx);
          fDepXmax = std::max(fDepXmax, dx);
          fDepYmin = std::min(fDepYmin, dy);
          fDepYmax = std::max(fDepYmax, dy);
          fDepZmin = std::min(fDepZmin, dz);
          fDepZmax = std::max(fDepZmax, dz);
          ++fDepCount;
          int prim = rootPrimaryPDG(std::abs(d.TrackID()));
          switch (depositCategory(prim))
          {
          case 0:
            fE_dep_mu += e;
            break;
          case 1:
            fE_dep_p += e;
            break;
          case 2:
            fE_dep_n += e;
            break;
          case 3:
            fE_dep_pi += e;
            break;
          case 4:
            fE_dep_em += e;
            break;
          case 5:
            fE_dep_nuc += e;
            break;
          default:
            fE_dep_other += e;
            break;
          }
        }
      }
      else
      {
        mf::LogWarning("MyEnergyAnalysis")
            << "SimEnergyDeposit '" << fEdepLabel.encode()
            << "' not found — E_dep will be 0. Fix the label (see header).";
      }

      // ---- (4) interaction vertices: binding + neutrino, per channel ----
      for (auto const &parentPair : children)
      {
        int parentID = parentPair.first;
        if (parentID == 0)
          continue; // skip the primary-list pseudo-mother
        auto pit = fPmap.find(parentID);
        if (pit == fPmap.end())
          continue;
        const simb::MCParticle &parent = *pit->second;

        // cluster this parent's daughters by (creation position, process)
        std::vector<Vertex> verts;
        for (int childID : parentPair.second)
        {
          auto cit = fPmap.find(childID);
          if (cit == fPmap.end())
            continue;
          const simb::MCParticle *d = cit->second;
          const TLorentzVector &v = d->Position(0);
          const std::string proc = d->Process();
          bool merged = false;
          for (auto &vx : verts)
          {
            if (vx.process == proc &&
                std::abs(vx.x - v.X()) < 0.01 &&
                std::abs(vx.y - v.Y()) < 0.01 &&
                std::abs(vx.z - v.Z()) < 0.01)
            {
              vx.daughters.push_back(d);
              merged = true;
              break;
            }
          }
          if (!merged)
            verts.push_back(Vertex{v.X(), v.Y(), v.Z(), v.T(), proc, {d}});
        }

        for (auto const &vtx : verts)
        {
          if (!interestingProcesses().count(vtx.process))
            continue; // skip pure-EM ionization

          // --- tally products ---
          double sumM_out = 0, sumKE_out_massive = 0, sumE_out_massless = 0, sumE_nu = 0;
          double sumCreatedMesonM = 0; // rest mass of created products (baryon-number-0 mesons and
                                       //   massive leptons, plus baryon-antibaryon pair mass). Paid
                                       //   for by the projectile and given back downstream
                                       //   (E_dep/E_escape), so it is NOT nuclear binding — removed
                                       //   in the decompose below.
          int sumA = 0, sumZ = 0;
          int nN = 0, nP = 0, nG = 0, nNu = 0, nNuc = 0;
          int residualPDG = 0;
          bool hasElectron = false;
          fV_out_pdg.clear();
          fV_out_KE.clear();
          std::vector<int> freeHadrons; // for the channel label
          std::vector<int> fragments;   // every nuclear fragment (residual + any alpha/d/t/...)

          for (const simb::MCParticle *d : vtx.daughters)
          {
            int dpdg = d->PdgCode();
            double dM = robustMass(dpdg, d->Mass());
            double dE = d->Momentum(0).E();
            double dKE = dE - dM;

            fV_out_pdg.push_back(dpdg);
            fV_out_KE.push_back(dKE);

            sumM_out += dM;
            if (dM > 1e-9)
              sumKE_out_massive += dKE;
            else
              sumE_out_massless += dE;
            // Created meson (baryon number 0, massive): its rest mass is new — paid for by the
            // projectile — and it leaves the vertex to deposit or escape downstream. Accumulate
            // so the decompose below can take it back out of the binding number.
            // A created antibaryon signals a baryon-antibaryon PAIR (baryon-number conservation
            // requires a partner baryon in the final state): both partners are new mass paid for
            // by the projectile, so remove 2x the antibaryon mass. Knocked-out target nucleons
            // (no accompanying antibaryon) are NOT removed — their mass is handled by the A/Z
            // target balance, which the pair leaves untouched (net baryon number 0).
            if (isCreatedParticle(dpdg) && baryonNumber(dpdg) == 0 && dM > 1e-9)
              sumCreatedMesonM += dM;
            else if (baryonNumber(dpdg) < 0)
              sumCreatedMesonM += 2.0 * dM;

            if (isNeutrino(dpdg))
            {
              sumE_nu += dE;
              ++nNu;
            }
            else if (dpdg == 22)
              ++nG;
            else if (dpdg == 2112)
            {
              ++nN;
              sumA += 1;
              freeHadrons.push_back(dpdg);
            }
            else if (dpdg == 2212)
            {
              ++nP;
              sumA += 1;
              sumZ += 1;
              freeHadrons.push_back(dpdg);
            }
            else if (isNucleus(dpdg))
            {
              ++nNuc;
              fragments.push_back(dpdg);
              // Residual_PDG branch keeps the heaviest fragment (the recoil); the channel
              // key below records ALL fragments so each bucket is one exact reaction.
              if (residualPDG == 0 || nuclearA(dpdg) > nuclearA(residualPDG))
                residualPDG = dpdg;
              sumA += nuclearA(dpdg);
              sumZ += nuclearZ(dpdg);
            }
            else
            {
              sumA += baryonNumber(dpdg);
              sumZ += chargeOf(dpdg);
              if (std::abs(dpdg) == 11)
                hasElectron = true;
              else
                freeHadrons.push_back(dpdg);
            }
          }

          const double KE_in = incomingKEAtVertex(parent, vtx.x, vtx.y, vtx.z);
          const bool vtxInside = inside(vtx.x, vtx.y, vtx.z);

          // --- reconstruct the implicit target nucleus (A,Z balance) ---
          int targetPDG = 0;
          if (!isNucleus(parent.PdgCode()))
          {
            int At = sumA - baryonNumber(parent.PdgCode());
            int Zt = sumZ - chargeOf(parent.PdgCode());
            if (At > 0 && Zt >= 0)
              targetPDG = makeNucleusPDG(Zt, At);
          }

          // --- binding energy, two routes ---
          double Eb_mass = 0, Eb_cons = 0, Enu_vtx = 0, Eb_nuclear = 0, Eb_meson = 0;
          bool isMichel = false;

          if (isDecayProcess(vtx.process) || isCaptureAtRest(vtx.process))
          {
            // Conversion of the projectile (mu- decay -> e nu nu; mu- capture -> nu n).
            // No binding. The escaping neutrino is a loss ONLY if produced inside the
            // active volume — a muon that left the detector and stopped/decayed/captured
            // outside already had its energy counted as charged escape (incl. rest mass,
            // see the escape loop), so counting these neutrinos again would double count.
            Enu_vtx = sumE_nu;
            isMichel = (std::abs(parent.PdgCode()) == 13 && hasElectron && nNu >= 1);
            if (vtxInside && isMichel)
            {
              ++fN_michel_inside;
              fMichel_nu_E += sumE_nu;
            }
            // (E_neutrino_inside itself is summed globally in stage 6, which already
            //  applies the inside test, so nothing is added to the ledger here.)
          }
          else if (vtx.process == "hadElastic")
          {
            // Elastic scatter: the projectile survives (not in the daughter list), only
            // a nuclear recoil is created. No nucleon liberation => binding ~ 0, and the
            // A/Z target reconstruction is invalid (it would miss the surviving parent
            // and mis-identify the target by one nucleon). Record the channel, no binding.
            Eb_mass = 0.0;
            Eb_cons = 0.0;
            targetPDG = residualPDG; // recoil nucleus is (essentially) the target
          }
          else if (isNuclearTargetProcess(vtx.process))
          {
            // mass route (primary): ground-state mass change = -Q.
            // Target stays TABLE-ONLY: if the (almost always argon) target isn't in the
            // table, abandon the mass route (NaN) and use the conservation route, which
            // needs no target mass — mixing a crude target mass with table product
            // masses would corrupt the difference.
            double M_target = (targetPDG ? getMassFromPDG(targetPDG) : -1.0);
            double M_in = robustMass(parent.PdgCode(), parent.Mass());
            if (isNucleus(parent.PdgCode()))
            {
              Eb_mass = sumM_out - M_in; // incoming nucleus IS the target
            }
            else if (M_target > 0)
            {
              Eb_mass = sumM_out - M_in - M_target;
            }
            else
            {
              Eb_mass = std::nan(""); // target not in table
            }
            // conservation route (cross-check): energy not carried by product KE
            Eb_cons = KE_in - sumKE_out_massive - sumE_out_massless;

            // Pick the trustworthy TOTAL rest-mass change.
            // Eb_mass and Eb_cons both measure the total rest-mass change of the system; they
            // share no inputs (mass route: table only; conservation route: kinematics only), so
            // agreement is a genuine closure test — they matched to 0.1 MeV on the clean event.
            // Prefer the exact route when they agree; else the conservation route (no target-
            // table dependence); the fallback may be NaN and is caught by the guard below.
            const double kAgreeGeV = 0.05; // routes agree to 50 MeV
            bool agree = std::isfinite(Eb_mass) && std::isfinite(Eb_cons) &&
                         std::abs(Eb_mass - Eb_cons) < kAgreeGeV;
            double Eb_total; // total rest-mass change at this vertex
            if (agree)
              Eb_total = Eb_mass;
            else if (std::isfinite(Eb_cons))
              Eb_total = Eb_cons;
            else
              Eb_total = Eb_mass;
            //
            // sumCreatedMesonM is the OUTGOING created-meson mass. If the projectile is itself a
            // meson (pi/K inelastic), its mass was already in the initial state (subtracted via
            // M_in), so only the NET created meson mass should come out — subtract the incoming
            // meson mass too. For a nucleon/nucleus projectile there is no incoming meson.
            int parentPDG = parent.PdgCode();
            bool parentIsMeson = !isNucleus(parentPDG) && baryonNumber(parentPDG) == 0 &&
                                 !isNeutrino(parentPDG) && parentPDG != 22 &&
                                 std::abs(parentPDG) != 11 && std::abs(parentPDG) != 13 &&
                                 std::abs(parentPDG) != 15 &&
                                 robustMass(parentPDG, parent.Mass()) > 1e-9;
            double Mmeson_in = parentIsMeson ? robustMass(parentPDG, parent.Mass()) : 0.0;
            Eb_meson = sumCreatedMesonM - Mmeson_in; // net created-meson rest mass
            Eb_nuclear = std::isfinite(Eb_total) ? (Eb_total - Eb_meson) : std::nan("");

            // Guard on the nuclear remainder (not the total). It cannot physically exceed the
            // target's total binding (~344 MeV for Ar-40; 0.5 GeV is a generous ceiling). This
            // used to cap the total, which rejected every meson vertex and discarded the nuclear
            // piece with it. After the decompose a violation is rare and means the vertex is
            // mis-reconstructed (wrong target, missing product) — skip only those.
            const double kNucCapGeV = 0.5;
            double Eb_use = 0.0;
            if (std::isfinite(Eb_nuclear) && std::abs(Eb_nuclear) < kNucCapGeV)
            {
              Eb_use = Eb_nuclear;
            }
            else
            {
              // mf::LogWarning("MyEnergyAnalysis")
              //   << "Nuclear binding out of range after meson removal (total=" << Eb_total
              //   << " meson=" << Eb_meson << " nuclear=" << Eb_nuclear << " GeV) process "
              //   << vtx.process << " in_pdg " << parentPDG << " — skipped from sum.";
              Eb_use = 0.0;
            }
            fE_binding_total += (vtxInside ? Eb_use : 0.0);
            // Binding is only a detector loss if it happens inside the active volume.
            // A daughter that left the volume already had its KE counted as escape; its
            // subsequent interactions out in the cryostat are not the detector's loss.
            // (The vertex is still recorded in VertexTree with its Inside flag, so the
            //  per-channel table can include or exclude out-of-volume vertices.)
            // Positive E_binding = energy locked into rest mass (lost); negative =
            // exothermic release (e.g. nCapture), which reappears as gammas.
          }

          // --- build channel key: "n+Ar40[neutronInelastic]->Ar38+2n" ---
          auto nucName = [](int pdg) -> std::string
          {
            if (!isNucleus(pdg))
              return std::to_string(pdg);
            return "Z" + std::to_string(nuclearZ(pdg)) + "A" + std::to_string(nuclearA(pdg));
          };
          std::map<int, int> outCount;
          for (int h : freeHadrons)
            ++outCount[h];
          for (int g : fragments)
            ++outCount[g];
          std::string chan = std::to_string(parent.PdgCode());
          if (targetPDG)
            chan += "+" + nucName(targetPDG);
          chan += "[" + vtx.process + "]->";
          bool first = true;
          for (auto const &kv : outCount)
          {
            if (!first)
              chan += "+";
            if (kv.second > 1)
              chan += std::to_string(kv.second);
            chan += nucName(kv.first);
            first = false;
          }

          // Stamp every daughter of this vertex with the channel that produced it, so
          // escaping tracks (stage 5) can be attributed to their production interaction.
          for (const simb::MCParticle *d : vtx.daughters)
            birthChannel[d->TrackId()] = chan;

          // --- fill VertexTree ---
          fV_event = fEvent;
          fV_x = vtx.x;
          fV_y = vtx.y;
          fV_z = vtx.z;
          fV_t = vtx.t;
          fV_in_pdg = parent.PdgCode();
          fV_in_trk = parent.TrackId();
          fV_in_KE = KE_in;
          fV_process = vtx.process;
          fV_inside = vtxInside;
          fV_target_pdg = targetPDG;
          fV_residual_pdg = residualPDG;
          fV_nOut = (int)vtx.daughters.size();
          fV_nNeutron = nN;
          fV_nProton = nP;
          fV_nGamma = nG;
          fV_nNeutrino = nNu;
          fV_nNucleus = nNuc;
          fV_Ebind_mass = Eb_mass;
          fV_Ebind_cons = Eb_cons;
          fV_Ebind_diff = (std::isfinite(Eb_mass) ? Eb_mass - Eb_cons : std::nan(""));
          fV_Ebind_meson = Eb_meson;
          fV_Ebind_nuclear = Eb_nuclear;
          fV_Enu_vtx = Enu_vtx;
          fV_isMichel = isMichel;
          fV_channel = chan;
          fV_channel_hash = std::hash<std::string>{}(chan);
          fVertexTree->Fill();
        }
      }

      // ---- (5) per-track escape (charged + neutral) ----
      // Escape KE = KE at the LAST inside->outside crossing. Contained tracks (final
      // point inside) contribute 0. Neutrinos are excluded here — they are handled by
      // the neutrino term so they are not double counted.
      for (auto const &p : *ph)
      {
        int pdg = p.PdgCode();
        if (isNeutrino(pdg))
          continue;

        const unsigned int N = p.NumberTrajectoryPoints();
        if (N < 2)
          continue;

        bool wasInside = inside(p.Position(0));
        bool everInside = wasInside;
        double exitKE = -1, exitE = -1;
        unsigned int exitIdx = 0;
        for (unsigned int i = 1; i < N; ++i)
        {
          bool nowInside = inside(p.Position(i));
          if (wasInside && !nowInside)
          { // record LAST inside->outside
            exitKE = pointKE(p, i);
            exitE = p.Momentum(i).E();
            exitIdx = i;
          }
          if (nowInside)
            everInside = true;
          wasInside = nowInside;
        }
        if (!everInside)
          continue; // never in the detector at all
        if (exitKE < 0)
          continue; // ended inside -> contained
        if (exitKE < 0)
          exitKE = 0;

        // Energy that becomes invisible when the track leaves the active volume.
        // Created particles (mu, e, pi, gamma, ...) take their rest mass with them and
        // it came from the neutrino energy -> count total E. Pre-existing nucleons/nuclei
        // carry mass that was already in the target nucleus -> count KE only.
        double loss = isCreatedParticle(pdg) ? exitE : exitKE;

        bool charged = (chargeOf(pdg) != 0);
        if (charged)
          fE_escape_charged += loss;
        else
          fE_escape_neutral += loss;

        fX_event = fEvent;
        fX_trk = p.TrackId();
        fX_pdg = pdg;
        fX_birthKE = pointKE(p, 0);
        fX_exitKE = loss; // ExitKE branch now holds the loss
        fX_exitX = p.Position(exitIdx).X();
        fX_exitY = p.Position(exitIdx).Y();
        fX_exitZ = p.Position(exitIdx).Z();
        fX_endProcess = p.EndProcess();
        fX_charged = charged;
        {
          auto bc = birthChannel.find(p.TrackId());
          fX_birthChannel = (bc != birthChannel.end()) ? bc->second : std::string();
        }
        fEscapeTree->Fill();
      }

      // ---- (6) secondary (decay) neutrinos created inside the active volume ----
      // PRIMARY neutrinos (the GENIE outgoing NC neutrino) are counted from MCTruth in
      // stage 1; some largeant configs also keep a primary-neutrino stub here, so skip
      // Mother()==0 neutrinos to avoid double counting. Decay/Michel neutrinos (Mother!=0)
      // are not in the generator truth and are counted here.
      for (auto const &p : *ph)
      {
        if (!isNeutrino(p.PdgCode()))
          continue;
        if (p.Mother() == 0)
          continue; // primary neutrino -> handled via MCTruth
        if (inside(p.Vx(), p.Vy(), p.Vz()))
          fE_neutrino_inside += p.Momentum(0).E();
      }

      // ---- (7) close the ledger; residual is the QA number at truth level ----
      fE_residual = fGen_nu_E - (fE_dep_total + fE_binding_total + fE_escape_neutral + fE_escape_charged + fE_neutrino_inside);

      // Containment flag: fraction of the neutrino energy that left the active volume.
      // Not a cut — stored so the analysis can threshold offline (e.g. EscapeFrac < 0.3
      // for a well-contained subset). High-energy DIS events leak large fractions.
      fEscapeFrac = (fGen_nu_E > 0)
                        ? (fE_escape_charged + fE_escape_neutral) / fGen_nu_E
                        : 0.0;

      fEventTree->Fill();
    }

    DEFINE_ART_MODULE(MyEnergyAnalysis)

  } // namespace example
} // namespace lar

// share instructions for grids
// estimator for neutrons learn about neutron energy. Look at the ones that leave the detector and understand what they to and find the energy of the interactions.
// when neutron comes in look at the highest energy and plot the difference
// look at the variables for plots