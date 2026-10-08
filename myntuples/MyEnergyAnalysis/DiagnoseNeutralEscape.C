// Diagnose the event-level neutral escape bookkeeping.
//
// Run with:
//   root -l -b -q 'DiagnoseNeutralEscape.C("MiloEnergy_merged_500.root",20)'
//
// Important: in the current module, EscapeTree::ExitKE is filled with the
// variable named "loss".  For neutrons that value is their exit kinetic
// energy, but for some other particle species it can have different semantics.

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
TTree *findDiagnosticTree(TDirectory *dir, const char *wanted)
{
  if (!dir) return nullptr;
  if (auto *tree = dynamic_cast<TTree *>(dir->Get(wanted))) return tree;

  TIter next(dir->GetListOfKeys());
  while (auto *key = dynamic_cast<TKey *>(next())) {
    if (std::string(key->GetClassName()).find("TDirectory") == std::string::npos)
      continue;
    if (auto *tree = findDiagnosticTree(
            dynamic_cast<TDirectory *>(key->ReadObj()), wanted))
      return tree;
  }
  return nullptr;
}

struct EventRow {
  Long64_t entry = -1;
  int run = 0;
  int subrun = 0;
  int event = 0;
  double neutralLedger = 0.0;
};

struct EscapeRow {
  Long64_t entry = -1;
  int event = 0;
  int track = 0;
  int pdg = 0;
  double birthKE = 0.0;
  double exitField = 0.0;
  bool charged = false;
};
}

void DiagnoseNeutralEscape(const char *fileName, int maxEventEntries = 20)
{
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "ERROR: cannot open " << fileName << '\n';
    return;
  }

  TTree *eventTree = findDiagnosticTree(&input, "EventTree");
  TTree *escapeTree = findDiagnosticTree(&input, "EscapeTree");
  if (!eventTree || !escapeTree) {
    std::cerr << "ERROR: file must contain EventTree and EscapeTree.\n";
    return;
  }

  // Read all EventTree IDs so duplicate IDs in a merged file can be detected.
  int event = 0, run = 0, subrun = 0;
  double neutralLedger = 0.0;
  eventTree->SetBranchAddress("Event", &event);
  eventTree->SetBranchAddress("Run", &run);
  eventTree->SetBranchAddress("SubRun", &subrun);
  eventTree->SetBranchAddress("E_escape_neutral", &neutralLedger);

  std::map<int, int> eventIdMultiplicity;
  std::vector<EventRow> selectedEvents;
  const Long64_t nEventEntries = eventTree->GetEntries();
  for (Long64_t i = 0; i < nEventEntries; ++i) {
    eventTree->GetEntry(i);
    ++eventIdMultiplicity[event];
    if ((int)selectedEvents.size() < maxEventEntries)
      selectedEvents.push_back({i, run, subrun, event, neutralLedger});
  }

  // Read every EscapeTree contributor.  This tree currently has Event only,
  // not Run/SubRun, so duplicate Event IDs make the association ambiguous.
  int xEvent = 0, track = 0, pdg = 0;
  double birthKE = 0.0, exitField = 0.0;
  bool charged = false;
  escapeTree->SetBranchAddress("Event", &xEvent);
  escapeTree->SetBranchAddress("TrackID", &track);
  escapeTree->SetBranchAddress("PDG", &pdg);
  escapeTree->SetBranchAddress("BirthKE", &birthKE);
  escapeTree->SetBranchAddress("ExitKE", &exitField);
  escapeTree->SetBranchAddress("Charged", &charged);

  std::map<int, std::vector<EscapeRow>> escapeByEventId;
  for (Long64_t i = 0; i < escapeTree->GetEntries(); ++i) {
    escapeTree->GetEntry(i);
    escapeByEventId[xEvent].push_back(
        {i, xEvent, track, pdg, birthKE, exitField, charged});
  }

  std::ofstream summary("neutral_escape_diagnostic_first20.csv");
  summary << "event_tree_entry,run,subrun,event,event_id_multiplicity,ambiguous_event_id,"
             "E_escape_neutral,max_neutron_ExitKE_field,max_neutron_track,"
             "E_escape_neutral_minus_max_neutron,sum_neutral_ExitKE_fields,"
             "E_escape_neutral_minus_sum_fields,n_neutral_escape_rows,n_neutron_escape_rows\n";

  std::ofstream contributors("neutral_escape_contributors_first20.csv");
  contributors << "event_tree_entry,run,subrun,event,ambiguous_event_id,"
                  "escape_tree_entry,track,pdg,birthKE,ExitKE_field,charged\n";

  int negativeCount = 0;
  int negativeUnambiguousCount = 0;
  int ambiguousCount = 0;
  constexpr double tolerance = 1.0e-9;

  std::cout << std::fixed << std::setprecision(9);
  std::cout << "\nChecking the first " << selectedEvents.size()
            << " EventTree entries\n";
  std::cout << "NOTE: EscapeTree::ExitKE is the module's stored 'loss' field.\n\n";

  for (const EventRow &ev : selectedEvents) {
    const bool ambiguous = eventIdMultiplicity[ev.event] > 1;
    if (ambiguous) ++ambiguousCount;

    double maxNeutronExit = -std::numeric_limits<double>::infinity();
    int maxNeutronTrack = -1;
    double sumNeutralFields = 0.0;
    int nNeutral = 0;
    int nNeutron = 0;

    const auto found = escapeByEventId.find(ev.event);
    if (found != escapeByEventId.end()) {
      for (const EscapeRow &x : found->second) {
        if (x.charged) continue;
        ++nNeutral;
        sumNeutralFields += x.exitField;
        if (x.pdg == 2112) {
          ++nNeutron;
          if (x.exitField > maxNeutronExit) {
            maxNeutronExit = x.exitField;
            maxNeutronTrack = x.track;
          }
        }
        contributors << ev.entry << ',' << ev.run << ',' << ev.subrun << ','
                     << ev.event << ',' << (ambiguous ? 1 : 0) << ','
                     << x.entry << ',' << x.track << ',' << x.pdg << ','
                     << std::setprecision(17) << x.birthKE << ','
                     << x.exitField << ',' << (x.charged ? 1 : 0) << '\n';
      }
    }

    const bool hasNeutron = maxNeutronTrack >= 0;
    const double difference = hasNeutron
                                  ? ev.neutralLedger - maxNeutronExit
                                  : std::numeric_limits<double>::quiet_NaN();
    const double sumDifference = ev.neutralLedger - sumNeutralFields;
    const bool negative = hasNeutron && difference < -tolerance;
    if (negative) {
      ++negativeCount;
      if (!ambiguous) ++negativeUnambiguousCount;
    }

    std::cout << "entry=" << ev.entry
              << "  run/subrun/event=" << ev.run << '/' << ev.subrun << '/'
              << ev.event
              << "  Eneutral=" << ev.neutralLedger;
    if (hasNeutron)
      std::cout << "  maxNExit=" << maxNeutronExit
                << " (track " << maxNeutronTrack << ')'
                << "  diff=" << difference;
    else
      std::cout << "  no escaping neutron";
    std::cout << "  neutralSum=" << sumNeutralFields
              << "  ledger-sum=" << sumDifference;
    if (negative) std::cout << "  <-- NEGATIVE";
    if (ambiguous)
      std::cout << "  <-- AMBIGUOUS: Event ID appears "
                << eventIdMultiplicity[ev.event] << " times";
    std::cout << '\n';

    summary << ev.entry << ',' << ev.run << ',' << ev.subrun << ','
            << ev.event << ',' << eventIdMultiplicity[ev.event] << ','
            << (ambiguous ? 1 : 0) << ',' << std::setprecision(17)
            << ev.neutralLedger << ',';
    if (hasNeutron)
      summary << maxNeutronExit << ',' << maxNeutronTrack << ',' << difference;
    else
      summary << "nan,-1,nan";
    summary << ',' << sumNeutralFields << ',' << sumDifference << ','
            << nNeutral << ',' << nNeutron << '\n';
  }

  std::cout << "\nSummary\n"
            << "  selected EventTree entries: " << selectedEvents.size() << '\n'
            << "  entries with ambiguous Event-only matching: " << ambiguousCount << '\n'
            << "  negative E_escape_neutral - max neutron ExitKE: "
            << negativeCount << '\n'
            << "  negative and unambiguous: " << negativeUnambiguousCount << '\n'
            << "\nCreated:\n"
            << "  neutral_escape_diagnostic_first20.csv\n"
            << "  neutral_escape_contributors_first20.csv\n";

  if (ambiguousCount > 0) {
    std::cout << "\nIMPORTANT: EscapeTree lacks Run/SubRun, so duplicated Event IDs "
                 "cannot be matched safely after hadd. Add Run and SubRun branches "
                 "to EscapeTree before interpreting those rows.\n";
  }
}
