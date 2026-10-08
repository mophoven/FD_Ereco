// Plot corrected neutral escape quantities from the per-file diagnostic CSV.
// root -l -b -q 'PlotCorrectedNeutralEscape.C("neutral_escape_diagnostic_all_files.csv")'

#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::string> splitCSV(const std::string& line) {
  std::vector<std::string> out; std::string field; bool quoted=false;
  for(size_t i=0;i<line.size();++i){
    char c=line[i];
    if(c=='"'){if(quoted&&i+1<line.size()&&line[i+1]=='"'){field+='"';++i;}else quoted=!quoted;}
    else if(c==','&&!quoted){out.push_back(field);field.clear();}
    else field+=c;
  }
  out.push_back(field); return out;
}
}

void PlotCorrectedNeutralEscape(const char* csv="neutral_escape_diagnostic_all_files.csv") {
  std::ifstream in(csv); if(!in){std::cerr<<"ERROR: cannot open "<<csv<<'\n';return;}
  std::string line; std::getline(in,line);
  TH2D hNeutral("hNeutral",
    "Neutral Escape Energy vs Highest-Energy Escaping Neutron;Max escaping neutron ExitKE [GeV];E_{escape,neutral} [GeV]",
    100,0,5,100,0,5);
  TH2D hDifference("hDifference",
    "Neutral Escape Energy Minus Leading Neutron ExitKE;Max escaping neutron ExitKE [GeV];E_{escape,neutral} - Max neutron ExitKE [GeV]",
    100,0,5,100,0,5);
  long long plotted=0;
  while(std::getline(in,line)){
    auto f=splitCSV(line); if(f.size()<18||f[9]=="nan")continue;
    const double neutral=std::stod(f[8]);
    const double leading=std::stod(f[9]);
    const double difference=std::stod(f[11]);
    hNeutral.Fill(leading,neutral);
    hDifference.Fill(leading,difference);
    ++plotted;
  }
  if(!plotted){std::cerr<<"ERROR: no usable rows\n";return;}
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.045);
  gStyle->SetPalette(kViridis);
  auto draw=[&](TH2D& h,const char* name){
    TCanvas c(name,h.GetTitle(),1100,850);
    c.SetLeftMargin(0.13); c.SetRightMargin(0.15); c.SetBottomMargin(0.12); c.SetTopMargin(0.10);
    h.GetXaxis()->SetTitleSize(0.045); h.GetYaxis()->SetTitleSize(0.045);
    h.GetXaxis()->SetLabelSize(0.038); h.GetYaxis()->SetLabelSize(0.038);
    h.GetYaxis()->SetTitleOffset(1.30); h.Draw("COLZ");
    c.SaveAs((std::string(name)+".png").c_str());
  };
  draw(hNeutral,"neutral_escape_vs_leading_neutron_corrected");
  draw(hDifference,"neutral_escape_minus_leading_neutron_corrected");
  std::cout<<"Plotted "<<plotted<<" events with escaping neutrons as corrected TH2D density maps.\n";
}
