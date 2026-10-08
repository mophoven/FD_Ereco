// Plot corrected neutral escape quantities from the per-file diagnostic CSV.
// root -l -b -q 'PlotCorrectedNeutralEscape.C("neutral_escape_diagnostic_all_files.csv")'

#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TStyle.h>
#include <algorithm>
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
  std::vector<double> leading,neutral,difference;
  while(std::getline(in,line)){
    auto f=splitCSV(line); if(f.size()<18||f[9]=="nan")continue;
    neutral.push_back(std::stod(f[8])); leading.push_back(std::stod(f[9]));
    difference.push_back(std::stod(f[11]));
  }
  if(leading.empty()){std::cerr<<"ERROR: no usable rows\n";return;}
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.045);

  auto draw=[&](const char* name,const char* title,const char* ytitle,const std::vector<double>& y,double ymin,double ymax){
    TCanvas c(name,title,1100,760); c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetBottomMargin(0.12);
    double xmax=*std::max_element(leading.begin(),leading.end())*1.04;
    TH2D frame("frame",title,100,0,xmax,100,ymin,ymax);
    frame.GetXaxis()->SetTitle("Leading escaping neutron energy (GeV)"); frame.GetYaxis()->SetTitleOffset(1.35);
    frame.GetYaxis()->SetTitle(ytitle);
    frame.Draw(); TGraph g(leading.size(),leading.data(),y.data()); g.SetMarkerStyle(20); g.SetMarkerSize(0.35); g.SetMarkerColorAlpha(kAzure+2,0.38); g.Draw("P SAME");
    TLine zero(0,0,xmax,0); zero.SetLineColor(kRed+1); zero.SetLineWidth(2); zero.SetLineStyle(2); zero.Draw();
    c.SaveAs((std::string(name)+".png").c_str()); c.SaveAs((std::string(name)+".pdf").c_str());
  };
  double neutralMax=*std::max_element(neutral.begin(),neutral.end())*1.06;
  double diffMax=*std::max_element(difference.begin(),difference.end())*1.06;
  draw("neutral_escape_vs_leading_neutron_corrected","Corrected neutral escape energy","Total escaping neutral energy (GeV)",neutral,0,neutralMax);
  draw("neutral_escape_minus_leading_neutron_corrected","Corrected neutral escape minus leading neutron energy","Neutral escape minus leading neutron energy (GeV)",difference,0,diffMax);
  std::cout<<"Plotted "<<leading.size()<<" events with escaping neutrons.\n";
}
