#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TGLViewer.h>
#include <TVirtualViewer3D.h>

TFile* file = TFile::Open("test5files.root");
TTree* tree = (TTree*)file->Get("MyEnergyAnalysis/MyTree");

void ParticleTracker(){


int Event = 414;

bool desiredEventProcessed = false;

TTreeReader reader(tree);

TTreeReaderValue<std::vector<float>> xPosReader(reader, "P_vtx_x");
TTreeReaderValue<std::vector<float>> yPosReader(reader, "P_vtx_y");
TTreeReaderValue<std::vector<float>> zPosReader(reader, "P_vtx_z");

while (reader.Next()){
        if (!reader.Next()) break;
        if (reader.GetCurrentEntry() == Event && !desiredEventProcessed){
                desiredEventProcessed = true;

        const std::vector<float>& xPos = *xPosReader;
        const std::vector<float>& yPos = *yPosReader;
        const std::vector<float>& zPos = *zPosReader;

        TCanvas* canvas = new TCanvas("canvas", "Particle Tracker", 800, 600);
        TVirtualPad* myPad = canvas->GetPad(1);
        TGLViewer viewer(myPad);

        canvas -> SetViewer3D(&viewer);

        TPolyLine3D* nutrack = new TPolyLine3D(2);
                nutrack->SetPoint(0, xPos[0], yPos[0], zPos[0]);
                nutrack->SetPoint(1, xPos[2], yPos[2], zPos[2]);
                nutrack->SetLineColor(kBlue);
                nutrack->SetLineWidth(2);
                nutrack->Draw();

        TPolyLine3D* neutrontrack = new TPolyLine3D(2);
                neutrontrack->SetPoint(0, xPos[2], yPos[2],zPos[2]);
                neutrontrack->SetPoint(1, xPos[16], yPos[16], zPos[16]);
                neutrontrack->SetLineColor(5);
                neutrontrack->SetLineWidth(2);
                neutrontrack->Draw();

        TPolyLine3D* p1 = new TPolyLine3D(2);
                p1->SetPoint(0, xPos[1], yPos[1], zPos[1]);
                p1->SetPoint(1, xPos[4], yPos[4], zPos[4]);
                p1->SetLineColor(kRed);
                p1->SetLineWidth(2);
                p1->Draw();

        TPolyLine3D* p2 = new TPolyLine3D(2);
                p2->SetPoint(0, xPos[2], yPos[2], zPos[2]);
                p2->SetPoint(1, xPos[5], yPos[5], zPos[5]);
                p2->SetLineColor(kRed);
                p2->SetLineWidth(2);
                p2->Draw();

        TPolyLine3D* p3 = new TPolyLine3D(2);
                p3->SetPoint(0, xPos[5], yPos[5], zPos[5]);
                p3->SetPoint(1, xPos[6], yPos[6], zPos[6]);
                p3->SetLineColor(kRed);
                p3->SetLineWidth(2);
                p3->Draw();

        TPolyLine3D* p4 = new TPolyLine3D(2);
                p4->SetPoint(0, xPos[5], yPos[5], zPos[5]);
                p4->SetPoint(1, xPos[7], yPos[7], zPos[7]);
                p4->SetLineColor(kRed);
                p4->SetLineWidth(2);
                p4->Draw();

        TPolyLine3D* p5 = new TPolyLine3D(2);
                p5->SetPoint(0, xPos[6], yPos[6], zPos[6]);
                p5->SetPoint(1, xPos[8], yPos[8], zPos[8]);
                p5->SetLineColor(5);
                p5->SetLineWidth(2);
                p5->Draw();

        TPolyLine3D* p6 = new TPolyLine3D(2);
                p6->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p6->SetPoint(1, xPos[9], yPos[9], zPos[9]);
                p6->SetLineColor(kRed);
                p6->SetLineWidth(2);
                p6->Draw();

        TPolyLine3D* p7 = new TPolyLine3D(2);
                p7->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p7->SetPoint(1, xPos[10], yPos[10], zPos[10]);
                p7->SetLineColor(kRed);
                p7->SetLineWidth(2);
                p7->Draw();

        TPolyLine3D* p8 = new TPolyLine3D(2);
                p8->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p8->SetPoint(1, xPos[11], yPos[11], zPos[11]);
                p8->SetLineColor(kRed);
                p8->SetLineWidth(2);
                p8->Draw();

        TPolyLine3D* p9 = new TPolyLine3D(2);
                p9->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p9->SetPoint(1, xPos[12], yPos[12], zPos[12]);
                p9->SetLineColor(kRed);
                p9->SetLineWidth(2);
                p9->Draw();

        TPolyLine3D* p10 = new TPolyLine3D(2);
                p10->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p10->SetPoint(1, xPos[13], yPos[13], zPos[13]);
                p10->SetLineColor(kRed);
                p10->SetLineWidth(2);
                p10->Draw();

        TPolyLine3D* p11 = new TPolyLine3D(2);
                p11->SetPoint(0, xPos[8], yPos[8], zPos[8]);
                p11->SetPoint(1, xPos[14], yPos[14], zPos[14]);
                p11->SetLineColor(kRed);
                p11->SetLineWidth(2);
                p11->Draw();

        TPolyLine3D* p12 = new TPolyLine3D(2);
                p12->SetPoint(0, xPos[7], yPos[7], zPos[7]);
                p12->SetPoint(1, xPos[15], yPos[15], zPos[15]);
                p12->SetLineColor(6);
                p12->SetLineWidth(2);
                p12->Draw();

        canvas->Draw();
        }
}
}
