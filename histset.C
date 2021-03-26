#ifndef HISTS
#define HISTS
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "recosim.C"
#include "PCTools.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <map>
#include <tuple>
#include <iomanip>
#include "Hungarian.h"
using MyTH1D = ROOT::TThreadedObject<TH1D>;
using MyTH2D = ROOT::TThreadedObject<TH2D>;

// struct for derived quantities of each conversion


class histset{
    public:
       double PI =4.0*atan(1.0);
       histset();
       void init();
       void setweightoption();
       void AnalyzeEntry(recosim& s);
       #include "Enums.h"
// make a big vector and load enumerated histograms onto the vector
       std::vector<MyTH1D*>  TH1Manager{};
       std::vector<MyTH2D*>  TH2Manager{};
// locate the histogram and perform pointer copying
       void FillTH1(int index, double x, double w);
       void FillTH2(int index, double x, double y, double w);
       void WriteHist();
};

histset::histset(){
    std::vector<MyTH1D*>  Manager1(numTH1Hist);
    TH1Manager=Manager1;
    std::vector<MyTH2D*>  Manager2(numTH2Hist);
    TH2Manager=Manager2;
    init();
    setweightoption();
}

void histset::setweightoption(){
    for(int i=0; i<numTH1Hist; i++){
        auto hptr = TH1Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }
    for(int i=0; i<numTH2Hist; i++){
        auto hptr = TH2Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }
}

#include "Hists.h"     //Put the histogram declarations in one separate file

void histset::FillTH1(int index, double x, double w=1.0){
	//we must make pointer copies for performance reasons when trying to fill a histogram
	auto myhist = TH1Manager.at(index)->Get();
	myhist->Fill(x,w);
}

void histset::FillTH2(int index, double x, double y, double w=1.0){
	auto myhist = TH2Manager.at(index)->Get();
	myhist->Fill(x,y,w);
}

void histset::WriteHist(){

	TFile* outfile = new TFile("Outfile.root", "RECREATE");

	for(int i=0; i<numTH1Hist; i++){
	//do a check for entries, merge isn't safe for empty histograms
        auto hptr = TH1Manager.at(i)->Get();
	   if(hptr->GetEntries() > 0){
           auto histmerged = TH1Manager.at(i)->Merge();
           TH1D* h = (TH1D*) histmerged->Clone();
		   outfile->WriteObject(h, h->GetName() );
        }
        else{
           auto h = TH1Manager.at(i)->Get()->Clone();
           outfile->WriteObject(h, h->GetName() );
        }
	}

	for(int i=0; i<numTH2Hist; i++){
	//do a check for entries, merge isn't safe for empty histograms
        auto hptr = TH2Manager.at(i)->Get();
	    if(hptr->GetEntries() > 0){
           auto histmerged = TH2Manager.at(i)->Merge();
           TH2D* h = (TH2D*) histmerged->Clone();
		   outfile->WriteObject(h, h->GetName() );
        }
        else{
           auto h = TH2Manager.at(i)->Get()->Clone();
           outfile->WriteObject(h, h->GetName() );
        }
	}
	outfile->Close();
}

void histset::AnalyzeEntry(recosim& s){
	double w = 1.0;
	#include "localTreeMembers.h"     //All the variable incantations needed

	FillTH1( id_numpcHist, numberOfPC,w);

	for(int i=0; i<numberOfPC; i++){
		double PC_pt = sqrt(PC_Px[i]*PC_Px[i] + PC_Py[i] * PC_Py[i]);
		FillTH1( id_pzHist, PC_Pz[i],w);
		FillTH1( id_ptHist, PC_pt,w);
		FillTH2( id_xyHist, PC_x[i], PC_y[i], w);
		FillTH2( id_xywideHist, PC_x[i],PC_y[i], w);
		double PC_r = sqrt(PC_x[i]*PC_x[i] + PC_y[i]*PC_y[i]);
		double PC_phi = atan2(PC_y[i],PC_x[i]);	
		FillTH2(id_rphiHist,PC_r, PC_phi,w);
	}

	//get sim pc mask
	sim_pc SPC;
	SPC = GetSimPC(s);

	int np14 = SPC.p14_key.size();
	FillTH1(id_numSPCHist, np14, w);
	//loop over sim
	int gidx;
	for(int i=0; i<np14; i++ ){
		gidx = SPC.p14_key[i];
		FillTH1(id_ptSPCHist,SimTrk_pt[gidx], w);	

	}
	
	//calculate common variables
	std::vector<CommonVars> CVs = GetCommonVars(s,false);
	 	
	//get cut mask (check numerator cuts)
	std::vector<bool> cutmask = GetCutMask(s,CVs);
	int npcCut=0;
	for(int i=0; i<cutmask.size(); i++){
		if(cutmask[i]) npcCut++;
	}
	FillTH1(id_numpccutHist, npcCut, w);
	
        //disambiguation (only do hungarian algorithm on pc's that pass cutmask)
	hgn_pc HGN = pc_disambiguation(s,cutmask);
	FillTH1( id_numHGNPCHist, HGN.vsel.size(), w);


}  // End of sub-program
#endif
