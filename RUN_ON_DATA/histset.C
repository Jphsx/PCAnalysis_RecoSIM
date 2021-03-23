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
/////////////////////////////////////////?SIM things
	//get sim pc mask
/*	sim_pc SPC;
	SPC = GetSimPC(s);

	int np14 = SPC.p14_key.size();
	FillTH1(id_numSPCHist, np14, w);
	//loop over sim
	int gidx,t1idx,t2idx;
	for(int i=0; i<np14; i++ ){
		gidx = SPC.p14_key[i];
		
		FillTH1(id_ptSPCHist,SimTrk_pt[gidx], w);	
		
	}

	double sx,sy,sz;
	double ptp,ptm;
	for(int i=0; i<nSimVtx; i++){//loop because mask
	//fill denoms based on acceptance cuts
		//does this ptype14 pass acceptance?
		//acceptance => |zpc|<25cm |cost_g|<0.85 min(pt1,pt2)>0.2GeV
		//matching is considered within 0.5cm in vtxdl
		if( SPC.sim_mask[i] == 14 ){
			gidx = SPC.p14_g[i];
			t1idx = SPC.p14_t1[i];
			t2idx = SPC.p14_t2[i];
			if(abs(SimVtx_z[i]) < 25 ){
			if(SimTrk_pt[t1idx]>0.2 && SimTrk_pt[t2idx]>0.2){
			//cos0.85 => eta 1.25615 approx
			if(abs(SimTrk_eta[gidx]) < 1.25615){ 		
			//restrict to R <20
			sx = SimVtx_x[i];
			sy = SimVtx_y[i];
			sz = SimVtx_z[i];
			if( sqrt(sx*sx + sy*sy ) < 20. ){
			//	
				//this sim conv has passed acceptance
				FillTH1(id_effPtD, SimTrk_pt[gidx], w);
				FillTH1(id_effRD, sqrt(sx*sx + sy*sy ) , w);
				
				if(SimTrk_pdgId[t1idx] > 0){
					ptm = SimTrk_pt[t1idx];
					ptp = SimTrk_pt[t2idx];
				}
				else{
					ptm = SimTrk_pt[t2idx];
					ptp = SimTrk_pt[t1idx];
				}
	
				FillTH1(id_effXPD, ptp/(ptm+ptp), w);
			}}}}
		}
	}
	
//////////////////////////////////////////end SIM THINGS
	//calculate common variables
	std::vector<CommonVars> CVs = GetCommonVars(s,false);//note this is mc
	 	
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


	////////////More sim stuff (matching for efficiency and purity)
	int pcidx;
	int matchidx;
	for(int i=0; i<HGN.vsel.size(); i++){
		pcidx = HGN.vsel[i];
		//if we have a match
		if(Conv_vtxdl[pcidx] < 0.5){
		matchidx = Conv_convVtxIdx[pcidx];
		if(SimVtx_processType[matchidx] == 14){		
			FillTH1(id_effPtN, CVs[pcidx].pt, w);
			FillTH1(id_effRN, CVs[pcidx].radius, w);
			FillTH1(id_effXPN, CVs[pcidx].xplus, w);
			FillTH1(id_purityPtN, CVs[pcidx].pt, w);
			FillTH1(id_purityRN, CVs[pcidx].radius, w);
			FillTH1(id_purityXPN, CVs[pcidx].xplus, w);

			FillTH1(id_nconvPt, CVs[pcidx].pt, w);
			FillTH1(id_nconvR, CVs[pcidx].radius, w);
			FillTH1(id_nconvXP, CVs[pcidx].xplus, w);
		}}
		FillTH1(id_purityPtD, CVs[pcidx].pt, w);
		FillTH1(id_purityRD, CVs[pcidx].radius, w);
		FillTH1(id_purityXPD, CVs[pcidx].xplus, w);

		if(Conv_vtxdl[pcidx] > 0.5){
			FillTH1(id_nconvPt_fake, CVs[pcidx].pt, w);
			FillTH1(id_nconvR_fake, CVs[pcidx].radius, w);
			FillTH1(id_nconvXP_fake, CVs[pcidx].xplus, w);
		}
		if(Conv_vtxdl[pcidx] < 0.5){
		matchidx = Conv_convVtxIdx[pcidx];
		if(SimVtx_processType[matchidx] != 14){
			FillTH1(id_nconvPt_fake, CVs[pcidx].pt, w);
			FillTH1(id_nconvR_fake, CVs[pcidx].radius, w);
			FillTH1(id_nconvXP_fake, CVs[pcidx].xplus, w);
		}}	
		
	}

	/////end matching
*/

////////////////////////DATA RUN//////////
	//calculate common variables
	std::vector<CommonVars> CVs = GetCommonVars(s,true);//note this is mc
	 	
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
	int pcidx;
        int matchidx;
        for(int i=0; i<HGN.vsel.size(); i++){
                pcidx = HGN.vsel[i];
	      if(CVs[pcidx].radius < 20){//restrict radius for now
	       FillTH1(id_nconvPt, CVs[pcidx].pt, w);
               FillTH1(id_nconvR, CVs[pcidx].rho, w);
               FillTH1(id_nconvXP, CVs[pcidx].xplus, w); 	
		}
	}
	

}  // End of sub-program
#endif
