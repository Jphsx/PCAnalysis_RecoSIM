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
	
	double w; 
	#include "localTreeMembers.h"     //All the variable incantations needed
//	double mc2data_w = 0.883927;
//	double npvweight = weight;
//	w = npvweight * mc2data_w;
//	w=12.05311;
	w=1.;
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

///////////////////////////////////////sim photons
//
	sim_g SG;
	SG =GetSimG(s);
	//s_g vars
	double sg_pt,sg_eta,sg_phi, sg_P, sg_gx, sg_gy, sg_gz, sg_Rxy_origin, sg_ex, sg_ey, sg_ez, sg_Rxy_endpt;

	int pc_eidx, pc_pidx;
	int exyz_ptype;
	//min pt acceptance cut
	double g_minpt = 0.2;//GeV
	//nbins = 40 from 0 to 20
	int nGbins = 20.;
	double gBinEnd = 20.;
	double gBinWidth = gBinEnd/ ((double)nGbins);
	double gBinStart = 0.;
	double binLow,bincenter;
	binLow = gBinStart;

	Float_t Rbins[] = { 1,5,9,13,18,20 };
	Int_t  binnum = 5;
	for(int i=0; i< SG.g_idx.size(); i++){	
		sg_pt = SimTrk_pt[ SG.g_idx[i] ];
		sg_eta = SimTrk_eta[ SG.g_idx[i] ];
		sg_phi = SimTrk_phi[ SG.g_idx[i] ];
		sg_P = sg_pt*cosh(sg_eta);
		sg_gx = SimVtx_x[ SG.parent_simvtx_idx[i] ];
		sg_gy = SimVtx_y[ SG.parent_simvtx_idx[i] ];
		sg_gz = SimVtx_z[ SG.parent_simvtx_idx[i] ];	
		sg_Rxy_origin = sqrt( sg_gx*sg_gx + sg_gy*sg_gy );
		if( SG.endpoint_vtx_idx[i] != -1){
		sg_ex = SimVtx_x[ SG.endpoint_vtx_idx[i] ];
		sg_ey = SimVtx_y[ SG.endpoint_vtx_idx[i] ];
		sg_ez = SimVtx_z[ SG.endpoint_vtx_idx[i] ];
		exyz_ptype = SimVtx_processType[ SG.endpoint_vtx_idx[i] ];
		sg_Rxy_endpt = sqrt( sg_ex*sg_ex + sg_ey*sg_ey );
		}else{
			sg_ex=0;
			sg_ey=0;
			sg_ez=0;
			sg_Rxy_endpt = -1;
		}
		
		//look at photons with *no* cuts
		if(sg_Rxy_origin < 0.5){
		if(abs( sg_gz ) < 1.){
		//highly transverse
		if(abs(sg_eta) < 0.5){
			
			FillTH1(id_Ng_BP1,0.5,w);
			FillTH1(id_Ng_BP2,0.5,w);
			//does this photon convert in bpix 1 or 2??
			//BPIX1 1,5   BPIX2 5,9	
			if( sg_Rxy_endpt > 1 && sg_Rxy_endpt < 5 && exyz_ptype == 14){
				FillTH1(id_Nc_BP1,0.5,w);
			}
			if( sg_Rxy_endpt > 5 && sg_Rxy_endpt < 9 && exyz_ptype == 14){
				FillTH1(id_Nc_BP2,0.5,w);
			}
		
		}}}
		
//
		//does this sim photon satisfy acceptance cuts?
		//min energy?
		if(sg_pt > 2*sqrt( g_minpt*g_minpt + GV.MASS_ELECTRON * GV.MASS_ELECTRON ) ){
		//origin less than z<25?
		if( abs(sg_gz) < 25.){
		//cos theta< 0.85?
		if( abs( sg_eta ) < 1.25615 ){
		//if( sg_Rxy_origin < 0.5 ){
			//now check if the photon is created before radial point and ends after (if it has an endpoint)
			//if( sg_Rxy_origin < 1.0 ){
			//	FillTH1(id_ngPrompt,0.5,w);
			//}
		
			
			for(int j=1; j<=nGbins; j++){
				bincenter = binLow + gBinWidth;
				//does the ph oton start before and end after this bin center?
				if(sg_Rxy_origin < binLow ){
		//		if(sg_Rxy_origin < 1.0){
				if(sg_Rxy_endpt == -1 || sg_Rxy_endpt > bincenter){
					FillTH1(id_trueGeom, bincenter, w);	
				}}
				binLow += gBinWidth;
			}	
			binLow = gBinStart;
			bincenter=0;
// Float_t Rbins[] = { 1,5,9,13,18,20 };
  //    Int_t  binnum = 5;

			//coarse binning
			for(int j=0; j<binnum; j++){
				binLow = Rbins[j];
				bincenter = Rbins[j+1];
				if( sg_Rxy_origin < binLow ){
				if( sg_Rxy_endpt == -1 || sg_Rxy_endpt > bincenter ){
					FillTH1(id_trueGeom_Coarse, (binLow+bincenter)/2. , w);
				}}
			}
			binLow = gBinStart;
			bincenter=0;
			

		if( sg_Rxy_endpt != -1 ){
		if( exyz_ptype==14 ){
			FillTH1(id_trueConv, sg_Rxy_endpt, w);
			FillTH1(id_trueConv_Coarse, sg_Rxy_endpt, w);
		}}

		//check tracks and compare
		if( sg_Rxy_endpt != -1 ){
		if( exyz_ptype==14){
			pc_eidx = SG.pc_tke_idx[i];
			pc_pidx = SG.pc_tkp_idx[i];
			std::cout<<"pceidx pcpidx "<<pc_eidx<<" "<<pc_pidx<<std::endl;
			std::cout<<"pts "<<SimTrk_pt[pc_eidx] <<" "<<SimTrk_pt[pc_pidx]<<std::endl;
			if(SimTrk_pt[pc_eidx] < 0.2 || SimTrk_pt[pc_pidx] < 0.2){
				FillTH1(id_trueConv_gPass_trkFail, sg_Rxy_endpt,w);
				FillTH1(id_trueConv_gPass_trkFail_Coarse, sg_Rxy_endpt,w);
				if( SimTrk_pt[pc_eidx] < SimTrk_pt[pc_pidx]){
					FillTH1(id_minpt_gtrkpc, SimTrk_pt[pc_eidx], w);
				}
				else{
					FillTH1(id_minpt_gtrkpc, SimTrk_pt[pc_pidx], w);
				}

			}	
		}}
			
		}}}//}//end g acceptance

	}

//
/////////////////////////////////////////

/////////////////////////////////////////?SIM things
//new approach to simvtx thing

/*

	//get sim pc mask
	sim_pc SPC;
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
				FillTH1(id_trueConv, sqrt(sx*sx + sy*sy), w);		
			
				if(SimTrk_pdgId[t1idx] > 0){
					ptm = SimTrk_pt[t1idx];
					ptp = SimTrk_pt[t2idx];
				}
				else{	ptm
					ptm = SimTrk_pt[t2idx];
					ptp = SimTrk_pt[t1idx];
				}
	
				FillTH1(id_effXPD, ptp/(ptm+ptp), w);
			}}}}
		}
	}
*/
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
	double matchcut=0.5;
	for(int i=0; i<HGN.vsel.size(); i++){
		pcidx = HGN.vsel[i];
		//if we have a match
		if(Conv_vtxdl[pcidx] < matchcut){
		matchidx = Conv_convVtxIdx[pcidx];
		if(SimVtx_processType[matchidx] == 14){		
			FillTH1(id_effPtN, CVs[pcidx].pt, w);
			FillTH1(id_effRN, CVs[pcidx].radius, w);
			FillTH1(id_effXPN, CVs[pcidx].xplus, w);
			FillTH1(id_purityPtN, CVs[pcidx].pt, w);
			FillTH1(id_purityRN, CVs[pcidx].radius, w);
			FillTH1(id_purityXPN, CVs[pcidx].xplus, w);
			FillTH1(id_nconvR_match, CVs[pcidx].radius, w);
			FillTH1(id_nconvR_match_Coarse, CVs[pcidx].radius, w);
		//	FillTH1(id_nconvPt, CVs[pcidx].pt, w);
		//	FillTH1(id_nconvR, CVs[pcidx].radius, w);
		//	FillTH1(id_nconvXP, CVs[pcidx].xplus, w);
		}}
		FillTH1(id_purityPtD, CVs[pcidx].pt, w);
		FillTH1(id_purityRD, CVs[pcidx].radius, w);
		FillTH1(id_purityXPD, CVs[pcidx].xplus, w);
		FillTH1(id_nconvPt, CVs[pcidx].pt, w);
                FillTH1(id_nconvR, CVs[pcidx].radius, w);
                FillTH1(id_nconvXP, CVs[pcidx].xplus, w);
		FillTH1(id_nconvR_all, CVs[pcidx].radius, w);
		FillTH1(id_nconvR_all_Coarse, CVs[pcidx].radius, w);

		if(Conv_vtxdl[pcidx] > matchcut){
			FillTH1(id_nconvPt_fake, CVs[pcidx].pt, w);
			FillTH1(id_nconvR_fake, CVs[pcidx].radius, w);
			FillTH1(id_nconvXP_fake, CVs[pcidx].xplus, w);
		}
		if(Conv_vtxdl[pcidx] < matchcut){
		matchidx = Conv_convVtxIdx[pcidx];
		if(SimVtx_processType[matchidx] != 14){
			FillTH1(id_nconvPt_fake, CVs[pcidx].pt, w);
			FillTH1(id_nconvR_fake, CVs[pcidx].radius, w);
			FillTH1(id_nconvXP_fake, CVs[pcidx].xplus, w);
		}}	
		
	}

	/////end matching



	

}  // End of sub-program
#endif
