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
#include <bitset>
#include <boost/dynamic_bitset.hpp>

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
		if(cutmask[i]){
			double PC_pt = sqrt(PC_Px[i]*PC_Px[i] + PC_Py[i] * PC_Py[i]);
			npcCut++;
			FillTH1(id_ptCutHist, PC_pt, w);
			FillTH1(id_pzCutHist, PC_Pz[i], w);
			FillTH2(id_xywideCPCHist, CVs[i].x, CVs[i].y, w);

		}

	}
	FillTH1(id_numpccutHist, npcCut, w);
	
        //disambiguation (only do hungarian algorithm on pc's that pass cutmask)
	hgn_pc HGN = pc_disambiguation(s,cutmask);
	std::vector<bool> HGNmask(numberOfPC,false);
	FillTH1( id_numHGNPCHist, HGN.vsel.size(), w);
	int cidx;
	for(int i=0; i<HGN.vsel.size(); i++){
		cidx = HGN.vsel[i];
		HGNmask[cidx] = true;	
		double PC_pt = sqrt(PC_Px[cidx]*PC_Px[cidx] + PC_Py[cidx]*PC_Py[cidx]);
		FillTH1( id_ptHCutHist, PC_pt, w);
		FillTH1( id_pzHCutHist, PC_Pz[cidx], w);
		FillTH2( id_xywideHGNPCHist, CVs[cidx].x, CVs[cidx].y, w); 	
	}

////////////////////////end example stuff


	//do some sim stuff
	//plot sim asym and min pt tracks
	int t0idx, t1idx, vidx;
	double xplus, ptasym;
	double pt0,pt1,q0, gpt,gpz, simz, simthetag;
	std::vector<double> xplus_v(SPC.p14_key.size());
	std::vector<double> minpt_v(SPC.p14_key.size());
	//define mask2 minptTk > 0.2
	//define mask3 0.1<xplus<0.9
	//sim mask 0 allow for all to pass
	std::vector<bool> sim_mask0(nSimVtx,true);
	//sim mask 2 minpt pair cut
        std::vector<bool> sim_mask2(nSimVtx,false);
        //sim mask 3 xplus cut
        std::vector<bool> sim_mask3(nSimVtx,false);
	//sim mask 4 minpt pair cut and z+costg cut (same as standard numerator cuts)
	std::vector<bool> sim_mask4(nSimVtx,false);

        std::vector<bool> _simmask;
        //uncomment the one being used
	
	
	//vidx is index of pc on simvtx
	//gidx, tXidx are index of SimTrack that correspond to simvtx_i
	for(int i=0; i<SPC.p14_key.size(); i++){
		vidx= SPC.p14_key[i];
		gidx= SPC.p14_g[vidx];
		t0idx= SPC.p14_t1[vidx];
		t1idx= SPC.p14_t2[vidx];
		
		pt0 = SimTrk_pt[t0idx];
		pt1 = SimTrk_pt[t1idx];
		q0 = -( SimTrk_pdgId[t0idx]/abs(SimTrk_pdgId[t0idx]) );		
		ptasym = (pt0-pt1)/(pt0+pt1);
		if (q0<0) ptasym = -ptasym;
        	xplus = (1.0 + ptasym)/2.0;
		xplus_v[i] = xplus;
		minpt_v[i] = std::min(pt0,pt1);
		FillTH1(id_xplusSPCHist, xplus, w);
		FillTH1(id_minTkPtSPCHist, std::min(pt0,pt1), w);
	
		gpt = SimTrk_pt[gidx];
		gpz = gpt* sinh( SimTrk_eta[gidx]);
		simz = SimVtx_z[vidx];
		simthetag = atan2(gpt,gpz);


		if( xplus < 0.9 && xplus > 0.1 ) sim_mask3[vidx] = true;
		if( std::min(pt0,pt1) > 0.2 ) sim_mask2[vidx] = true;
		if( std::min(pt0,pt1) > 0.2 && abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT ) sim_mask4[vidx] = true;		                           
//		if( std::min(pt0,pt1) > 0.2 && abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT && gpt>3) sim_mask4[vidx] = true;  /////REMember to change (for pt>5 run)	

		//if( _simmask[vidx] == true ){
	        //    FillTH1(id_ptden1, gpt, w);
	        //    FillTH1(id_xpden1, xplus, w);
                //    FillTH1(id_minTkden1, std::min(pt0,pt1), w)
	//	} 

	}	


	//uncomment the one being used
//      _simmask=sim_mask0;//////////////////////////////////////////////// remember to change
//        _simmask=sim_mask2;
//       _simmask=sim_mask3;
	_simmask=sim_mask4;
	double rsim;
	for(int i=0; i<nSimVtx; i++){
		if(_simmask[i]){
	            t0idx = SPC.p14_t1[i];
		    t1idx = SPC.p14_t2[i];
//		    gidx= SPC.p14_g[vidx];
		    gidx = SPC.p14_g[i];
		    gpt = SimTrk_pt[gidx];
 		    pt0 = SimTrk_pt[t0idx];
                    pt1 = SimTrk_pt[t1idx];
                    q0 = -( SimTrk_pdgId[t0idx]/abs(SimTrk_pdgId[t0idx]) );
                    ptasym = (pt0-pt1)/(pt0+pt1);
                    if (q0<0) ptasym = -ptasym;
                    xplus = (1.0 + ptasym)/2.0;
				     
		    rsim = sqrt(SimVtx_x[i]*SimVtx_x[i] + SimVtx_y[i]*SimVtx_y[i]);

                //    xplus_v[i] = xplus;
                //    minpt_v[i] = std::min(pt0,pt1);
		    FillTH1(id_ptden1, gpt, w);
                    FillTH1(id_xpden1, xplus, w);
                    FillTH1(id_minTkden1, std::min(pt0,pt1), w);
		    FillTH1(id_Rden1, rsim, w);
		    FillTH1(id_Rwideden1,rsim, w);
		   if( rsim <= 5){
			FillTH1(id_rden_s1, 0, w);
		   }
		   if( rsim > 5 && rsim <= 9 ){
			FillTH1(id_rden_s2, 0, w);
		   }
		   if( rsim > 9 && rsim <= 14 ){
                        FillTH1(id_rden_s3, 0, w);
                   }
		   if( rsim > 14 && rsim <= 18 ){
                        FillTH1(id_rden_s4, 0, w);
                   }
		   if( rsim > 18 && rsim <= 25 ){
                        FillTH1(id_rden_s5, 0, w);
                   }
		
 

		}
	}
/////////////////////
	//study some cuts at sim level pc
	



////////////////////////
	//do some matching with reco do it both for HGNpc and PC with no cuts
	//2d plot distance between pc and bg svtx and asymmetry
	double sx,sy,sz;
	double cx,cy,cz;
	double dL; //distance between two points in 3d
	double mindL_c, mindL_b; //the nearest p14 and nearest vtx != p14 respectively
	double mindL_x, mindL2_x; // the first and second nearest vtx (any process type)
	double idx_c, idx_b, idx_x, idx_x2;

	//use dL cut of 0.3 on dL_x to designate matching
	std::vector<int> recoType_v0(numberOfPC);//define 3 categories 0=signal 1=bkg 2=unmatched bkg
	int Type;
	
	std::vector<double> dL_c_v0(numberOfPC);
	std::vector<double> dL_b_v0(numberOfPC);
	std::vector<double> dL_x_v0(numberOfPC);
	std::vector<double> dL2_x_v0(numberOfPC);
	
        std::vector<double> dL_c_v0_sidx(numberOfPC);
        std::vector<double> dL_b_v0_sidx(numberOfPC);
        std::vector<double> dL_x_v0_sidx(numberOfPC);
        std::vector<double> dL2_x_v0_sidx(numberOfPC);
// i will rerun this program with different cutmasks and name the output files different by hand


	
	//make all min big
	mindL_c = 999;
	mindL_b = 999;
	mindL_x = 999;
	mindL2_x = 999;
	for(int i=0; i< numberOfPC; i++){
		cx = PC_x[i];
		cy = PC_y[i];
		cz = PC_z[i];
		for(int j=0; j<nSimVtx; j++){
			sx = SimVtx_x[j];
			sy = SimVtx_y[j];
			sz = SimVtx_z[j];
			dL = sqrt( (sx-cx)*(sx-cx) + (sy-cy)*(sy-cy) + (sz-cz)*(sz-cz) );
			if( SimVtx_processType[j] == 14  && _simmask[j]==true ){
				if( dL < mindL_c ){
				    mindL_c = dL;
				    idx_c = j; 	 
				}
			}
			if( SimVtx_processType[j] != 14){
				if( dL < mindL_b ){
				    mindL_b = dL;
				    idx_b = j;
				}
			}
			if( dL < mindL_x ){
			    mindL_x = dL;
			    idx_x = j;
			}
			if( dL > mindL_x && dL < mindL2_x ){
			    mindL2_x = dL;	
			    idx_x2 = j;
			}
				
		}
		dL_c_v0[i] = mindL_c;
		dL_b_v0[i] = mindL_b;
		dL_x_v0[i] = mindL_x;
		dL2_x_v0[i] = mindL2_x;

		dL_c_v0_sidx[i] = idx_c;
                dL_b_v0_sidx[i] = idx_b;
                dL_x_v0_sidx[i] = idx_x;
                dL2_x_v0_sidx[i] = idx_x2;

		
		FillTH1(id_mindLc0, mindL_c, w);
		FillTH1(id_mindLb0, mindL_b, w);
		FillTH1(id_mindLx0, mindL_x, w);
		FillTH1(id_mindL2x0, mindL2_x , w);
		
		if( mindL_c < mindL_b) Type =0;
		if( mindL_b < mindL_c) Type =1;
	//	if( mindL_x > 0.3) Type =2;
		if( std::min(mindL_c, mindL_b) > 0.5) Type=2;
		recoType_v0[i] = Type;
		FillTH1(id_comp0, Type, w);
	}

	std::vector<double> dL_c_v1(HGN.vsel.size());
        std::vector<double> dL_b_v1(HGN.vsel.size());
        std::vector<double> dL_x_v1(HGN.vsel.size());
        std::vector<double> dL2_x_v1(HGN.vsel.size());

        std::vector<double> dL_c_v1_sidx(HGN.vsel.size());
        std::vector<double> dL_b_v1_sidx(HGN.vsel.size());
        std::vector<double> dL_x_v1_sidx(HGN.vsel.size());
        std::vector<double> dL2_x_v1_sidx(HGN.vsel.size());
		
	std::vector<int> recoType_v1(HGN.vsel.size());//no selection on sim

	mindL_c = 999;
        mindL_b = 999;
        mindL_x = 999;
        mindL2_x = 999;
	for(int i=0; i< HGN.vsel.size(); i++){
		cx = PC_x[i];
                cy = PC_y[i];
                cz = PC_z[i];
                for(int j=0; j<nSimVtx; j++){
                        sx = SimVtx_x[j];
                        sy = SimVtx_y[j];
                        sz = SimVtx_z[j];
                        dL = sqrt( (sx-cx)*(sx-cx) + (sy-cy)*(sy-cy) + (sz-cz)*(sz-cz) );
                        if( SimVtx_processType[j] == 14 && _simmask[j] == true ){
                                if( dL < mindL_c ){ 
                                    mindL_c = dL;
                                    idx_c = j;
                                }
                        }
                        if( SimVtx_processType[j] != 14){
                                if( dL < mindL_b ){ 
                                    mindL_b = dL;
                                    idx_b = j;
                                }
                        }
                        if( dL < mindL_x ){
                            mindL_x = dL;
                            idx_x = j;
                        }
                        if( dL > mindL_x && dL < mindL2_x ){ 
                            mindL2_x = dL; 
                            idx_x2 = j;
                        } 
                }//end nSimVtx loop
		dL_c_v0[i] = mindL_c;
                dL_b_v0[i] = mindL_b;
                dL_x_v0[i] = mindL_x;
                dL2_x_v0[i] = mindL2_x;

                dL_c_v0_sidx[i] = idx_c;
                dL_b_v0_sidx[i] = idx_b;
                dL_x_v0_sidx[i] = idx_x;
                dL2_x_v0_sidx[i] = idx_x2;

		
		FillTH1(id_mindLc1, mindL_c, w);
                FillTH1(id_mindLb1, mindL_b, w);
                FillTH1(id_mindLx1, mindL_x, w);
                FillTH1(id_mindL2x1, mindL2_x , w);

		if( mindL_c < mindL_b) Type =0;
                if( mindL_b < mindL_c) Type =1;
         //       if( mindL_x > 0.3) Type =2;
         	if( std::min(mindL_c, mindL_b) > 0.5) Type=2;
                recoType_v1[i] = Type;
		FillTH1(id_comp1, Type, w);
	}// end HGN vsel loop
	

//need to set sim mask and reco mask every compile, they should be consistent
	double recoidx,recoxplus,recominTk,recogpt, rpt0, rpt1;
	std::vector<bool> reco_mask0(numberOfPC,true);
	std::vector<bool> reco_mask2(numberOfPC,false);
	std::vector<bool> reco_mask3(numberOfPC,false);
	std::vector<bool> _recomask;


//do a loop first and set recomask
	for(int i=0; i<HGN.vsel.size(); i++){
		recoidx = HGN.vsel[i];
                recogpt = sqrt(PC_x[recoidx]*PC_x[recoidx] + PC_y[recoidx]*PC_y[recoidx]);
                rpt0 = sqrt(CVs[recoidx].px0p * CVs[recoidx].px0p + CVs[recoidx].py0p * CVs[recoidx].py0p);
                rpt1 = sqrt(CVs[recoidx].px1p * CVs[recoidx].px1p + CVs[recoidx].py1p * CVs[recoidx].py1p);
                recominTk = std::min(rpt0,rpt1);
                recoxplus = CVs[recoidx].xplus;
                if( recoxplus<0.9 && recoxplus>0.1){
                        reco_mask3[recoidx] = true;
                }
                if( recominTk > 0.2 ){
                //if( recominTk > 0.2 && recogpt>3){///////////////////////////remember to change (for pt>5 radial eff run)
                        reco_mask2[recoidx] = true;
                }

	}



//	_recomask=reco_mask0; ////////////////////////////////////////remember to change
	_recomask=reco_mask2;
//	_recomask=reco_mask3;
//	std::vector<int> TYPEmask(numberOfPC,-1);
	for(int i=0; i< HGN.vsel.size(); i++){//HGN num eff loop
	    //if(recoType_v1[i] == 0){//matched to sim conversion
		recoidx = HGN.vsel[i];	
		recogpt = sqrt(PC_x[recoidx]*PC_x[recoidx] + PC_y[recoidx]*PC_y[recoidx]);
		rpt0 = sqrt(CVs[recoidx].px0p * CVs[recoidx].px0p + CVs[recoidx].py0p * CVs[recoidx].py0p);
	        rpt1 = sqrt(CVs[recoidx].px1p * CVs[recoidx].px1p + CVs[recoidx].py1p * CVs[recoidx].py1p);
		recominTk = std::min(rpt0,rpt1);
		recoxplus = CVs[recoidx].xplus;
		

	     if(recoType_v1[i] == 0 && _recomask[recoidx]){
		FillTH1(id_ptnum1, recogpt, w);
		FillTH1(id_xpnum1, recoxplus, w);
		FillTH1(id_minTknum1, recominTk, w);
		FillTH1(id_Rnum1 , CVs[recoidx].radius, w);
		FillTH1(id_Rwidenum1, CVs[recoidx].radius, w);
		
		   if( CVs[recoidx].radius <= 5){
                        FillTH1(id_rnum_s1, 0, w);
                   }
                   if( CVs[recoidx].radius > 5 && CVs[recoidx].radius <= 9 ){
                        FillTH1(id_rnum_s2, 0, w);
                   }
                   if( CVs[recoidx].radius > 9 && CVs[recoidx].radius <= 14 ){
                        FillTH1(id_rnum_s3, 0, w);
                   }
                   if( CVs[recoidx].radius > 14 && CVs[recoidx].radius <= 18 ){
                        FillTH1(id_rnum_s4, 0, w);
                   }
                   if( CVs[recoidx].radius > 18 && CVs[recoidx].radius <= 25 ){
                        FillTH1(id_rnum_s5, 0, w);
                   }
		

	     }
	     if(recoType_v1[i] == 1){
		FillTH1(id_ptnumB1, gpt, w);
		FillTH1(id_xpnumB1, recoxplus, w);
		FillTH1(id_minTknumB1, recominTk, w);
		FillTH1(id_RnumB1 , CVs[recoidx].radius, w);
                FillTH1(id_RwidenumB1, CVs[recoidx].radius, w);
	     }
	     if(recoType_v1[i] == 2){
		FillTH1(id_ptnumU1, gpt, w);
	        FillTH1(id_xpnumU1, recoxplus, w);
	        FillTH1(id_minTknumU1, recominTk, w);
		FillTH1(id_RnumU1 , CVs[recoidx].radius, w);
                FillTH1(id_RwidenumU1, CVs[recoidx].radius, w);
             }

	     

		
	}//end HGN vsel loop for eff num
////////////////////////CUTFLOW//////////////////////////////////////	
        enum cuts{ _base= 0,_r=1, _z=2, _cost=3, _prob=4, _nb40=5, _nb41=6, _HGN=7, _minpt0=8, _minpt1=9, _numcuts=10};
	boost::dynamic_bitset<uint16_t> bitS(_numcuts);	
        boost::dynamic_bitset<uint16_t> bitB(_numcuts);
	boost::dynamic_bitset<uint16_t> bitU(_numcuts);
	for(int i=0; i<numberOfPC; i++){
		if(recoType_v0[i]==0) bitS.set(_base);
                if(recoType_v0[i]==1) bitB.set(_base);
                if(recoType_v0[i]==2) bitU.set(_base);

		if( CVs[i].rerr < GV.RERRCUT){
			if(recoType_v0[i]==0) bitS.set(_r);
			if(recoType_v0[i]==1) bitB.set(_r);
			if(recoType_v0[i]==2) bitU.set(_r);
		}
		if( abs(CVs[i].z) < GV.ZCUT){
			if(recoType_v0[i]==0) bitS.set(_z);
			if(recoType_v0[i]==1) bitB.set(_z);
			if(recoType_v0[i]==2) bitU.set(_z);
		}
		if( abs(cos( CVs[i].theta)) < GV.COSTCUT){
			if(recoType_v0[i]==0) bitS.set(_cost);
			if(recoType_v0[i]==1) bitB.set(_cost);
			if(recoType_v0[i]==2) bitU.set(_cost);	
		}
		if( CVs[i].pfit > GV.FITPROBCUT ){
			if(recoType_v0[i]==0) bitS.set(_prob);
			if(recoType_v0[i]==1) bitB.set(_prob);
			if(recoType_v0[i]==2) bitU.set(_prob);
		}
		if( PC_vTrack0_nBefore[i] == 0){
			if(recoType_v0[i]==0) bitS.set(_nb40);
			if(recoType_v0[i]==1) bitB.set(_nb40);
			if(recoType_v0[i]==2) bitU.set(_nb40);
		}
		if( PC_vTrack1_nBefore[i] == 0){
			if(recoType_v0[i]==0) bitS.set(_nb41);
			if(recoType_v0[i]==1) bitB.set(_nb41);
			if(recoType_v0[i]==2) bitU.set(_nb41);
		}
		if( HGNmask[i] ){//NOTE this cut has all of the preceding cuts built into it!!!!!
			if(recoType_v0[i]==0) bitS.set(_HGN);
			if(recoType_v0[i]==1) bitB.set(_HGN);
			if(recoType_v0[i]==2) bitU.set(_HGN);
		}
		rpt0 = sqrt(CVs[i].px0p * CVs[i].px0p + CVs[i].py0p * CVs[i].py0p);
                rpt1 = sqrt(CVs[i].px1p * CVs[i].px1p + CVs[i].py1p * CVs[i].py1p);
	
		if( rpt0 > 0.2 ){
			if(recoType_v0[i]==0) bitS.set(_minpt0);
			if(recoType_v0[i]==1) bitB.set(_minpt0);
			if(recoType_v0[i]==2) bitU.set(_minpt0);
		} 
		if( rpt1 > 0.2 ){
			if(recoType_v0[i]==0) bitS.set(_minpt1);
			if(recoType_v0[i]==1) bitB.set(_minpt1);
			if(recoType_v0[i]==2) bitU.set(_minpt1);
		}
	}

	//TYPE 0
	bool flag = true;
        for(int i=0; i<bitS.size(); i++){
                if( !bitS.test(i) ) flag = false;
                if( flag) FillTH1(id_pcScutflow,i,w);
        }
	//TYPE 1
	flag = true;
	for(int i=0; i<bitB.size(); i++){ 
                if( !bitB.test(i) ) flag = false;
                if( flag) FillTH1(id_pcBcutflow,i,w);
        }	
	//TYPE 2
	flag = true;
	for(int i=0; i<bitU.size(); i++){ 
                if( !bitU.test(i) ) flag = false;
                if( flag) FillTH1(id_pcUcutflow,i,w);
        }


	
}  // End of sub-program
#endif
