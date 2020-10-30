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
		}

	}
	FillTH1(id_numpccutHist, npcCut, w);
	
        //disambiguation (only do hungarian algorithm on pc's that pass cutmask)
	hgn_pc HGN = pc_disambiguation(s,cutmask);
	FillTH1( id_numHGNPCHist, HGN.vsel.size(), w);
	int cidx;
	for(int i=0; i<HGN.vsel.size(); i++){
		cidx = HGN.vsel[i];
		double PC_pt = sqrt(PC_Px[cidx]*PC_Px[cidx] + PC_Py[cidx]*PC_Py[cidx]);
		FillTH1( id_ptHCutHist, PC_pt, w);
		FillTH1( id_pzHCutHist, PC_Pz[cidx], w);
	
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
/*	
	std::vector<double> dL_c_v0(numberOfPC);
	std::vector<double> dL_b_v0(numberOfPC);
	std::vector<double> dL_x_v0(numberOfPC);
	std::vector<double> dL2_x_v0(numberOfPC);
	
        std::vector<double> dL_c_v0_sidx(numberOfPC);
        std::vector<double> dL_b_v0_sidx(numberOfPC);
        std::vector<double> dL_x_v0_sidx(numberOfPC);
        std::vector<double> dL2_x_v0_sidx(numberOfPC);
*/
// i will rerun this program with different cutmasks and name the output files different by hand


	
	//make all min big
	mindL_c = 999;
	mindL_b = 999;
	mindL_x = 999;
	mindL2_x = 999;
/*	for(int i=0; i< numberOfPC; i++){
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
		if( std::min(mindL_c, mindL_b) > 0.3) Type=2;
		recoType_v0[i] = Type;
		FillTH1(id_comp0, Type, w);
	}
*/
/*	std::vector<double> dL_c_v1(HGN.vsel.size());
        std::vector<double> dL_b_v1(HGN.vsel.size());
        std::vector<double> dL_x_v1(HGN.vsel.size());
        std::vector<double> dL2_x_v1(HGN.vsel.size());

        std::vector<double> dL_c_v1_sidx(HGN.vsel.size());
        std::vector<double> dL_b_v1_sidx(HGN.vsel.size());
        std::vector<double> dL_x_v1_sidx(HGN.vsel.size());
        std::vector<double> dL2_x_v1_sidx(HGN.vsel.size());
*/		
	std::vector<int> recoType_v1(HGN.vsel.size());//no selection on sim

	mindL_c = 999;
        mindL_b = 999;
        mindL_x = 999;
        mindL2_x = 999;
/*	for(int i=0; i< HGN.vsel.size(); i++){
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
         	if( std::min(mindL_c, mindL_b) > 0.3) Type=2;
                recoType_v1[i] = Type;
		FillTH1(id_comp1, Type, w);
	}// end HGN vsel loop
*/	

//need to set sim mask and reco mask every compile, they should be consistent
	double recoidx,recoxplus,recominTk,recogpt, rpt0, rpt1;
	std::vector<bool> reco_mask0(numberOfPC,true);
	std::vector<bool> reco_mask2(numberOfPC,false);
	std::vector<bool> reco_mask3(numberOfPC,false);
	std::vector<bool> _recomask;


//do a loop first and set recomask
/*	for(int i=0; i<HGN.vsel.size(); i++){
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
                        reco_mask2[recoidx] = true;
                }

	}
*/

//	_recomask=reco_mask0; ////////////////////////////////////////remember to change
	_recomask=reco_mask2;
//	_recomask=reco_mask3;
/*	for(int i=0; i< HGN.vsel.size(); i++){//HGN num eff loop
	    //if(recoType_v1[i] == 0){//matched to sim conversion
		recoidx = HGN.vsel[i];	
		recogpt = sqrt(PC_x[recoidx]*PC_x[recoidx] + PC_y[recoidx]*PC_y[recoidx]);
		rpt0 = sqrt(CVs[recoidx].px0p * CVs[recoidx].px0p + CVs[recoidx].py0p * CVs[recoidx].py0p);
	        rpt1 = sqrt(CVs[recoidx].px1p * CVs[recoidx].px1p + CVs[recoidx].py1p * CVs[recoidx].py1p);
		recominTk = std::min(rpt0,rpt1);
		recoxplus = CVs[recoidx].xplus;
		

	     if(recoType_v1[i] == 0 && _recomask[recoidx]){
		FillTH1(id_ptnum1, gpt, w);
		FillTH1(id_xpnum1, recoxplus, w);
		FillTH1(id_minTknum1, recominTk, w);
		FillTH1(id_Rnum1 , CVs[recoidx].radius, w);
		FillTH1(id_Rwidenum1, CVs[recoidx].radius, w);
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
*/
	//find all photons parent and endpoint
	//divide sectors
	//fluxing photon means parent is vtx before or in sector and endpoint farther than sector outer edge
	//converting photon means parent vtx is before or in sector and endpoint is in sector
	//look at pt dist of fluxing photons,
	
	//gmask, store endpoint, origin, processtype endpoint, PT
	std::vector<int> gparentIdx(nSimTrk,-1);//simvtx idx
	std::vector<int> gmask(nSimTrk,0);//mask of simtrks, pdg==22 and passes simcuts for photon (pt?) 0=not photon, 1=photon, 2=photon + and passes simcuts
	std::vector<int> gchildIdx(nSimTrk,-1);//simvtx idx
	std::vector<int> gtop14(nSimTrk,0);//0=no p14, 1=converts, 2= converts and conversion passes simcuts

	std::map<int,int> ptid_cvtx;
	//loop to map child tid and simvtx
	for(int i=0; i<nSimVtx; i++){
		//map every child and parent index together?
		ptid_cvtx.insert(pair<int,int>( SimVtx_simtrk_parent_tid[i], i));
		
	}

	//loop to establish parent
	for(int i=0; i<nSimTrk;i++){
		if(SimTrk_pdgId[i] == 22){
			gmask[i] = 1;
			if(SimTrk_pt[i] > 0.2){

				gmask[i] = 2;	
			}
			gparentIdx[i] = SimTrk_simvtx_Idx[i];

			//use the map to get child sim vtx index
			auto search =  ptid_cvtx.find( SimTrk_trackId[i] );
			gchildIdx[i] = search->second;	
			if( SimVtx_processType[ gchildIdx[i]  ] == 14 ){
				gtop14[i] = 1;
				//todo if check mask of conv to pass? for 2
			}
			//std::cout<<search->first<<" "<<search->second<<std::endl;	
		}	
		
	}	
	//loop again to establish child
	//int printonce=0;
	int SIMVTXCIDX;
	if(SPC.p14_key.size() >0 ){
	  SIMVTXCIDX = SPC.p14_key[0];
	}

	int startIdx;
	int endIdx;
	double start_x,start_y,start_z, start_r;
	double end_x,end_y,end_z, end_r;
	double GPT,THETAG, GPZ;
	//check to make sure mapping went okay
	//
	//define sectors
	//0<r<=5
	//5<r<=9
	//9<r<=14
	//14<r<=18
	//18<r<25
	for( int k=0; k<nSimTrk; k++){//gmask is size of sim trk
	  
            if(gmask[k] >0){
		 // std::cout<<"gparentvtxIdx "<<gparentIdx[k]<<" gmaskidx "<< k<< " gchildvtxIdx "<< gchildIdx[k]<<" mask "<<gmask[k]<<std::endl;
		startIdx = gparentIdx[k];
		endIdx = gchildIdx[k];
		start_x = SimVtx_x[startIdx];
		start_y = SimVtx_y[startIdx];
		start_z = SimVtx_z[startIdx];
		end_x = SimVtx_x[endIdx];
		end_y = SimVtx_y[endIdx];
		end_z = SimVtx_z[endIdx];
		//std::cout<<"g start position "<<start_x<<" "<<start_y<<" "<<start_z<<" "<<std::endl;
		//std::cout<<"g end position "<<end_x<<" "<<end_y<<" "<<end_z<<" converted? "<<gtop14[k]<<" ptype "<<SimVtx_processType[endIdx]<<std::endl;
		start_r = sqrt( start_x*start_x + start_y*start_y);
		end_r = sqrt(end_x*end_x + end_y*end_y);
		
			
		//sector 1 flux
		GPT = SimTrk_pt[k];
                GPZ = gpt* sinh( SimTrk_eta[k]);
               
                THETAG = atan2(GPT,GPZ);
		
		if( abs(start_z) < GV.ZCUT && abs(cos(THETAG)) < GV.COSTCUT){

		if( start_r < 5 ){
			FillTH1(id_gfluxpt_s1, GPT, w);
			FillTH1(id_ngflux_s1, 0, w);
			if( end_r < 5 && abs(end_z) <GV.ZCUT ){
				FillTH1(id_ngflux_pcall_s1, 0, w);	
			}
			if( end_r < 5 && sim_mask4[ gchildIdx[k] ] ){
				FillTH1(id_ngflux_pccut_s1,0,w);
			}
		}
		if( start_r < 9 && end_r > 5){
			FillTH1(id_gfluxpt_s2, GPT, w);
                        FillTH1(id_ngflux_s2, 0, w);
                        if( end_r < 9 && abs(end_z) <GV.ZCUT ){
                                FillTH1(id_ngflux_pcall_s2, 0, w);
                        }
                        if( end_r < 9 && sim_mask4[ gchildIdx[k] ] ){
                                FillTH1(id_ngflux_pccut_s2,0,w);
                        }

		}
		if( start_r < 14 && end_r > 9){
			FillTH1(id_gfluxpt_s3, GPT, w);
                        FillTH1(id_ngflux_s3, 0, w);
                        if( end_r < 14 && abs(end_z) <GV.ZCUT ){
                                FillTH1(id_ngflux_pcall_s3, 0, w);
                        }
                        if( end_r < 14 && sim_mask4[ gchildIdx[k] ] ){
                                FillTH1(id_ngflux_pccut_s3,0,w);
                        }

		}
		if( start_r < 18 && end_r > 14){
			FillTH1(id_gfluxpt_s4, GPT, w);
                        FillTH1(id_ngflux_s4, 0, w);
                        if( end_r < 18 && abs(end_z) <GV.ZCUT ){
                                FillTH1(id_ngflux_pcall_s4, 0, w);
                        }
                        if( end_r < 18 && sim_mask4[ gchildIdx[k] ] ){
                                FillTH1(id_ngflux_pccut_s4,0,w);
                        }

		}
		if( start_r < 25 && end_r > 18){
			FillTH1(id_gfluxpt_s5, GPT, w);
                        FillTH1(id_ngflux_s5, 0, w);
                        if( end_r < 25 && abs(end_z) <GV.ZCUT ){
                                FillTH1(id_ngflux_pcall_s5, 0, w);
                        }
                        if( end_r < 25 && sim_mask4[ gchildIdx[k] ] ){
                                FillTH1(id_ngflux_pccut_s5,0,w);
                        }

		}
		
		}//end z/cost conditional

	    }
	//    if(SPC.p14_key.size()>0 && gchildIdx[k] == SIMVTXCIDX){
	//	std::cout<<"found"<<std::endl;
	//	std::cout<<"idx of G from SPC "<<SPC.p14_g[SIMVTXCIDX]<< " SIMVTXCIDX "<< SIMVTXCIDX<<std::endl;
	  // }
	}

	

	
}  // End of sub-program
#endif
