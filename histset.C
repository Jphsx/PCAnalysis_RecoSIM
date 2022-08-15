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
#include <map>
//#include <bitset>
//#include <boost/dynamic_bitset.hpp>

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
	//double w_ndof = w_pvndof0_BC;
	double w_sumtk = w_pvtk_BC;
	//double w_evt = 54.7157747; //with ndof 100 cut
	//double w_evt = 28.3953880;
//	double w_evt = 12.1005464; //with sumtrk 40 cut	
//	double w_evt = 17.079357593;//fixed from 12, normalized to weight dist integrals corresponds to the 40 cut
//	double w_evt = 17.0789640;// corresponds to 20 cut
//	double w_evt = 16.83468670;// corresponds to 8-20 cut
	double w_evt = 20.23366; //bug fix reweighting
	//if(nPV == 1 && PV_ndof[0] >50) w_ndof=1.0;
	//w = w_ndof * w_evt;
	w = w_sumtk* w_evt;
	
	FillTH1(id_npv, nPV,w); 
        //if(nPV != 1 || PV_ndof[0] > 50) return;
       
//	if(nPV==0) return;
//	if(nPV>0 && abs( PV_Z[0] ) > 1) return;

 
	if(nPV == 1)
        FillTH1(id_PVndof, PV_ndof[0],w);			
	
	double sumtkw=0;
	double sum_sumtkw=0;
	//loop over vertices get sumwi
	for(int v=0; v<nPV; v++){
		sumtkw = ((PV_ndof[v]+3.)/2.);
		sum_sumtkw += sumtkw;
		FillTH1(id_pvz, PV_Z[v], w);
		FillTH1(id_pvtrksum, sumtkw, w);
		if(nPV>1 && v>0)
		FillTH1(id_pvdz, PV_Z[0] - PV_Z[v], w);
		
	}	
	//FillTH1(id_sumpvtrksum, sum_sumtkw, w);

/////Built in sum_sumtkw <= 40 built into everything!!
	//if( sum_sumtkw > 40 ) return;
	//if( sum_sumtkw > 20 ) return; //testing the tight cut
	if( sum_sumtkw > 20 || sum_sumtkw < 8 ) return; //tight cut with low cut
//	if( sum_sumtkw > 40 || sum_sumtkw < 10 ) return;

	FillTH1(id_sumpvtrksum, sum_sumtkw, w);	

        //plot "raw" conv stuff
        FillTH1( id_numpcHist, numberOfPC,w);
std::vector<CommonVars> CVs = GetCommonVars(s,false);

        for(int i=0; i<numberOfPC; i++){
                double PC_pt = sqrt(PC_Px[i]*PC_Px[i] + PC_Py[i] * PC_Py[i]);
                FillTH1( id_pzHist, PC_Pz[i],w);
                FillTH1( id_ptHist, PC_pt,w);
                FillTH2( id_xyHist, PC_x[i], PC_y[i], w);
                FillTH2( id_xywideHist, PC_x[i],PC_y[i], w);
                double PC_r = sqrt(PC_x[i]*PC_x[i] + PC_y[i]*PC_y[i]);
                double PC_phi = atan2(PC_y[i],PC_x[i]);
                FillTH2(id_rphiHist,PC_r, PC_phi,w);
		FillTH1(id_r25RHist , CVs[i].radius,w);
		FillTH1(id_etaHist, CVs[i].etaphys, w);
                FillTH2(id_ndof_pcReta, PV_ndof[0] ,CVs[i].etaphys,w);
                FillTH2(id_ndof_pcRpt, PV_ndof[0], PC_pt,w);

	
        }

        //calculate common variables
//        std::vector<CommonVars> CVs = GetCommonVars(s,false);

        //get cut mask (check numerator cuts)
        std::vector<bool> cutmask = GetCutMask(s,CVs);
        int npcCut=0;
        for(int i=0; i<cutmask.size(); i++){
                if(cutmask[i]){
                        double PC_pt = sqrt(PC_Px[i]*PC_Px[i] + PC_Py[i] * PC_Py[i]);
                        npcCut++;
                        FillTH1(id_ptCutHist, PC_pt, w);
                        FillTH1(id_pzCutHist, PC_Pz[i], w);
                       // FillTH2(id_xywideCPCHist, CVs[i].x, CVs[i].y, w);
			FillTH1(id_r25CHist , CVs[i].radius,w);
			FillTH1(id_etaCutHist, CVs[i].etaphys, w);
			
                }

        }
        FillTH1(id_numpccutHist, npcCut, w);

        //fill "hgn" conv stuff

        hgn_pc HGN = pc_disambiguation(s,cutmask);
        std::vector<bool> HGNmask(numberOfPC,false);
	double cutdL = 0.5;
	//std::map<int, std::pair<int, double> > pcMatchColl = getPCMatchingColl(s,cutdL); // get matches for eff numerator
        std::map<int, std::vector<double> > pcMatchColl = getPCMatchingColl(s,cutdL);
	FillTH1( id_numHGNPCHist, HGN.vsel.size(), w);
        int cidx;
	double ptp;
	//std::pair<int,double> match_criteria;
	std::vector<double>  match_criteria;

	int leadfound,subfound,leadlost,sublost,leadqual,subqual;
        for(int i=0; i<HGN.vsel.size(); i++){
                cidx = HGN.vsel[i];
                HGNmask[cidx] = true;
                double PC_pt = sqrt(PC_Px[cidx]*PC_Px[cidx] + PC_Py[cidx]*PC_Py[cidx]);
                FillTH1( id_ptHCutHist, PC_pt, w);
                FillTH1( id_pzHCutHist, PC_Pz[cidx], w);
         //       FillTH2( id_xywideHGNPCHist, CVs[cidx].x, CVs[cidx].y, w);
        	FillTH1(id_r25HHist , CVs[cidx].radius,w);
		FillTH1(id_r25Hist_b2p5, CVs[cidx].radius, w);
		FillTH1(id_r25Hist_b2p5_nowt, CVs[cidx].radius, 1);
		FillTH1(id_r25coarse, CVs[cidx].radius, w);
		FillTH1(id_etaHGNHist, CVs[cidx].etaphys,w);
                FillTH2(id_ndof_pcHeta, PV_ndof[0] ,CVs[cidx].etaphys,w);
                FillTH2(id_ndof_pcHpt, PV_ndof[0], PC_pt,w);

		FillTH1(id_pc_chi2ndof,Conv_vtx_normalizedChi2[cidx], w); 
		FillTH2(id_nchi2_r,CVs[cidx].radius,Conv_vtx_normalizedChi2[cidx], w);
 		FillTH2(id_nchi2_dr,CVs[cidx].rerr,Conv_vtx_normalizedChi2[cidx],w);

		
		FillTH1(id_rerrHGNHist,CVs[cidx].rerr ,w);
		FillTH2(id_pt_rerr, CVs[cidx].rerr, PC_pt,w);
		FillTH2(id_r_rerr, CVs[cidx].rerr, CVs[cidx].radius,w);

		FillTH2(id_r_phi_p6, CVs[cidx].radius, CVs[cidx].phi , w);
	
		FillTH2(id_reta_pc, CVs[cidx].etaphys, CVs[cidx].radius, w);

		if(CVs[cidx].radius < 8)
		FillTH1( id_sumtksum_rlo, sum_sumtkw, w);

		if(CVs[cidx].radius > 8 && CVs[cidx].radius < 25 )
		FillTH1( id_sumtksum_rhi, sum_sumtkw, w);

		match_criteria = pcMatchColl[cidx];
		FillTH1(id_matchdR, match_criteria.at(1), w);

		
		if(match_criteria.at(0) ==0 && match_criteria.at(1) < cutdL){
			//FillTH1(id_eRnum, CVs[cidx].radius, w);
			//FillTH1(id_eRnum_b2p5, CVs[cidx].radius, w);
			FillTH1(id_eRnum, match_criteria.at(2), w);
			FillTH1(id_eRnum_nowt, match_criteria.at(2), 1);
			FillTH1(id_eRnum_b2p5, match_criteria.at(2), w);
			FillTH1(id_eRnum_b2p5_nowt, match_criteria.at(2), 1);
                	FillTH1(id_ePtnum,CVs[cidx].pt,w);
                	FillTH1(id_etksumnum, sum_sumtkw, w);
			FillTH2(id_effptr_num,match_criteria.at(2), CVs[cidx].pt, w);
			FillTH2(id_effptr_num_nowt, match_criteria.at(2), CVs[cidx].pt, 1);
			FillTH2(id_effrphiN, match_criteria.at(2), CVs[cidx].phi, w);	
			FillTH2(id_effrphiN_nowt, match_criteria.at(2), CVs[cidx].phi, 1);		
			
			FillTH2(id_reta_effN, GetSimGfromTID(s,match_criteria.at(3) ).eta_physics , match_criteria.at(2), w);
			FillTH2(id_reta_effN_nowt, GetSimGfromTID(s,match_criteria.at(3) ).eta_physics , match_criteria.at(2), 1);

			//PC_vTrack0_pt	
			if(PC_vTrack0_pt[cidx] >= PC_vTrack1_pt[cidx]){
				leadfound= PC_vTrack0_found[cidx];
				leadlost=  PC_vTrack0_lost[cidx];
				leadqual=  PC_vTrack0_quality[cidx];
				subfound=  PC_vTrack1_found[cidx];
				sublost=   PC_vTrack1_lost[cidx];
				subqual=   PC_vTrack1_quality[cidx];
			}else{
				leadfound= PC_vTrack1_found[cidx];
                                leadlost=  PC_vTrack1_lost[cidx];
                                leadqual=  PC_vTrack1_quality[cidx];
                                subfound=  PC_vTrack0_found[cidx];
                                sublost=   PC_vTrack0_lost[cidx];
                                subqual=   PC_vTrack0_quality[cidx];
			}	

			//track analysis by pt matched
			FillTH2(id_pt_leadfound, leadfound, CVs[cidx].pt, w); 
			FillTH2(id_pt_subfound, subfound, CVs[cidx].pt, w);
			FillTH2(id_pt_leadlost, leadlost, CVs[cidx].pt, w);
			FillTH2(id_pt_sublost, sublost, CVs[cidx].pt, w);
			FillTH2(id_pt_leadqual, leadqual, CVs[cidx].pt, w);
			FillTH2(id_pt_subqual, subqual, CVs[cidx].pt, w);
			FillTH2(id_pt_shared, PC_nSharedHits[cidx], CVs[cidx].pt, w); 
			//
		}
		if(match_criteria.at(0) !=0 || match_criteria.at(1) >= cutdL){
			//fill bg eff radial stuff
			FillTH1(id_eRnumf, CVs[cidx].radius, w);
                        FillTH1(id_eRnumf_b2p5, CVs[cidx].radius, w);
                        FillTH1(id_ePtnumf,CVs[cidx].pt,w);
                        FillTH1(id_etksumnumf, sum_sumtkw, w);

		}
	}


/////////////////// Efficiency stuff
	sim_pc SPC = GetSimPC(s);//need to apply cuts by hand
//	std::vector< std::pair<int, int> > GColl = getGParentColl(s);
//	std::vector< std::pair<int, double> > pcMatchColl = getPCMatchingColl(s, 0.5);
	//debug GColl
/*	std::cout<<" EVENT \n";
	std::cout<<" NUMPC =" << numberOfPC <<"\n";
	std::pair<std::vector<double>, std::vector<double> > points;
	for(int i=0; i<GColl.size(); i++){
		std::cout<<"G"<<i<<" idx: "<<GColl[i].first<<" Type: "<<GColl[i].second<<" ";
		points = getGEndpoints( s, GColl[i].first);
		std::cout<<"(sx,sy,sz)= "<<points.first[0]<<" "<<points.first[1]<<" "<<points.first[2]<<" ";
		std::cout<<"(ex,ey,ez)= "<<points.second[0]<<" "<<points.second[1]<<" "<<points.second[2]<<"\n";
	}
	for(int i=0; i<pcMatchColl.size(); i++){
		std::cout<<"PC"<<i<<" Mtype: "<<pcMatchColl[i].first<<" dR: "<<pcMatchColl[i].second<<"\n";
	}
*/

//////////////////////testing masking from early build

        std::vector<bool> sim_mask4(nSimVtx,false);

//        std::vector<bool> _simmask;
        //uncomment the one being used

	 double gpt,gpz,geta, tpt1,tpt2, costg, zpos;
	 int gidx,t1idx,t2idx;
	 double simr, sx,sy,simphi;
	int vidx,t0idx;
	double pt0,pt1,q0,simz, simthetag;
        //vidx is index of pc on simvtx
        //gidx, tXidx are index of SimTrack that correspond to simvtx_i
        for(int i=0; i<SPC.p14_key.size(); i++){
                vidx= SPC.p14_key[i];
                gidx= SPC.p14_g[vidx];
                t0idx= SPC.p14_t1[vidx];
                t1idx= SPC.p14_t2[vidx];

                pt0 = SimTrk_pt[t0idx];
                pt1 = SimTrk_pt[t1idx];
                //q0 = -( SimTrk_pdgId[t0idx]/abs(SimTrk_pdgId[t0idx]) );
                //ptasym = (pt0-pt1)/(pt0+pt1);
                //if (q0<0) ptasym = -ptasym;
                //xplus = (1.0 + ptasym)/2.0;
                //xplus_v[i] = xplus;
                //minpt_v[i] = std::min(pt0,pt1);
                //FillTH1(id_xplusSPCHist, xplus, w);
                //FillTH1(id_minTkPtSPCHist, std::min(pt0,pt1), w);

                gpt = SimTrk_pt[gidx];
                gpz = gpt* sinh( SimTrk_eta[gidx]);
                simz = SimVtx_z[vidx];
                simthetag = atan2(gpt,gpz);


         //       if( xplus < 0.9 && xplus > 0.1 ) sim_mask3[vidx] = true;
           //     if( std::min(pt0,pt1) > 0.2 ) sim_mask2[vidx] = true;
                //if( std::min(pt0,pt1) > 0.2 && abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT ) sim_mask4[vidx] = true;
//              if( std::min(pt0,pt1) > 0.2 && abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT && gpt>3) sim_mask4[vidx] = true;  /////REMember to change (for pt>5 run)   
//		if( abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT && gpt>0.4  && abs(SimTrk_eta[gidx]) < 0.1  ) sim_mask4[vidx] = true;// changing this for graham
		if( abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT && gpt>0.4 && abs(SimTrk_eta[gidx]) < 0.1 ) sim_mask4[vidx] = true;
	//	if( abs(simz) <GV.ZCUT && abs(cos(simthetag))<GV.COSTCUT && gpt>1.0 && abs(SimTrk_eta[gidx]) < 0.1 ) sim_mask4[vidx] = true;          
      //if( _simmask[vidx] == true ){
                //    FillTH1(id_ptden1, gpt, w);
                //    FillTH1(id_xpden1, xplus, w);
                //    FillTH1(id_minTkden1, std::min(pt0,pt1), w)
        //      } 

        }




///////////////////
	
	//loop over spc mask, put stuff that passes cuts into denom
	//for this restrict to r<25 for pt and 
//	double gpt,gpz,geta, tpt1,tpt2, costg, zpos;
//	int gidx,t1idx,t2idx;
//	double simr, sx,sy;
	double simE;
	for(int i=0; i<nSimVtx; i++){
		//debug section
		if(SimVtx_processType[i] == 14){
			sx = SimVtx_x[i];
	                sy = SimVtx_y[i];
        	        simr = sqrt(sx*sx + sy*sy);

			FillTH1(id_debug1, simr, 1);
			FillTH1(id_debug2, simr, 1);
		}

		//if( SPC.sim_mask[i] != 14 ) continue;
		if(sim_mask4[i]){	

		gidx = SPC.p14_g[i];
		t1idx = SPC.p14_t1[i];
		t2idx = SPC.p14_t2[i];
		
		gpt = SimTrk_pt[gidx];
		tpt1 = SimTrk_pt[t1idx];
		tpt2 = SimTrk_pt[t2idx];

		geta = SimTrk_eta[gidx];
		gpz = gpt*sinh(geta);
		costg = cos( atan2(gpt,gpz) );
		simE = gpt*cosh(geta);

		zpos = SimVtx_z[i];
		sx = SimVtx_x[i];
		sy = SimVtx_y[i];
		simr = sqrt(sx*sx + sy*sy);
		

		simphi=	atan2(sy,sx);
		if (simphi < 0) { simphi += 2 * M_PI; }
		//does it pass cuts? Require R<25 -- currently looking at context central bpix only
	//	if( abs(zpos) < GV.ZCUT && abs(costg) < GV.COSTCUT && tpt1 > GV.MINPT && tpt2 > GV.MINPT && simr<25  ){
			FillTH1(id_eRden, simr, w);
			FillTH1(id_eRden_b2p5, simr, w);
			FillTH1(id_eRden_nowt, simr, 1);
			FillTH1(id_eRden_b2p5_nowt, simr, 1);
			FillTH1(id_ePtden,gpt,w);
			FillTH1(id_etksumden, sum_sumtkw, w);
			FillTH2(id_effptr_den,simr, gpt, w);
			FillTH2(id_effptr_den_nowt,simr,gpt,1);
			FillTH2(id_effER_den, simr, simE, w);

			
			FillTH2(id_effrphiD, simr,simphi,w); 
			FillTH2(id_effrphiD_nowt, simr,simphi,1);			

			FillTH2(id_reta_effD, geta, simr, w);
			FillTH2(id_reta_effD_nowt, geta, simr, 1);
//		}
		}//end mask4 check	
	}
///fresh flux plot stuff
//
	std::vector<int> gidxlist =  getSimpleGidxList(s);
	std::pair<std::vector<double>, std::vector<double> > g_endpoints{};
	double Sx,Sy,Sz,Ex,Ey,Ez;
	double Rorigin,Rfinal,Rdrift, xdrift, ydrift, tdrift;
	double g_eta_physics;
	//double gpt_forcut;
	double gE;
	//create upper bin edge vector
	std::vector<double> rbinedge{};
	std::vector<double> rbincenter{};
	int nbinsx;
	auto gfluxhist = TH1Manager.at(id_gflux)->Get();
        

	nbinsx = gfluxhist->GetNbinsX();
	for(int i=1; i<nbinsx; i++){
		rbinedge.push_back( gfluxhist->GetBinCenter(i) + (gfluxhist->GetBinWidth(i)/2.) );
		rbincenter.push_back(gfluxhist->GetBinCenter(i));
	}

	//debug what are bin edges?
	//std::cout<<"edge\n";
	//for(int i=0; i<rbinedge.size(); i++){
	//	std::cout<<rbinedge.at(i)<<" ";
	//}
	//std::cout<<std::endl;
	
	//loop over valid photons
	for(int i=0; i< gidxlist.size(); i++){
		//collect the endpoints for every valid photon
		g_endpoints = getGEndpoints(s, gidxlist[i]);


		//if photon origin is outside of viewing range dont bother with it
		Sx = g_endpoints.first[0];
                Sy = g_endpoints.first[1];
                Sz = g_endpoints.first[2];

 		Rorigin = sqrt(Sx*Sx + Sy*Sy);
                if(Rorigin > 25.) continue;

		//coordinate geometry time
		//first address only candidates that have valid endpoints
		if( g_endpoints.second[0] != -9999 && g_endpoints.second[1] != -9999 && g_endpoints.second[2] != -9999 ){
		
			Ex = g_endpoints.second[0];
			Ey = g_endpoints.second[1];
			Ez = g_endpoints.second[2];

			Rfinal = sqrt(Ex*Ex + Ey*Ey);
	                tdrift = (GV.ZCUT - Sz)/(Ez - Sz);
        	        xdrift = Sx + tdrift*(Ex - Sx);
               		ydrift = Sy + tdrift*(Ey - Sy);
                	Rdrift = sqrt( xdrift*xdrift + ydrift*ydrift );

		    }//end valid endpoint check
		   else{
			//std::cout<<"No valid endpoint!! -- projecting arbitrary endpoint from origin and momentum vector\n"; 
			//calculate Ex,Ex,Ez through scaling momentum
			Rdrift = getDriftFromPxPyPz(s, gidxlist[i] );
			Rfinal = 25.1; //set rfinal just outside of R "viewing range" because it'll pass through everything if it doesnt drift out in z
			}

		//std::cout<<"photon Rfinal="<<Rfinal<<" Rdrift="<<Rdrift<<std::endl;
		//if Rdrift is less than conversion point, then we drifted out of z .. so set rfinal=rdrift
		if( Rdrift < Rfinal ) Rfinal = Rdrift;
	 	g_eta_physics = SimTrk_eta.At( gidxlist[i] );
		FillTH1(id_gflux_eta, g_eta_physics, w);
		if( abs(g_eta_physics) > 0.1 ) continue; //central cut REMOVING Central cut for graham slides
//		gpt_forcut = SimTrk_pt.At(gidxlist[i] );
		gE = SimTrk_pt.At(gidxlist[i])*cosh( SimTrk_eta.At(gidxlist[i]) );
//		if( gpt_forcut < 1.0) continue;	
			//fill every bin up to the point of conversion or drift	
		for(int j=0; j<rbinedge.size(); j++){
			if( Rfinal >= rbinedge.at(j) ){
				FillTH1(id_gflux, rbincenter.at(j) , w);
				FillTH1(id_gflux_nowt, rbincenter.at(j), 1);
				
				FillTH2(id_reta_ng, g_eta_physics, rbincenter.at(j), w);	
				FillTH2(id_gfluxE, rbincenter.at(j), gE, w);	
			}	
			else{
				break;
			}
		}//end bin edge loop					


	}//end gidxlist loop
	
/*	
	//do flux plots
	std::vector< std::pair<int, int> > GColl = getGParentColl(s);
	std::pair<int,int> coll;
	std::pair<std::vector<double>, std::vector<double> > points;
	//int gidx,type;
	const Int_t nbinz = 5;
        Double_t edgez[nbinz + 1] = {1.0, 5.0, 9.0, 13.5, 18.0, 25.0};
	int type;
	double Sx,Sy,Sz,Ex,Ey,Ez,SR,ER;
	for(int i=0; i<GColl.size(); i++){
		coll = GColl[i];
		gidx = coll.first;
		type = coll.second;
		FillTH1(id_fluxcomp, type, w);
		
		//get this photons endpoints
		points = getGEndpoints(s, gidx);
		//1cm bins in R0-25	
		//for(double binedge=0; binedge<25; binedge=binedge+0.25){
		double binstart,binend;
		for(int I=0; I<nbinz; I++){
			binstart= edgez[I];
			binend= edgez[I+1];
			//does it pass cuts
			//pt>.4? costheta? zpos?
			 gpt = SimTrk_pt[gidx];
			 geta = SimTrk_eta[gidx];
                	 gpz = gpt*sinh(geta);
                         costg = cos( atan2(gpt,gpz) );
			 Sx = points.first[0]; 
			 Sy = points.first[1];
			 Sz = points.first[2];
			 Ex = points.second[0];
			 Ey = points.second[1];
			 Ez = points.second[2];
			if( abs(Sz) < GV.ZCUT && abs(costg) < GV.COSTCUT && gpt > GV.MINPT*2. ){	
				//
				SR = sqrt(Sx*Sx + Sy*Sy);
				ER = sqrt(Ex*Ex + Ey*Ey);

				//if(SR < binedge && Ex==-9999 && Ey==-9999){
				if(SR < binstart && Ex == -9999 && Ey == -9999){
				//	FillTH1(id_radflux, double(binedge)-0.1, w);
					FillTH1(id_radfluxcoarse, binstart+0.1, w);
				}
				//else if(SR < binedge && ER >= binedge && Ex!=-9999 && Ey!=-9999){
				else if(SR < binstart && ER >= binend && Ex!=-9999 && Ey!=-9999){
					FillTH1(id_radfluxcoarse, binstart+0.1, w);
				//	FillTH1(id_radflux,double(binedge)-0.1, w);
				}
			}
			
		}//end bin edge loop	
	}
	
*/
}  // End of sub-program
#endif
