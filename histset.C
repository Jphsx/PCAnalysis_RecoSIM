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

	 int valid_nPV=0;
        for(int v=0; v< nPV_mask.size(); v++){
                if(nPV_mask.at(v)) valid_nPV++;
        }


	//FillTH1(id_npv, nPV,w);
	FillTH1(id_npv,valid_nPV);
	
	double sumndof=0;
	double tkwsum=0;
	double PVR=-1.;
	for(int v = 0; v<nPV; v++){
		if(nPV_mask.at(v)){

		sumndof+= PV_ndof[v];
		FillTH1(id_pvz, PV_Z[v],w);
		//lets check out PV_R distribution
		PVR= sqrt(PV_X[v] *PV_X[v] + PV_Y[v]*PV_Y[v]);
		FillTH1(id_pvr, PVR,w);
		if(PV_ndof[v] >= 4){
		FillTH1(id_ndofall, PV_ndof[v],w);
		tkwsum = (PV_ndof[v] + 3.)/2.;
		FillTH1(id_pvtrk, tkwsum, w);
		}
		//if(nPV == 1 && PV_ndof[v] >= 4){
		//if(valid_nPV == 1 && PV_ndof[v] >= 4){
		if(valid_nPV == 1){
		FillTH1(id_ndofpv1, PV_ndof[v], w);
		tkwsum =(PV_ndof[v] + 3.)/2.;
		//FillTH1(id_pvtrkpv1, tkwsum, w); //fill this after cuts for weighting ref dist
		}
		//if(nPV > 1 && v>0){
		if(valid_nPV > 1 && v>0){
			if( nPV_mask.at(v))
			FillTH1(id_pvdz, PV_Z[0]-PV_Z[v], w);
		}
		}//END NPVMASK CHECK
	}
	if(sumndof >= 4)
	FillTH1(id_sumndof, sumndof,w);
 
       // if(nPV != 1 || !nPV_mask[0] || PV_ndof[0]<4 || PV_ndof[0]>50) return;
        //if(nPV != 1 || !nPV_mask[0] || PV_ndof[0]<4 || tkwsum>40) return;
//	if(nPV != 1 || !nPV_mask[0] || PV_ndof[0]<4 || tkwsum>20) return; //trying the tight cut
       // if(nPV != 1 || !nPV_mask[0] || PV_ndof[0]<4 || tkwsum>20 || tkwsum<8 ) return; //tight cut with lower cut
	if(valid_nPV != 1 || PV_ndof[0]<4 || tkwsum>20 || tkwsum<8 ) return;
	//if(nPV != 1 || !nPV_mask[0] || tkwsum < 9 || tkwsum > 40) return;
	
	FillTH1(id_pvtrkpv1, tkwsum, w);
        FillTH1(id_PVndof, PV_ndof[0],w);			
/////EVERYTHING HAS NPV =1 and tksum cut EMBEDDED
	
	//find the pvz and pvr of the first and only valid pv, load this into pv cut for weighting
	for(int v =0; v<nPV; v++){
		if(nPV_mask.at(v)){
			PVR= sqrt(PV_X[v] *PV_X[v] + PV_Y[v]*PV_Y[v]);
			FillTH1(id_pvr_cut, PVR,w);
			FillTH1(id_pvz_cut, PV_Z[v],w);	
		}
	}		


	//plot "raw" conv stuff
	FillTH1( id_numpcHist, numberOfPC,w);
std::vector<CommonVars> CVs = GetCommonVars(s,true);

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
		FillTH1(id_rho25RHist, CVs[i].rho,w);
		FillTH1(id_rps25RHist, CVs[i].rps,w);
		FillTH1(id_etaHist, CVs[i].etaphys, w);
		FillTH2(id_ndof_pcReta, PV_ndof[0] ,CVs[i].etaphys,w);
		FillTH2(id_ndof_pcRpt, PV_ndof[0], PC_pt,w);
		FillTH1(id_Rtheta, CVs[i].theta, w);
		FillTH1(id_Rphi, CVs[i].phi, w);
        }
	
	//calculate common variables
//        std::vector<CommonVars> CVs = GetCommonVars(s,true);

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
			FillTH1(id_rho25CHist, CVs[i].rho,w);
                	FillTH1(id_rps25CHist, CVs[i].rps,w);
			FillTH1(id_etaCutHist, CVs[i].etaphys, w);
                }

        }
        FillTH1(id_numpccutHist, npcCut, w);

	//fill "hgn" conv stuff

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
         //       FillTH2( id_xywideHGNPCHist, CVs[cidx].x, CVs[cidx].y, w);
         	FillTH1(id_r25HHist , CVs[cidx].radius,w);
		FillTH1(id_rho25HHist, CVs[cidx].rho,w);
                FillTH1(id_rps25HHist, CVs[cidx].rps,w);
		FillTH1(id_etaHGNHist, CVs[cidx].etaphys,w);
		FillTH2(id_ndof_pcHeta, PV_ndof[0] ,CVs[cidx].etaphys,w);
                FillTH2(id_ndof_pcHpt, PV_ndof[0], PC_pt,w);
		FillTH1(id_Htheta, CVs[cidx].theta, w);
                FillTH1(id_Hphi, CVs[cidx].phi, w);
		FillTH1(id_r25HHist_5bin, CVs[cidx].radius, w);		

		if(CVs[cidx].radius < 8)
		FillTH1( id_sumtksum_rlo, tkwsum, w);

		if(CVs[cidx].radius > 8 && CVs[cidx].radius < 25)
		FillTH1( id_sumtksum_rhi, tkwsum, w);

		FillTH2( id_xyHistHGN, PC_x[cidx], PC_y[cidx], w);
                FillTH2( id_xywideHGN, PC_x[cidx], PC_y[cidx], w);

		FillTH1(id_rps25_b2p5, CVs[cidx].rps,w);
		FillTH1(id_rho25_b2p5, CVs[cidx].rho,w);
		FillTH1(id_r25_b2p5, CVs[cidx].radius,w);
		FillTH2(id_rpt, CVs[cidx].radius, PC_pt, w);
        }


}  // End of sub-program
#endif
