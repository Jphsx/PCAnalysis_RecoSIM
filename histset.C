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

/*	FillTH1( id_numpcHist, numberOfPC,w);

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
*/

	//get sim pc mask
	
//	sim_pc SPC;
//	SPC = GetSimPC(s);

//	std::cout<<"CUBE TEST"<<std::endl;
//	Cube* c;
//	c = buildCube(0, 0, 0, 1);
//	c->PrintCube();

/////////////////////////////////////////////////////////////////////////////
//cartesian box approach
sim_g SG;
SG = GetSimG(s);

//get a set of cubes in R
//choose theta=pi/2 phi=pi/2
//get r from 0,10 in 10 bins
double pi = 4.*atan2(1,1);
std::vector<Cube*> RCubes;
double x,y,z;
double _r,_theta,_phi;
_theta=pi/2.;
_phi = pi/2.;
double L = 1.;
for( double i=10.; i<=100.; i=i+10.0){
	_r =i;
	x = _r* sin(_theta)*cos(_phi);
	y = _r* sin(_theta)*sin(_phi);
	z = _r* cos(_theta);

	RCubes.push_back( buildCube(x,y,z,10.) );
}

std::vector<Cube*> TCubes;
for( double i=0.; i<=pi; i=i+0.01){
	_r = 50.;
	x = _r* sin(i)*cos(_phi);
	y = _r* sin(i)*sin(_phi);
	z = _r* cos(i);

	TCubes.push_back( buildCube(x,y,z,10.) );
}

double pt,eta,phi,P,gx,gy,gz,R_origin;
double ex,ey,ez,R_end;
bool validEndpoint=false;
double px,py,pz;
int gptype;
double Rcube,Rcubex,Rcubey,Rcubez;
//double Rcube_xy;
double Tcube,Ecube;
for(int i=0; i< SG.g_idx.size(); i++){
		//std::cout<<"N photons : "<< SG.g_idx.size()<<std::endl;
		//std::cout<<"photon "<<i<<std::endl;
		//std::cout<<"|gidx "<<SG.g_idx[i]<<" |";
		//std::cout<<"|pvtxidx "<< SG.parent_simvtx_idx[i]<<" |";
		//std::cout<<"|cvtxidx "<< SG.endpoint_vtx_idx[i]<<" |";
		//std::cout<<"|ptype "<< SG.parent_simvtx_ptype[i]<<" |";
		//std::cout<<"|ppdg "<<SG.parent_pdg[i]<<" |"<<std::endl; 
		//std::cout<<"P"<<
		
		
		pt = SimTrk_pt[ SG.g_idx[i] ];
		eta = SimTrk_eta[ SG.g_idx[i] ];
		phi = SimTrk_phi[ SG.g_idx[i] ];
		P = pt*cosh(eta);
		gx = SimVtx_x[ SG.parent_simvtx_idx[i] ];
		gy = SimVtx_y[ SG.parent_simvtx_idx[i] ];
		gz = SimVtx_z[ SG.parent_simvtx_idx[i] ];	
		R_origin = sqrt( gx*gx + gy*gy + gz*gz );

		if( SG.endpoint_vtx_idx[i] != -1){
			validEndpoint=true;
			ex = SimVtx_x[ SG.endpoint_vtx_idx[i] ];
			ey = SimVtx_y[ SG.endpoint_vtx_idx[i] ];
			ez = SimVtx_z[ SG.endpoint_vtx_idx[i] ];
		}
		else{
			validEndpoint=false;
			R_end = -1;
			ex = 0;
			ey = 0;
			ez = 0;
		}

		R_end= sqrt( ex*ex + ey*ey + ez*ez );

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);	

		gptype = SG.parent_simvtx_ptype[i];

		//std::cout<<"P,pt,eta,phi "<<P<<" "<<pt<<" "<<eta<<" "<<phi<<std::endl;
		//std::cout<<"R_org x,y,z "<<R_origin<<" "<<gx<<" "<<gy<<" "<<gz<<std::endl;
		//std::cout<<"R_end x,y,z "<<R_end<<" "<<ex<<" "<<ey<<" "<<ez<<std::endl;

		
		for(int j=0; j<RCubes.size(); j++){
			//make sure same direction as cube
			
			//if(R_origin<1 && hasCubeIntersection(RCubes[j], gx, gy, gz,  px, py, pz )){
			//dot product must be positive and cube needs to be further out in r		
			Rcubex = RCubes[j]->center->x;
			Rcubey = RCubes[j]->center->y;
			Rcubez = RCubes[j]->center->z;
			Rcube = RCubes[j]->center->r;
			if( Rcube > R_origin ){
			if( validEndpoint && Rcube < R_end ){
			if( (Rcubex*px + Rcubey*py + Rcubez*pz) > 0 ){
			if( hasCubeIntersection(RCubes[j], gx, gy, gz,  px, py, pz ) ){
				//Rcube = RCubes[j]->centerRTP->r;
			
				FillTH1(id_Rgeom1_g, Rcube);
				if(gptype == 0){
				    FillTH1(id_Rgeom1_g0, Rcube, Rcube);
				}
				if(gptype == 201){
				    FillTH1(id_Rgeom1_g201, Rcube, Rcube);
				}
				if(gptype == 3){			
				    FillTH1(id_Rgeom1_g3, Rcube, Rcube);
				}

			}}}}//end r and intersection checks


		}//end rcubes loop

		//angular stuff

		for(int j=0; j<TCubes.size(); j++){
			Rcubex = TCubes[j]->center->x;
			Rcubey = TCubes[j]->center->y;
			Rcubez = TCubes[j]->center->z;
			Rcube = TCubes[j]->center->r;
			Tcube = TCubes[j]->center->theta;
			Ecube = TCubes[j]->center->eta;
			//std::cout<<"Tcube "<<Tcube<<" "<<Ecube<<std::endl;

			if( Rcube > R_origin ){
			if( validEndpoint && Rcube < R_end ){
			if( (Rcubex*px + Rcubey*py + Rcubez*pz) > 0 ){
			if( hasCubeIntersection(TCubes[j], gx, gy, gz,  px, py, pz ) ){
				//Rcube = RCubes[j]->centerRTP->r;
			
				FillTH1(id_thgeom1_g, Ecube);
				FillTH1(id_thgeom1a_g, cos(Tcube), sin(Tcube));
				if(gptype == 0){
				    FillTH1(id_thgeom1_g0, Ecube);
				    FillTH1(id_thgeom1a_g0, cos(Tcube), sin(Tcube));
				}
				if(gptype == 201){
				    FillTH1(id_thgeom1_g201, Ecube);
				    FillTH1(id_thgeom1a_g201, cos(Tcube), sin(Tcube) );
				}
				if(gptype == 3){			
				    FillTH1(id_thgeom1_g3, Ecube);
				    FillTH1(id_thgeom1a_g3, cos(Tcube), sin(Tcube));
				}
			}}}}

		}//end tcubes loop


}




//////////////////////////////////////////////////////////////////////////
//spherical approach
/****
	sim_g SG;
	SG =GetSimG(s);
	double R2,R1;
	double P,pt,eta,phi;
	double gx,gy,gz, R_origin;//photon creation point	
	//set up radial flux parameters
	double A=1.;
	double pi=4.*atan2(1,1);
	double phiStar = 0.0;
	
	std::vector<AngleLimits> RLims;
	for(double i=1.; i<=10.; i=i+1.0){
		RLims.push_back( GetAngleLimits(i,A,pi/2.));
		/*
		std::cout<<"RLims"<<std::endl;
		AngleLimits atest = GetAngleLimits(i,A,pi/2.);
		std::cout<<"(R, A, thetaStar): "<< atest.R<<" "<<atest.A<<atest.thetaStar<<std::endl;
		std::cout<<"(dsi, etaup, etalo) "<< atest.dphi <<" "<<atest.etaup<<" "<<atest.etalow<<std::endl; */


/*****

	}
	std::vector<AngleLimits> thLims;
	for(double i=0.; i<pi; i=i+0.1){//looop in theta but convert to eta later
		thLims.push_back( GetAngleLimits(2.0, A, i));

		/*std::cout<<"thLims"<<std::endl;
		AngleLimits atest = GetAngleLimits(2.0,A,i);
		std::cout<<"(R, A, thetaStar, etastar): "<< atest.R<<" "<<atest.A<<" "<<atest.thetaStar<<" "<<atest.etaStar<<std::endl;
		std::cout<<"(dsi, etaup, etalo) "<< atest.dphi <<" "<<atest.etaup<<" "<<atest.etalow<<std::endl;
 		*/
/****
	}

	double gptype;


	//thetastar variation
	

	for(int i=0; i< SG.g_idx.size(); i++){
		std::cout<<"N photons : "<< SG.g_idx.size()<<std::endl;
		std::cout<<"photon "<<i<<std::endl;
		std::cout<<"|gidx "<<SG.g_idx[i]<<" |";
		std::cout<<"|pvtxidx "<< SG.parent_simvtx_idx[i]<<" |";
		std::cout<<"|cvtxidx "<< SG.endpoint_vtx_idx[i]<<" |";
		std::cout<<"|ptype "<< SG.parent_simvtx_ptype[i]<<" |";
		std::cout<<"|ppdg "<<SG.parent_pdg[i]<<" |"<<std::endl; 
		//std::cout<<"P"<<
		
		
		pt = SimTrk_pt[ SG.g_idx[i] ];
		eta = SimTrk_eta[ SG.g_idx[i] ];
		phi = SimTrk_phi[ SG.g_idx[i] ];
		P = pt*cosh(eta);
		gx = SimVtx_x[ SG.parent_simvtx_idx[i] ];
		gy = SimVtx_y[ SG.parent_simvtx_idx[i] ];
		gz = SimVtx_z[ SG.parent_simvtx_idx[i] ];	
		R_origin = sqrt( gx*gx + gy*gy + gz*gz );		

		std::cout<<"P,pt,eta,phi "<<P<<" "<<pt<<" "<<eta<<" "<<phi<<std::endl;
		std::cout<<"R_org x,y,z "<<R_origin<<" "<<gx<<" "<<gy<<" "<<gz<<std::endl;

		FillTH1(id_gpt, pt, w);
		FillTH1(id_gP, P, w);
		FillTH1(id_geta, eta, w);
		//loop over Rlimits
		gptype = SG.parent_simvtx_ptype[i];
		for(int j=0; j< RLims.size(); j++){
			
			R2 = RLims[j].R * RLims[j].R;
			R1 = RLims[j].R;
			//if( phi>= phiStar-RLims[j].dphi && phi<= phiStar+RLims[j].dphi && eta >= RLims[j].etalow && eta <= RLims[j].etaup ){
			if( eta >= RLims[j].etalow && eta <= RLims[j].etaup ){
						
			//make sure photons are created before area segment
			//if( R_origin < RLims[j].R ){
			if(R_origin < 1.0){
			//if(eta >= RLims[j].etalow && eta <= RLims[j].etaup){
				//we flux through the area
				FillTH1(id_Rgeom1_g, RLims[j].R, R1);
				if(gptype == 0){
				    FillTH1(id_Rgeom1_g0, RLims[j].R, R1);
				}
				if(gptype == 201){
				    FillTH1(id_Rgeom1_g201, RLims[j].R, R1);
				}
				if(gptype == 3){			
				    FillTH1(id_Rgeom1_g3, RLims[j].R, R1);
				}
			}//end rorigin check
			}//end eta check
			
		}//end rlim loop
		for(int j=0; j< thLims.size(); j++){
			if( eta >= thLims[j].etalow && eta <= thLims[j].etaup ){
			if(R_origin < 1.0){
				FillTH1(id_thgeom1_g, thLims[j].etaStar, w);
				FillTH1(id_thgeom1a_g, thLims[j].thetaStar, w);
				if(gptype == 0){
				    FillTH1(id_thgeom1_g0, thLims[j].etaStar, w);
				    FillTH1(id_thgeom1a_g0, thLims[j].thetaStar, w);
				}
				if(gptype == 201){
				    FillTH1(id_thgeom1_g201, thLims[j].etaStar, w);
				    FillTH1(id_thgeom1a_g201, thLims[j].thetaStar, w);
				}
				if(gptype == 3){			
				    FillTH1(id_thgeom1_g3, thLims[j].etaStar, w);
				    FillTH1(id_thgeom1a_g3, thLims[j].thetaStar, w);
				}		
			}
			}
		}

		
		if(SG.parent_simvtx_ptype[i] == 0){
			FillTH1(id_g0pt, pt, w);
			FillTH1(id_g0P, P, w);
			FillTH1(id_g0eta, eta, w);
		}
		if(SG.parent_simvtx_ptype[i] == 201){
			FillTH1(id_g201pt, pt, w);
			FillTH1(id_g201P, P, w);
			FillTH1(id_g201eta, eta, w);
		}
		if(SG.parent_simvtx_ptype[i] == 3){
			FillTH1(id_g3pt, pt, w);
			FillTH1(id_g3P, P, w);
			FillTH1(id_g3eta, eta, w);
		}
		
		
	}
	
/*
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
*/	

}  // End of sub-program
#endif
