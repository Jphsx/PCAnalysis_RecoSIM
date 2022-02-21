#include <iomanip>
#include "Hungarian.h"
#include "TMath.h"
#include <set>

//////////////////////////////////////////////////////////////
//
//(1) hungarian disambiguation
//
struct hgn_pc{
	std::vector<int> vsel; //vector indices of selected conversions after removing duplicates
	int nedges;
	int Nm;
	int Np;
	int nm;
	int np;	

};
hgn_pc pc_disambiguation(recosim& s, std::vector<bool> cutmask){

    const bool lpr = false;       // print flag
    const bool lreduce = true;   // do problem reduction
    const bool lassign = true;   // Do assignment problem

    auto& PC_vTrack0_charge = s.Conv_Tk0_charge;
    auto& PC_vTrack1_charge = s.Conv_Tk1_charge;		
    auto& PC_x = s.Conv_vtx_X;
    auto& Tk0_chi2 = s.Conv_Tk0_chi2;
    auto& Tk1_chi2 = s.Conv_Tk1_chi2;
    auto& Tk0_ndof = s.Conv_Tk0_ndof;
    auto& Tk1_ndof = s.Conv_Tk1_ndof;
    auto numberOfPC = *(s.nConv); 
    auto& PC_vtx_chi2 = s.Conv_vtx_chi2; 
    auto eventNumber = *(s.event);
 
 //fill the cutmask elsewhere and pass it here, then we dont need to hardcode cuts here
    std::vector<bool> vcuts;
    std::vector<double> PosTkInfo;
    std::vector<double> NegTkInfo;
   // std::vector<double> PosPt;
   // std::vector<double> NegPt;
    std::vector<int> vcandidate;
    std::map<double, int> mNeg;
    std::multimap<double, int> mmNeg;
    std::map<double, int> mPos;
    std::multimap<double, int> mmPos;
    std::multimap<int, double> mmCandidate;
    std::set <std::pair<double,double> > trkPair;
    std::vector<int> vsel;   // Vector for indices of selected conversions after removing duplicates

    //replacing tup with just 1 fit probability vec
    std::vector<double> tup;
    int nassigned = 0;
    //only do disambiguation for conversions that pass set of cuts that defines cutmask
     //cut mask will be the size numPC, cutmask[i] is true for PC that passes all cuts
    for(int i=0; i<PC_x.GetSize(); i++){
	tup.push_back(-1.0);

        vcuts.push_back(false);
        PosTkInfo.push_back(-1.0);
        NegTkInfo.push_back(-1.0);
    //    PosPt.push_back(-1.0);
    //    NegPt.push_back(-1.0);

	int q0 = PC_vTrack0_charge[i];
        int q1 = PC_vTrack1_charge[i];
	double fitprob = TMath::Prob(PC_vtx_chi2[i], 3);

	   if( cutmask[i] ){
           	vcuts[i] = true;
	// Keep track of charge-signed chisq/dof of constituent tracks to later identify conversion 
	// candidates that use the same track (sign by charge of track)
            	if(q0 == 1 && q1 == -1){
               		PosTkInfo[i] =  Tk0_chi2[i]/double(Tk0_ndof[i]);
               		NegTkInfo[i] = -Tk1_chi2[i]/double(Tk1_ndof[i]);
            //   PosPt[i] = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
            //   NegPt[i] = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
//               PosTkInfo[i] =  abs(Tk0_sd0[i])*Tk0_chi2[i]/double(Tk0_ndof[i]);
//               NegTkInfo[i] = -abs(Tk1_sd0[i])*Tk1_chi2[i]/double(Tk1_ndof[i]);
            	}//end if q0 q1
            	else{
               		PosTkInfo[i] =  Tk1_chi2[i]/double(Tk1_ndof[i]);
               		NegTkInfo[i] = -Tk0_chi2[i]/double(Tk0_ndof[i]);
            //   PosPt[i] = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
            //   NegPt[i] = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
//               PosTkInfo[i] =  abs(Tk1_sd0[i])*Tk1_chi2[i]/double(Tk1_ndof[i]);
//               NegTkInfo[i] = -abs(Tk0_sd0[i])*Tk0_chi2[i]/double(Tk0_ndof[i]);
           	 }//end else
	// Need to make sure that each track pair is distinguishable from pairs already selected
            	if(!trkPair.insert(std::make_pair(NegTkInfo[i], PosTkInfo[i])).second){
               		if(lpr)std::cout << "INDISTINGUISHABLE EDGE " << i << " ignored event "<< eventNumber << std::endl;
	// good idea to keep a tally in some histogram bin
               vcuts[i] = false;
            	}//endif
		else{
	// Fill STL containers that help define the "matching problem"
               vcandidate.push_back(i);
	// For the key/value pair use charge-signed chi2/dof as key, and conversion index as edge id for value
               mNeg.insert(std::make_pair(NegTkInfo[i],i));
               mmNeg.insert(std::make_pair(NegTkInfo[i],i));
               mPos.insert(std::make_pair(PosTkInfo[i],i));
               mmPos.insert(std::make_pair(PosTkInfo[i],i));
	       tup[i] = fitprob;
	       }//end mapping else
	   }//end cutmask if
    }//end pc_x loop 
	

	// Characterize our matching problem for this event
	// Note this is unnecessary if there are no duplicates, ie. n- = n+ = nedges.

       if(std::min(mNeg.size(), mPos.size()) < vcandidate.size() && lassign ){
	// We actually have an assignment problem to worry about and we want to worry about it
       if(lpr){
          std::cout << " " << std::endl;
          std::cout << "Event " << eventNumber << " numberOfPC " << numberOfPC << std::endl;
     //     std::cout << "numberOfPC " << numberOfPC << std::endl;
	  std::cout << "nedges = " << vcandidate.size() << std::endl;
          std::cout << "N- = " << mmNeg.size() << std::endl;
          std::cout << "N+ = " << mmPos.size() << std::endl;
          std::cout << "n- = " << mNeg.size() << std::endl;
          std::cout << "n+ = " << mPos.size() << std::endl;
          std::cout << "Target maximum possible cardinality of solution = " << std::min(mNeg.size(),mPos.size()) << std::endl;
	// Prior solution with no arbitration
          std::cout << "Prior edge solution (no arbitration at all) : ";
          for(int i=0; i<numberOfPC; ++i){
              if(vcuts[i])std::cout << std::setw(3) << i;
          }
          std::cout << std::endl;
       }
	// Let's do some more characterization of the ambiguity complexity by investigating the multimaps
	// with a view to removing the non-ambiguities.
	// Note that in case of count==1 from the multimap, the map value is the unique edge ID for this polarity of track 
       if(lreduce){
          for (auto i = mNeg.begin(); i!= mNeg.end(); ++i) {
              auto tkInfo = i->first;
              auto tkEdge = i->second;
              auto negCount = mmNeg.count(tkInfo);
              if(lpr){
                  std::cout << negCount << " edge(s) for e- with first key,value=("
                            << tkInfo << "," << tkEdge << ") [edges: ";
	// Example from http://www.cplusplus.com/reference/map/multimap/count/
	// equal_range returns a pair with lower_bound and upper_bound positions.
                  for (auto it=mmNeg.equal_range(tkInfo).first; it!=mmNeg.equal_range(tkInfo).second; ++it){
                       std::cout << ' ' << (*it).second;
                  }
                  std::cout << " ] " << std::endl;
              }
              if(negCount==1){
	// Make candidate multimap when only one possibility for this electron. The key is the edge id. 
                 mmCandidate.insert(std::make_pair(tkEdge, tkInfo));
              }
          }//end mNeg auto loop
    	   for (auto i = mPos.begin(); i!= mPos.end(); ++i) {
              auto tkInfo = i->first;
              auto tkEdge = i->second;
              auto posCount = mmPos.count(tkInfo);
              if(lpr){
                  std::cout << posCount << " edge(s) for e+ with first key,value=("
                            << tkInfo << "," << tkEdge << ") [edges: ";
                  for (auto it=mmPos.equal_range(tkInfo).first; it!=mmPos.equal_range(tkInfo).second; ++it){
                       std::cout << ' ' << (*it).second;
                  }
                  std::cout << " ] " << std::endl;
              }
              if(posCount==1){
                 mmCandidate.insert(std::make_pair(tkEdge, tkInfo));
              }
          }//end Mpos auto loop
          if(lpr)std::cout << "mmCandidate multimap size " << mmCandidate.size() << std::endl; 
	  for (auto i = mmCandidate.begin(); i!= mmCandidate.end(); ++i) {
               auto tkEdge = i->first;
               auto tkInfo = i->second;
               if(lpr)std::cout << "mmCandidate key= " << tkEdge << " multiplicity "
                                << mmCandidate.count(tkEdge) << " (value " << tkInfo << " ) " << std::endl;
               if(mmCandidate.count(tkEdge) == 2 ){
// There is an electron-positron pairing where the degree of each vertex is 1, and the same edge is incident on both.
// So if we add this edge to the selected candidates we can erase this one and its constituents from the assignment problem.
                  if(tkInfo<0.0){
                     vsel.push_back(tkEdge);
                     vcandidate.erase(remove(vcandidate.begin(),vcandidate.end(),tkEdge),vcandidate.end());
                     mNeg.erase(tkInfo);
                  }
                  else{
                     mPos.erase(tkInfo);
                  }
               }
          }//end mmcandidate loop
          if(lpr){
              std::cout << "Already selected " << vsel.size() << " non-ambiguous pairings " << std::endl;
              std::cout << "After reduction, nedges = " << vcandidate.size() << std::endl;
              std::cout << "n- = " << mNeg.size() << std::endl;
              std::cout << "n+ = " << mPos.size() << std::endl;

          }
      }  // end of lreduce clause
	
    std::vector< std::vector <double> > costMatrix;
    std::vector< std::vector <int> > edgeMatrix;
    int irow = -1;
    for ( auto i = mNeg.begin(); i != mNeg.end(); ++i ){        // n- rows for each -ve track
         irow++;
         std::vector<double> v;
         std::vector<int> e;
         int nfound = 0;
         for ( auto j = mPos.begin(); j != mPos.end(); ++j ){   // n+ cols for each +ve track
	// Go through all the edge candidates and see if there is an edge corresponding to this pairing.
             double negInfo = i->first;
             double posInfo = j->first;
             bool found = false;
             for (auto iter = vcandidate.begin(); iter != vcandidate.end(); ++iter) {
                 unsigned int k = *iter;
                 if(NegTkInfo[k] == negInfo && PosTkInfo[k] == posInfo ) {
	// Found matching edge
                    if(lpr)std::cout << "Found matching edge " << k << std::endl;
	// Add cost value of corresponding chi-squared value for 1 dof.
                    v.push_back(TMath::ChisquareQuantile(1.0-tup[k],1.0));
                    e.push_back(k);
                    found=true;
                    nfound++;
                 }
             }
             if(!found){
                v.push_back(10000.0);
                e.push_back(-1);
             }
         }
         if(lpr)std::cout << "irow " << irow << " nedges = " << nfound << std::endl;
         costMatrix.push_back(v);
         edgeMatrix.push_back(e);
    }// end mNeg auto loop

   // Debug printing
    if(lpr){
       for (auto irow = costMatrix.begin(); irow != costMatrix.end(); ++irow) {
           std::cout << "Row weights:   ";
           for (auto pos = irow->begin(); pos != irow->end(); ++pos) {
               std::cout << std::setw(10) << *pos << " ";
           }
           std::cout << std::endl;
       }
       for (auto irow = edgeMatrix.begin(); irow != edgeMatrix.end(); ++irow) {
           std::cout << "Edges      :   ";
           for (auto pos = irow->begin(); pos != irow->end(); ++pos) {
               std::cout << std::setw(10) << *pos << " ";
           }
           std::cout << std::endl;
       }
    }// end debug print
    
    // Use Hungarian Algorithm to find the minimum cost assignment of e- to e+. 
    // The total cost will be the total chi-squared of all assignments (each with 1 dof).
    // Need to take care also of the case where there is no one-sided perfect matching, 
    // and algorithmically, the extra assignments are assigned to the fictional 
    // high-cost edges with nominal weight of 10000.
    HungarianAlgorithm HungAlgo;
    std::vector<int> assignment;
    double cost = HungAlgo.Solve(costMatrix, assignment);
    if(lpr)std::cout << "costMatrix.size() " << costMatrix.size() << std::endl;
        for (unsigned int x = 0; x < costMatrix.size(); x++){
    // Need to check if this row is assigned (may not be if nRows > nCols)
        if(assignment[x] >= 0){
           if(edgeMatrix[x][assignment[x]] == -1){
              cost-=10000.0;
           }
           else{
              vsel.push_back(edgeMatrix[x][assignment[x]]);
              nassigned++;
           }
           if(lpr)std::cout << x << "," << assignment[x]
                     << " " << costMatrix[x][assignment[x]] << " "
                     << edgeMatrix[x][assignment[x]] << std::endl;
        }
    }
    if(nassigned > 0){
       if(lpr)std::cout << "Assigned " << nassigned << " initially ambiguous pairings" << std::endl;
       double psel = TMath::Prob(cost, nassigned);
           if(lpr)std::cout << "Minimized total chisq: " << cost << " ( " << nassigned << " ) " << " p-value " << psel <<std::endl;
    }
    if(vsel.size() > 0){
       std::sort(vsel.begin(), vsel.end());
       if(lpr){
          std::cout << "Selected conversions (" << vsel.size() << ")" ;
          for(int i=0; i<vsel.size(); ++i){
              std::cout << std::setw(3) << " " << vsel[i];
          }
          std::cout << std::endl;
       }
    }
 }
  else{
// No ambiguities - so no assignment problem to solve 
     if(vcandidate.size() > 0){
        vsel = vcandidate;
        std::sort(vsel.begin(), vsel.end());
        if(lpr){
           std::cout << " " << std::endl;
           std::cout << "Selected conversions S: " << eventNumber << "  (" << vsel.size() << ")";
           for(int i=0; i<vsel.size(); ++i){
              std::cout << std::setw(3) << " " << vsel[i];
           }
           std::cout << std::endl;
        }
     }
  }


    hgn_pc HGN;
     HGN.vsel= vsel; //vector indices of selected conversions after removing duplicates
    //    int nedges;
    //    int Nm;
     //   int Np;
     //   int nm;
     //   int np;
    return HGN;

}
/////////////////////////////////////////////////////////////////////////////////////////////
//(2) sim identification
//
struct sim_photon{
//	double px;
//	double py;
//	double pz;'
	double pt;
	double eta_physics;
	double phi;
	
};
sim_photon GetSimGfromTID(recosim& s, int tid){
	sim_photon SG;
	auto& SimTrk_trackId = s.SimTrk_trackId;
	auto& SimTrk_eta = s.SimTrk_eta;
	auto& SimTrk_phi = s.SimTrk_phi;
	auto& SimTrk_pt = s.SimTrk_pt;
	SG.pt=0.;
	SG.eta_physics=0.;
	SG.phi=0.;
	for(int i=0; i<SimTrk_trackId.GetSize(); i++){
		if(SimTrk_trackId.At(i) == tid){
			SG.pt = SimTrk_pt.At(i);
			SG.eta_physics = SimTrk_eta.At(i);
			SG.phi = SimTrk_phi.At(i);
			return SG;
		}
	} 

	return SG;
		

}
struct sim_pc{
	
	std::vector<int> sim_mask;
	std::vector<int> p14_t1;
	std::vector<int> p14_t2;
	std::vector<int> p14_g;
	std::vector<int> p14_key;//list of simvtx indices (intended to access conversions without fully looping sim vtx and checking mask
	
	

};


//Return struct of sim mask, currently masks process 14 and finds associated e/p pair and returns indices
sim_pc GetSimPC(recosim& s){
	sim_pc SPC;

	int nSimVtx = (s.SimVtx_processType).GetSize();
        int nSimTrk = (s.SimTrk_simvtx_Idx).GetSize();

	auto& SimVtx_processType = s.SimVtx_processType;
	auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
    	auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;
    	auto& SimTrk_trackId = s.SimTrk_trackId;

	std::vector<int> childholder;
	std::vector<int> sim_mask(nSimVtx); //contains a 14 if simvtx[i] is a pc
	std::vector<int> p14_t1(nSimVtx); //contains the index of simtrack of first child of simvtx[i]
	std::vector<int> p14_t2(nSimVtx); // contains the index of simtrack of the second child of simvtx[i]	
	std::vector<int> p14_g(nSimVtx); //contains the index of simtrack of parent photon to simvtx[i]
	std::vector<int> p14_key;//list of simvtx indices
	
	int numchild=0;
	int ptid;
	int ptid_idx;
	for(int i=0; i<nSimVtx; i++){
 	 	sim_mask[i]= -1;
  		p14_t1[i] = -1;
 		p14_t2[i] = -1;
		p14_g[i] = -1;
		ptid = SimVtx_simtrk_parent_tid[i];
  	if( SimVtx_processType[i] != 14){
        	continue;
  	}
  	else{
        	numchild=0;
		ptid_idx=-1;
        	for(int j=0; j<nSimTrk; j++){
                	if( SimTrk_simvtx_Idx[j] == i){
                        	numchild++;
                        	childholder.push_back(j);
                	}
			if( SimTrk_trackId[j] == ptid ){
				ptid_idx = j;
			}
			//if you find everything stop looping
			if(numchild == 2 && ptid_idx != -1) break;
               	/*	if(numchild == 2){
                        	vtxmask[i] = true;
                        	vtxc1[i] = childholder[0];
                        	vtxc2[i] = childholder[1];
                        	childholder.clear();
				
                        	break;
                	}*/
			
		}
		if(numchild == 2){
			sim_mask[i] = 14;
			p14_t1[i] = childholder[0];
			p14_t2[i] = childholder[1];
			p14_g[i] = ptid_idx;	
			p14_key.push_back(i);	
		}
       		childholder.clear();

  	}	

}
	SPC.sim_mask = sim_mask;
	SPC.p14_t1 = p14_t1;
	SPC.p14_t2 = p14_t2;
	SPC.p14_g = p14_g;
	SPC.p14_key = p14_key;
	return SPC;

}

////////////////////////////////////////////////////////////////////
//
//(3) calculation of common variables PER conversion
//
struct CommonVars{
	
    double radius;
    double rerr;
    double pfit;
    double phierr;
    double zerr;
    double rps;
    double rnominal;
    double phi;
    double rho;
    double pt;
    double E;
    double phip;
    double xplus;
    double theta;
    double etaphys;
    double x;
    double y;
    double z;
    int id;

    //track 0and1 params swim from inner vtx to pc vtx
    double q0;
    double q1;
    double px0p;
    double py0p;
    double pz0p;
    double x0p; 
    double y0p; 
    double z0p; 

    double px1p;
    double py1p;
    double pz1p;
    double x1p;
    double y1p;
    double z1p;    
	
};
struct GlobalValues{
    const double RERRCUT = 0.25;
    const double COSTCUT = 0.85;
    const double ZCUT = 25.0;
    const double FITPROBCUT = 0.010;
    const double MASSCUT = 0.15;
    const double MINPT = 0.2;

    const double MASS_ELECTRON = 0.5109989461e-3;
    const double MASS_PION = 139.57061e-3;
    const double MASS_KAON = 493.677e-3;
    const double MASS_PROTON = 938.272081e-3;

    // We now have various "centers" to compare to for radial coordinates.
// beam pipe displacement (in cm) from Anna's DPF2019 talk
    const double x0bpdata =  0.171;
    const double y0bpdata = -0.176;
    const double x0bpmc = 0.0;
    const double y0bpmc = 0.0;

    //BPIX center displacement (in cm) (Anna's 16-Dec-2019 talk page 3)
    const double x0data =  0.086;
    const double y0data = -0.102;
    const double x0mc = 0.0;
    const double y0mc = 0.0;

    //Pixel support displacement (in cm) (Anna's December talk page 4)
    const double x0psdata = -0.080;
    const double y0psdata = -0.318;
    const double x0psmc = 0.0;
    const double y0psmc = 0.0;

};
//////////////////GLOBAL DECLARATION OF CUTS AND OTHER CONST VALUES////////
GlobalValues GV;
/////////////////////////////////////////////////////////////////////////

std::vector<CommonVars> GetCommonVars(recosim& s, bool isRealData){//pass in bool for mc or data, this isn't stored in tree, keep track manually in histset analysis


    auto numberOfPC = *(s.nConv);  
    auto& PC_x = s.Conv_vtx_X;
    auto& PC_y = s.Conv_vtx_Y;
    auto& PC_z = s.Conv_vtx_Z;
    auto& PC_vtx_chi2 = s.Conv_vtx_chi2;
    auto& PC_Px = s.Conv_refittedPair4Momentum_Px;
    auto& PC_Py = s.Conv_refittedPair4Momentum_Py;
    auto& PC_Pz = s.Conv_refittedPair4Momentum_Pz;
    auto& PC_E = s.Conv_refittedPair4Momentum_E;
    auto& PC_vtx_sigmaxx = s.Conv_vtx_cov_00;
    auto& PC_vtx_sigmaxy = s.Conv_vtx_cov_01;
    auto& PC_vtx_sigmayy = s.Conv_vtx_cov_11;
    auto& PC_vtx_sigmazz = s.Conv_vtx_cov_22;

    std::vector<CommonVars> pc_comm(numberOfPC);

    double x0bp,y0bp,x0,y0,x0ps,y0ps;
    if(isRealData){
       x0bp = GV.x0bpdata;
       y0bp = GV.y0bpdata;
       x0 = GV.x0data;
       y0 = GV.y0data;
       x0ps = GV.x0psdata;
       y0ps = GV.y0psdata;
    }
    else{// MC  
       x0bp = GV.x0bpmc;
       y0bp = GV.y0bpmc;
       x0 = GV.x0mc;
       y0 = GV.y0mc;
       x0ps = GV.x0psmc;
       y0ps = GV.y0psmc;
    }
    double x,y,z,pt,pz,E,r;	
    double vxx, vxy, vyy, vzz;
    double phi,theta,eta,rho,phip,rps,rnominal;
    double cphi,sphi,varsum_r,varsum_phi, rerr,phierr,zerr,fitprob;

    double pt0,tanl0,phi0,qR0,px0p,py0p,pz0p,x0p,y0p,z0p,q0;
    double pt1,tanl1,phi1,qR1,px1p,py1p,pz1p,x1p,y1p,z1p,q1;
    for(int i=0; i<numberOfPC; i++){	
	CommonVars CVi;
	x=PC_x[i];
	y=PC_y[i];
	z=PC_z[i];
	pt = sqrt(PC_Px[i]*PC_Px[i] + PC_Py[i]*PC_Py[i]);
	pz = PC_Pz[i];
   	r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
  	E = PC_E[i];
        vxx = PC_vtx_sigmaxx[i];
        vxy = PC_vtx_sigmaxy[i];
        vyy = PC_vtx_sigmayy[i];
        vzz = PC_vtx_sigmazz[i];	
   
	phi = atan2(y-y0, x-x0);
        if (phi < 0) { phi += 2 * M_PI; }
        theta = atan2(pt,pz);
        eta = -log(tan(theta/2.0));
        rho  =  sqrt( (x-x0bp)*(x-x0bp) + (y-y0bp)*(y-y0bp)) ;
        phip =  atan2(y-y0bp, x-x0bp);
        rps = sqrt( (x-x0ps)*(x-x0ps) + (y-y0ps)*(y-y0ps) );
        rnominal = sqrt( x*x + y*y );
        vxx = PC_vtx_sigmaxx[i];
        vxy = PC_vtx_sigmaxy[i];
        vyy = PC_vtx_sigmayy[i];
        vzz = PC_vtx_sigmazz[i];
        cphi = cos(phi);
        sphi = sin(phi);
	// This is the correct one
        varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
        varsum_phi = sphi*sphi*vxx - 2.0*sphi*cphi*vxy + cphi*cphi*vyy;
        rerr = sqrt(varsum_r);
        phierr = sqrt(varsum_phi)/r;
        zerr = sqrt(vzz);
        fitprob = TMath::Prob(PC_vtx_chi2[i], 3);
	
	CVi.radius = r;
	CVi.rerr= rerr;
	CVi.pfit = fitprob;
	CVi.phierr = phierr;
	CVi.zerr = zerr;	
        CVi.rps= rps;
	CVi.rho= rho;
	CVi.phip = phip;
	CVi.rnominal = rnominal;
	CVi.phi = phi;
        CVi.theta = theta;
	CVi.pt = pt;
	CVi.E = E;
	CVi.etaphys = eta;
	CVi.id = i;
	CVi.x = x;
	CVi.y = y;
	CVi.z = z;

        auto& Tk0_px = s.Conv_tracksPin_Px_Tk0;
    	auto& Tk0_py = s.Conv_tracksPin_Py_Tk0;
    	auto& Tk0_pz = s.Conv_tracksPin_Pz_Tk0;
    	auto& Tk0_x  = s.Conv_tracksInnerPosition_X_Tk0;
    	auto& Tk0_y  = s.Conv_tracksInnerPosition_Y_Tk0;
    	auto& Tk0_z  = s.Conv_tracksInnerPosition_Z_Tk0;

    	auto& Tk1_px = s.Conv_tracksPin_Px_Tk1;
    	auto& Tk1_py = s.Conv_tracksPin_Py_Tk1;
    	auto& Tk1_pz = s.Conv_tracksPin_Pz_Tk1;
    	auto& Tk1_x  = s.Conv_tracksInnerPosition_X_Tk1;
    	auto& Tk1_y  = s.Conv_tracksInnerPosition_Y_Tk1;
    	auto& Tk1_z  = s.Conv_tracksInnerPosition_Z_Tk1;

        auto& PC_vTrack0_charge = s.Conv_Tk0_charge;
    	auto& PC_vTrack1_charge = s.Conv_Tk1_charge;

	//// Swim each track from inner hit to the conversion vertex - see convsel.C code
        pt0 = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
        tanl0 = Tk0_pz[i]/pt0;
        phi0 = atan2(Tk0_py[i],Tk0_px[i]);
        qR0 = double(PC_vTrack0_charge[i])*100.0*pt0/(0.2998*3.80); //in cm
        double A0 =  2.0*qR0*( (Tk0_y[i]-y)*cos(phi0) - (Tk0_x[i]-x)*sin(phi0) - qR0 );
        double B0 = -2.0*qR0*( (Tk0_y[i]-y)*sin(phi0) + (Tk0_x[i]-x)*cos(phi0) );
        double alp0 = atan2(B0/A0,1.0);
        px0p = pt0*cos(phi0 + alp0);
        py0p = pt0*sin(phi0 + alp0);
        pz0p = Tk0_pz[i];
	// Also swim the track position
        x0p = Tk0_x[i] + qR0*( (1.0-cos(alp0))*sin(phi0) - sin(alp0)*cos(phi0) );
        y0p = Tk0_y[i] - qR0*( (1.0-cos(alp0))*cos(phi0) + sin(alp0)*sin(phi0) );
        z0p = Tk0_z[i] - qR0*tanl0*alp0;
       

	// Swim each track from inner hit to the conversion vertex - see convsel.C code
        pt1 = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
        tanl1 = Tk1_pz[i]/pt1;
        phi1 = atan2(Tk1_py[i],Tk1_px[i]);
        qR1 = double(PC_vTrack1_charge[i])*100.0*pt1/(0.2998*3.80); //in cm
        double A1 =  2.0*qR1*( (Tk1_y[i]-y)*cos(phi1) - (Tk1_x[i]-x)*sin(phi1) - qR1 );
        double B1 = -2.0*qR1*( (Tk1_y[i]-y)*sin(phi1) + (Tk1_x[i]-x)*cos(phi1) );
        double alp1 = atan2(B1/A1,1.0);
        px1p = pt1*cos(phi1 + alp1);
        py1p = pt1*sin(phi1 + alp1);
        pz1p = Tk1_pz[i];
	// Also swim the track position
        x1p = Tk1_x[i] + qR1*( (1.0-cos(alp1))*sin(phi1) - sin(alp1)*cos(phi1) );
        y1p = Tk1_y[i] - qR1*( (1.0-cos(alp1))*cos(phi1) + sin(alp1)*sin(phi1) );
        z1p = Tk1_z[i] - qR1*tanl1*alp1;
        
	CVi.q0 = PC_vTrack0_charge[i];
        CVi.q1 = PC_vTrack1_charge[i];
        CVi.px0p = px0p;
	CVi.py0p = py0p;
        CVi.pz0p = pz0p;
	CVi.x0p = x0p;
  	CVi.y0p = y0p;
	CVi.z0p = z0p;

	CVi.px1p = px1p;
        CVi.py1p = py1p;
        CVi.pz1p = pz1p;
        CVi.x1p = x1p;
        CVi.y1p = y1p;
        CVi.z1p = z1p;

	double ptasym = (pt0-pt1)/(pt0+pt1);
        if (q0<0) ptasym = -ptasym;
        double xplus = (1.0 + ptasym)/2.0;
        CVi.xplus = xplus;



	pc_comm[i] = CVi;



    }//end loop over nconvs

return pc_comm;
}
double getDriftFromPxPyPz( recosim& s, int gidx){

	auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
	int svidx = SimTrk_simvtx_Idx[gidx];

        auto& SimVtx_x = s.SimVtx_x;
        auto& SimVtx_y = s.SimVtx_y;
        auto& SimVtx_z = s.SimVtx_z;
	double sx,sy,sz;
        sx = SimVtx_x[svidx];
        sy = SimVtx_y[svidx];
        sz = SimVtx_z[svidx];

        auto& SimTrk_pt = s.SimTrk_pt;
        auto& SimTrk_eta = s.SimTrk_eta;
	auto& SimTrk_phi = s.SimTrk_phi;
        double gpt,geta,gphi,gpz,gpx,gpy;


	gpt = SimTrk_pt[gidx];
        geta = SimTrk_eta[gidx];
        gpz = gpt*sinh(geta);
	gpx = gpt*cos(gphi);
	gpy = gpt*sin(gphi);

	double tdrift, xdrift, ydrift, Rdrift;
	tdrift = (GV.ZCUT - sz)/(gpz - sz);
        xdrift = sx + tdrift*(gpx - sx);
        ydrift = sy + tdrift*(gpy - sy);
        Rdrift = sqrt( xdrift*xdrift + ydrift*ydrift );
	
	return Rdrift;

}
std::pair<std::vector<double>, std::vector<double> > getGEndpoints(recosim& s, int gidx){
	//for simtrack photon at gidx origin and endpoints first=(sx,sy,sz) second=(ex,ey,ez)
	std::vector<double> spoints(3);
	std::vector<double> epoints(3);
	//get start point
	auto nSimVtx = *(s.nSimVtx);

	auto& SimTrk_simvtx_Idx	 = s.SimTrk_simvtx_Idx;
	double sx,sy,sz;
	int svidx = SimTrk_simvtx_Idx[gidx];
	
	auto& SimVtx_x = s.SimVtx_x;
	auto& SimVtx_y = s.SimVtx_y;
	auto& SimVtx_z = s.SimVtx_z;
	sx = SimVtx_x[svidx];
	sy = SimVtx_y[svidx];
	sz = SimVtx_z[svidx];

	auto& SimTrk_trackId = s.SimTrk_trackId;
	int gtid = SimTrk_trackId[gidx];
	double ex,ey,ez;
	//have to loop simvtx and find simvtx that has this photon as parent
	ex=-9999;
	ey=-9999;
	ez=-9999;
	auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;
	for(int i=0; i<nSimVtx; i++){
		if( SimVtx_simtrk_parent_tid[i] == gtid ){
			ex = SimVtx_x[i]; 
			ey = SimVtx_y[i];
			ez = SimVtx_z[i];
		}		
	}
	spoints[0]=sx;
	spoints[1]=sy;
	spoints[2]=sz;
	epoints[0]=ex;
	epoints[1]=ey;
	epoints[2]=ez;
	std::pair<std::vector<double>, std::vector<double> > points{};
	points.first = spoints;
	points.second= epoints;
	return points;
	
}
///////////////////////////////////new photon function
/*double getPhotonRend(recosim&s, int gidx){
 
	auto& SimTrk_simvtx_Idx  = s.SimTrk_simvtx_Idx;
        double sx,sy,sz,ex,ey,ez;
        int svidx = SimTrk_simvtx_Idx[gidx];//use track id to get origin

        auto& SimVtx_x = s.SimVtx_x;
        auto& SimVtx_y = s.SimVtx_y;
        auto& SimVtx_z = s.SimVtx_z;
        sx = SimVtx_x[svidx];
        sy = SimVtx_y[svidx];
        sz = SimVtx_z[svidx];
	auto& SimTrk_trackId = s.SimTrk_trackId;
        int gtid = SimTrk_trackId[gidx];
	


}*/
std::vector<int> getSimpleGidxList( recosim& s){
//find all sim photons 

//each element is the idx of the SimTrk that has a pdg of 22 (that passes sim acceptance cuts)
    std::vector<int> gidxlist{};

	auto& SimTrk_pdgId = s.SimTrk_pdgId;
	int nSimTrk = (s.SimTrk_pdgId).GetSize();
	
	auto& SimTrk_pt = s.SimTrk_pt;
	auto& SimTrk_eta = s.SimTrk_eta;
	double gpt,geta,gpz,costg;

	auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
	auto& SimVtx_z = s.SimVtx_z;
	int svidx;
	int nSimVtx = (SimVtx_z).GetSize();
	double Sz;
	

	for(int i=0; i<nSimTrk; i++){
		if(SimTrk_pdgId[i] == 22){

			gpt = SimTrk_pt[i];
                        geta = SimTrk_eta[i];
                        gpz = gpt*sinh(geta);
                        costg = cos( atan2(gpt,gpz) );
			//its a photon but does it pass cuts (check kinematics, if it passes check origin)
			if(abs(costg) < GV.COSTCUT && gpt > GV.MINPT*2.){
				//we passed some cuts now check origin
				svidx = SimTrk_simvtx_Idx[i];
                                Sz = SimVtx_z[svidx]; 
								
				if(abs(Sz) < GV.ZCUT){
					//this photon is good, save it
					gidxlist.push_back(i);													

				}
			} 
		}
	}	
    return gidxlist;

}

std::vector<std::pair<int,int> > getGParentColl(recosim& s){
//first photon sim track idx
//second:
//0=prompt pi0 
//1=non prompt pi0
//2=electron/brem
//3=prompt other
//4=non prompt other
//5=non prompt no parent found
	std::vector< std::pair<int,int> > Gcoll{};
	
	int nSimVtx = (s.SimVtx_processType).GetSize();
        int nSimTrk = (s.SimTrk_simvtx_Idx).GetSize();

        auto& SimVtx_processType = s.SimVtx_processType;
        auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
        auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;
        auto& SimTrk_trackId = s.SimTrk_trackId;
	auto& SimTrk_pdgId = s.SimTrk_pdgId;
	auto& SimTrk_eta = s.SimTrk_eta;
	auto& SimTrk_phi = s.SimTrk_phi;
        auto& SimTrk_pt = s.SimTrk_pt;
	

	int sidx;
	int gidx;
	int gtype;
	double px1,py1,pz1,P1;
	double px2,py2,pz2,P2;
	double mggpi0;
	int ptid;
	for(int i=0; i<nSimTrk; i++){
		
		//is this sim trk 22?
		if(SimTrk_pdgId[i] != 22) continue;
		gtype=-1;
		//at this point we have a photon
		//if it isnt the last sim trk:
		sidx = SimTrk_simvtx_Idx[i];
		if( SimVtx_processType[ sidx ] == 0){ // this photon is prompt
			//check around the photon neighbors to see if we have a pi0
			px1=SimTrk_pt[i]*cos(SimTrk_phi[i]);
			py1=SimTrk_pt[i]*sin(SimTrk_phi[i]);
			pz1=SimTrk_pt[i]*sinh(SimTrk_eta[i]);
			P1=SimTrk_pt[i]*cosh(SimTrk_eta[i]);
			//extra cases, if first or last simtrk
			if(i==0){
				//check pi0 mass with next photon
				if( SimTrk_pdgId[i+1] == 22){
					px2=SimTrk_pt[i+1]*cos(SimTrk_phi[i+1]);
                                        py2=SimTrk_pt[i+1]*sin(SimTrk_phi[i+1]);
                                        pz2=SimTrk_pt[i+1]*sinh(SimTrk_eta[i+1]);
                                        P2=SimTrk_pt[i+1]*cosh(SimTrk_eta[i+1]);
					mggpi0 = 2*(P1*P2 - (px1*px2+py1*py2+pz1*pz2) );
                                        mggpi0 = sqrt(mggpi0);
					//if neighbor makes pi0 type 0 else type3
					if(mggpi0 > 0.130 && mggpi0 < 0.140){
						gtype = 0;
					}
					else{
						gtype = 3;
					}
				}
				else{ //no neighbor photon so type 3
					gtype = 3;
				}	
			}else if(i == (nSimTrk-1) ){//end edgecase start edgecase
				//check pi0 mass with previous photon
				if( SimTrk_pdgId[i-1] == 22){
                                        px2=SimTrk_pt[i-1]*cos(SimTrk_phi[i-1]);
                                        py2=SimTrk_pt[i-1]*sin(SimTrk_phi[i-1]);
                                        pz2=SimTrk_pt[i-1]*sinh(SimTrk_eta[i-1]);
                                        P2=SimTrk_pt[i-1]*cosh(SimTrk_eta[i-1]);
                                        mggpi0 = 2*(P1*P2 - (px1*px2+py1*py2+pz1*pz2) );
                                        mggpi0 = sqrt(mggpi0);
                                        //if neighbor makes pi0 type 0 else type3
                                        if(mggpi0 > 0.130 && mggpi0 < 0.140){
                                                gtype = 0;
                                        }
                                        else{
                                                gtype = 3;
                                        }
                                }
                                else{ //no neighbor photon so type 3
                                        gtype = 3;
                                }
				

			}else{//end edgecase start standard case
				//check pi0 mass with neighbor photon
				if( SimTrk_pdgId[i-1] == 22){
                                        px2=SimTrk_pt[i-1]*cos(SimTrk_phi[i-1]);
                                        py2=SimTrk_pt[i-1]*sin(SimTrk_phi[i-1]);
                                        pz2=SimTrk_pt[i-1]*sinh(SimTrk_eta[i-1]);
                                        P2=SimTrk_pt[i-1]*cosh(SimTrk_eta[i-1]);
                                        mggpi0 = 2*(P1*P2 - (px1*px2+py1*py2+pz1*pz2) );
                                        mggpi0 = sqrt(mggpi0);
                                        //if neighbor makes pi0 type 0 else type3
                                        if(mggpi0 > 0.130 && mggpi0 < 0.140){
                                                gtype = 0; 
					//	continue;
                                        }
                                       // else{
                                       //         gtype = 3;
                                       // }
                                }//end left neighbor check
				if( SimTrk_pdgId[i+1] == 22){
                                        px2=SimTrk_pt[i+1]*cos(SimTrk_phi[i+1]);
                                        py2=SimTrk_pt[i+1]*sin(SimTrk_phi[i+1]);
                                        pz2=SimTrk_pt[i+1]*sinh(SimTrk_eta[i+1]);
                                        P2=SimTrk_pt[i+1]*cosh(SimTrk_eta[i+1]);
                                        mggpi0 = 2*(P1*P2 - (px1*px2+py1*py2+pz1*pz2) );
                                        mggpi0 = sqrt(mggpi0);
                                        //if neighbor makes pi0 type 0 else type3
                                        if(mggpi0 > 0.130 && mggpi0 < 0.140){
                                                gtype = 0;
                                        }
                                       // else{
                                       //         gtype = 3;
                                       // }
                                }//end right neighbor check
				if(  gtype == -1) gtype = 3; //if no neigbors or no pi0 gtype = -1

			}//end normal case

		}//end ptype 0 check
		else{//if its not prompt then it has a parent, is it a pi0, electron or something else?
			//loop over simtrk again
			ptid = SimVtx_simtrk_parent_tid[sidx];
			bool found=false;
			for(int j=0; j<nSimTrk; j++){
				if(SimTrk_trackId[j] == ptid){
					if(abs(SimTrk_pdgId[j]) == 111){//pi0?
						gtype = 1;
					}
					else if(abs(SimTrk_pdgId[j]) == 11){//electron?	
						gtype = 2;
					}else{
						gtype = 4;
					}
					found=true;
					break;
				}//end ptid match found
				if(!found) gtype=5;
			}//end simtrk ptid matching loop
		}//end non prompt parent checking			
		
		//create a pair, push the idx and type
		std::pair<int, int> photon;
		photon.first =i;
		photon.second = gtype;
		Gcoll.push_back(photon);		
		
		
	}//end sim trk loop	
	
	return Gcoll;

}
//std::vector<std::pair<int,double> > getPCMatchingColl(recosim& s, double cutdL){ //try to match all reco pc
//std::map<int, std::pair<int,double> > getPCMatchingColl(recosim& s, double cutdL){
std::map<int, std::vector<double> > getPCMatchingColl(recosim& s, double cutdL){
//mindL is the minimum to be matched
//loop svtx , label all reco convs
//0= matched , nearest dR
//1= NI , nearest dR
//2= unmatched, nearest dR
//
///
//the vector elements are (match flag, dL, simR .... sim parent track id)
	//std::map< int, std::pair<int,double> > pcmap{};
	std::map<int, std::vector<double> > pcmap{};
	auto& SimVtx_processType = s.SimVtx_processType;
	int nSimVtx = (s.SimVtx_processType).GetSize();
	auto& PC_x = s.Conv_vtx_X;
	auto& PC_y = s.Conv_vtx_Y;
	auto& PC_z = s.Conv_vtx_Z;
	int npc = PC_x.GetSize();
	auto& SimVtx_x = s.SimVtx_x;
    	auto& SimVtx_y = s.SimVtx_y;
    	auto& SimVtx_z = s.SimVtx_z;		

        auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;

		

		double mindL,dL;
		double matchType;
		double sx,sy,sz;
		double cx,cy,cz;
		double sR;
		int tid;
	for(int i=0; i<npc; i++){
		mindL = 999999;
		matchType= -1;
		tid=-1;		
		cx = PC_x[i];
                cy = PC_y[i];
                cz = PC_z[i];
		for(int j=0; j<nSimVtx; j++){
			sx = SimVtx_x[j];
                        sy = SimVtx_y[j];
                        sz = SimVtx_z[j];
						
			dL = sqrt( (sx-cx)*(sx-cx) + (sy-cy)*(sy-cy) + (sz-cz)*(sz-cz) );
			if( dL < mindL ){
				mindL = dL;
				if(SimVtx_processType[j] == 14 && mindL < cutdL){
					 matchType=0;
					 sR = sqrt(sx*sx + sy*sy);
					 tid = SimVtx_simtrk_parent_tid[j];
				}
				else if(SimVtx_processType[j] != 14 && mindL < cutdL){
					 matchType=1;
					sR = sqrt(sx*sx + sy*sy);
					tid = SimVtx_simtrk_parent_tid[j];
				}else{
					matchType =2;
					sR=-1.;
					tid=-1;

				}
				
			}// end dL check
		}//end simvtx loop
		//std::pair<int,double> pcmatch{};
		std::vector<double> pcmatch(4);
		//pcmatch.first = matchType;
		//pcmatch.second = mindL;
		pcmatch.at(0) = (double) matchType;
		pcmatch.at(1) = mindL;
		pcmatch.at(2) = sR;
		pcmatch.at(3) = tid;
		pcmap[i]=pcmatch;
	}//end npc loop

	return pcmap;

}

////////////////////////////////////////////////////////////////////////////////////////
//(0) create cutmask
//
std::vector<bool> GetCutMask(recosim& s, std::vector<CommonVars> cv ){

  auto& PC_x = s.Conv_vtx_X;
  auto& PC_vTrack0_nBefore = s.Conv_nHitsBeforeVtx_Tk0;
  auto& PC_vTrack1_nBefore = s.Conv_nHitsBeforeVtx_Tk1;

  std::vector<bool> cutmask(PC_x.GetSize());
  double rerr,z,theta,fitprob,nBefore0,nBefore1;
  double PT1,PT2;
  for(int i=0; i<cutmask.size(); i++){
        rerr = cv[i].rerr;
        z = cv[i].z;
        theta = cv[i].theta;
        fitprob = cv[i].pfit;
        nBefore0 = PC_vTrack0_nBefore[i];
        nBefore1 = PC_vTrack1_nBefore[i];
	PT1 = std::sqrt(cv[i].px0p*cv[i].px0p + cv[i].py0p*cv[i].py0p);
        PT2 = std::sqrt(cv[i].px1p*cv[i].px1p + cv[i].py1p*cv[i].py1p);
//	std::cout<< rerr <<" "<< abs(z)<<" "<<abs(cos(theta))<<" "<<fitprob<<" "<<std::max(nBefore0,nBefore1)<<std::endl;
        if( rerr < GV.RERRCUT && abs(z) < GV.ZCUT && abs(cos(theta)) < GV.COSTCUT && fitprob > GV.FITPROBCUT && std::max(nBefore0,nBefore1)==0 && PT1 > GV.MINPT && PT2 > GV.MINPT){
                cutmask[i] = true;
        }
        else{
                cutmask[i] = false;
        }
   }
  return cutmask;      

}





