#include <iomanip>
#include "Hungarian.h"
#include "TMath.h"
#include <set>
#include "TVector3.h"

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
//photon flux helper function

struct AngleLimits{
	double R;
	double A;
	double thetaStar;
	double etaStar;
	double etalow;
	double etaup;
	double dphi;

};

AngleLimits GetAngleLimits(double R, double A, double thetaStar){
	//return (eta,phi)
	//we compute theta,phi but convert theta to eta
	//this function assumes dtheta=dphi and the area segment is "square"
	AngleLimits al;
	al.R = R;
	al.A = A;
	al.thetaStar = thetaStar;
	double pi = 4*atan2(1,1);
	//double dpsi = sqrt( (3*A)/(4*R*sin(thetaStar)) );
	double dpsi = (3*A)/(4*pi*R*sin(thetaStar));	
	al.dphi=dpsi;
	double psi_dpsi_up = thetaStar+dpsi;
	double psi_dpsi_low = thetaStar-dpsi;
	al.etalow = -log( tan((thetaStar+dpsi)/2) );
	al.etaup =  -log( tan((thetaStar-dpsi)/2) );
	al.etaStar = -log( tan((thetaStar)/2));
	return al;
}

//photon flux method #2
class Point{
	public:
	
	double x;
	double y;
	double z;
	double r;
	double theta;
	double eta;
	double phi;
	Point(double xp, double yp, double zp){
		x=xp;
		y=yp;
		z=zp;
		PointRTP(xp,yp,zp);
	}
	void PointRTP(double xp, double yp, double zp){
		r = sqrt(xp*xp + yp*yp + zp*zp);
		phi = atan2(yp,xp);	
		double pi = 4*atan2(1,1);
		if (phi < 0) { phi += 2. * pi; }
		theta = atan2( sqrt(xp*xp + yp*yp) , zp);
		eta = -log(tan(theta/2.));
		 
	}
	void PrintPoint(){
		std::cout<<x<<" "<<y<<" "<<z<<" "<<std::endl;
	}
};
class CubeFace{
	public:
	Point* p0;//points that create face
	Point* p1;
	Point* p2;
	Point* p3;
	//constant for equation of plane
	double d;
	TVector3* v1;//vector 1 parallel to face
	TVector3* v2;//vector 2 parallel to face
	TVector3* n;//normal vector;
	CubeFace(Point* P0, Point* P1, Point* P2, Point* P3){
		p0=P0;
		p1=P1;
		p2=P2;
		p3=P3;
		//perform vector ops
		TVector3 _v1(p2->x - p0->x, p2->y - p0->y, p2->z - p0->z);
		TVector3 _v2(p1->x - p0->x, p1->y - p0->y, p1->z - p0->z);
		TVector3 _n = _v1.Cross(_v2);
		v1 = new TVector3(_v1);
		v2 = new TVector3(_v2);
		n = new TVector3(_n);
		d = n->x()*p0->x + n->y()*p0->y + n->z()*p0->z;
	
	}
	void PrintFace(){
		std::cout<<"Points: "<<std::endl;
		p0->PrintPoint();
		p1->PrintPoint();
		p2->PrintPoint();
		p3->PrintPoint();
		std::cout<<"Vecs "<<std::endl;
		v1->Print();
		v2->Print();
		n->Print();
	}
};
class Cube{
	public:
	std::vector<CubeFace*> faces;
	Point* corner; //corner is easiest to calculate	
	double L;//length of a face
	
	Point* center;
	//Point* centerRTP;

	Cube(double l, Point* p){
		L=l;
		corner = p;
		center = new Point( p->x + L/2., p->y + L/2. , p->z + L/2. );
		//centerRTP = new Point(0,0,0);
		//centerRTP->PointRTP(center->x, center->y, center->z);
		//centerRTP = new Point(center->x, center->y, center->z);
	}
	void PrintCube(){
		for(int i=0; i<faces.size(); i++){
			std::cout<<"Face: "<<i<<std::endl;
			faces[i]->PrintFace();
		}
	}
	
};

Cube* buildCube(double x0, double y0, double z0, double L){

	Point* p1 = new Point(x0,y0,z0);
	Cube* c = new Cube(L, p1);	

	//find all 8 points of the cube
	
	Point* p2 = new Point(x0,y0,z0+L);
	Point* p3 = new Point(x0,y0+L,z0);
	Point* p4 = new Point(x0,y0+L,z0+L);
	Point* p5 = new Point(x0+L, y0, z0);
	Point* p6 = new Point(x0+L,y0,z0+L);
	Point* p7 = new Point(x0+L, y0+L, z0);
	Point* p8 = new Point(x0+L, y0+L, z0+L);

	//define faces
	c->faces.push_back( new CubeFace(p1,p2,p3,p4) );
	c->faces.push_back( new CubeFace(p2,p6,p4,p8) );
	c->faces.push_back( new CubeFace(p6,p5,p8,p7) );
	c->faces.push_back( new CubeFace(p5,p1,p7,p3) );
	c->faces.push_back( new CubeFace(p1,p2,p6,p5) );
	c->faces.push_back( new CubeFace(p3,p4,p8,p7) );

	return c;
}
bool hasCubeIntersection(Cube* c, double gx, double gy, double gz, double px, double py, double pz ){
	

	double t;
	double ndp, ndr0; //n dot p and n dot r0 (r0 is photon vertex)
	double nx,ny,nz,d;
	double xint,yint,zint;
	
	//determine cube bounds
	double xlo,xup,ylo,yup,zlo,zup;
	xlo = c->faces[0]->p1->x;
	ylo = c->faces[0]->p1->y;
	zlo = c->faces[0]->p1->z;
	xup = xlo + c->L;
	yup = ylo + c->L;
	zup = zlo + c->L;
	
	for(int i=0; i < c->faces.size(); i++){

		
		nx = c->faces[i]->n->x();
		ny = c->faces[i]->n->y();
		nz = c->faces[i]->n->z();	
		d = c->faces[i]->d;
		ndp = nx*px + ny*py + nz*pz;
		if(ndp==0) continue; //no possible intersection go to next face
	
		ndr0 = nx*gx + ny*gy + nz*gz;
	
		t = (d-ndr0)/ndp; 
		xint = gx + px*t;
		yint = gy + py*t;
		zint = gz + pz*t;
	
		if( xlo <= xint && xup >= xint){
		if( ylo <= yint && yup >= yint){
		if( zlo <= zint && zup >= zint){
			//we have box intersection
			return true;
		}}}//end box conditions	
	}

	return false;

}


struct sim_g{
	std::vector<int> g_idx;//index on simtrack of this photon
	std::vector<int> parent_simvtx_idx; 
	std::vector<int> parent_simvtx_ptype;
	std::vector<int> parent_pdg;
	std::vector<int> endpoint_vtx_idx;
	
};
sim_g GetSimG(recosim& s){
	sim_g SG;

	std::vector<int> g_idx;
	std::vector<int> parent_simvtx_idx;
	std::vector<int> parent_simvtx_ptype;
	std::vector<int> parent_pdg;
	std::vector<int> endpoint_vtx_idx;
	
	int nSimTrk = (s.SimTrk_simvtx_Idx).GetSize();
	auto& SimTrk_pdg = s.SimTrk_pdgId;
	
	auto& SimVtx_processType = s.SimVtx_processType;
	auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
	auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;
    	auto& SimTrk_trackId = s.SimTrk_trackId;
	auto& SimTrk_pdgId = s.SimTrk_pdgId;

	int ptype;
	int pvtxidx;
	int ptid;
	int gtid;
	int gidx;

	//parent/child tracking variables
	int pidx;
	int cidx;
	int cvtxidx;
	int ppdg;

	//loop over sim tracks, is this pdg == 22?
	for(int i=0; i<nSimTrk; i++){
		if( SimTrk_pdg[i] == 22){ // found a photon
			
			gidx = i;
			pvtxidx = SimTrk_simvtx_Idx[i];
			gtid= SimTrk_trackId[i];
			ptid= SimVtx_simtrk_parent_tid[pvtxidx];
			ptype = SimVtx_processType[pvtxidx];
			cidx = -1;
			pidx = -1;
			cvtxidx = -1;
			ppdg=-999;
			//reloop over simtracks get parent pdg, and find child/endpoint at the same time
			for( int j=0; j<nSimTrk; j++){
				//find child if any
				if( gtid == SimVtx_simtrk_parent_tid[ SimTrk_simvtx_Idx[j] ] ){
					cvtxidx = SimTrk_simvtx_Idx[j]; //this is the index to the photon endpoint (of the sim vtx)
				}						
				//find the parent and store its pdg
				if( ptid == SimTrk_trackId[j] ){
					ppdg = SimTrk_pdgId[j];
				} 

			}//end j simtrk loop	
					
			

			

			g_idx.push_back( gidx );
			//map out photons ancestry
			parent_simvtx_idx.push_back( pvtxidx );
			parent_simvtx_ptype.push_back( ptype );
			parent_pdg.push_back( ppdg );
			endpoint_vtx_idx.push_back( cvtxidx );
			
			
		}//end 22 check
	}//end simtrkloop

	//populate struct and return
	SG.g_idx = g_idx;
	SG.parent_simvtx_idx = parent_simvtx_idx;
	SG.parent_simvtx_ptype = parent_simvtx_ptype;
	SG.parent_pdg = parent_pdg;
	SG.endpoint_vtx_idx = endpoint_vtx_idx;

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


////////////////////////////////////////////////////////////////////////////////////////
//(0) create cutmask
//
std::vector<bool> GetCutMask(recosim& s, std::vector<CommonVars> cv ){

  auto& PC_x = s.Conv_vtx_X;
  auto& PC_vTrack0_nBefore = s.Conv_nHitsBeforeVtx_Tk0;
  auto& PC_vTrack1_nBefore = s.Conv_nHitsBeforeVtx_Tk1;

  std::vector<bool> cutmask(PC_x.GetSize());
  double rerr,z,theta,fitprob,nBefore0,nBefore1;
  for(int i=0; i<cutmask.size(); i++){
        rerr = cv[i].rerr;
        z = cv[i].z;
        theta = cv[i].theta;
        fitprob = cv[i].pfit;
        nBefore0 = PC_vTrack0_nBefore[i];
        nBefore1 = PC_vTrack1_nBefore[i];
//	std::cout<< rerr <<" "<< abs(z)<<" "<<abs(cos(theta))<<" "<<fitprob<<" "<<std::max(nBefore0,nBefore1)<<std::endl;
        if( rerr < GV.RERRCUT && abs(z) < GV.ZCUT && abs(cos(theta)) < GV.COSTCUT && fitprob > GV.FITPROBCUT && std::max(nBefore0,nBefore1)==0){
                cutmask[i] = true;
        }
        else{
                cutmask[i] = false;
        }
   }
  return cutmask;      

}





