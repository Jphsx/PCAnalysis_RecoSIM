#include <iomanip>
#include "Hungarian.h"
#include "TMath.h"

//(0) create cutmask
//
std::vector<bool> GetCutMask(recosim& s){

  auto& PC_x = s.Conv_vtx_X;
  std::vector<bool> cutmask(PC_x.GetSize());
  for(int i=0; i<cutmask.size(); i++){
		cutmask[i] = true;
   }
  return cutmask;	
// rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT
 //              && fitprob > FITPROBCUT && std::max(nBefore0,nBefore1)==0 
}
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

    const bool lpr = true;       // print flag
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
