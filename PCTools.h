//TODO make tools
//
//(1) hungarian filtering




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

//typedef TTreeReaderArray<Int_t> TTRAI;
//typedef int ID;
//Return struct of sim mask, currently masks process 14 and finds associated e/p pair and returns indices
//sim_pc* GetSimPC(TTRAI SimVtx_processType, TTRAI SimTrk_simvtx_Idx, TTRAI SimVtx_simtrk_parent_tid, TTRAI SimTrk_trackId  ){
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
