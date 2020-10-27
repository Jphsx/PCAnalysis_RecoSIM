//       #include "histset2enums.h" 
// bookeeping enumeration: if we do this we don't need to worry about hist pointer copies and merging
       enum th1d_ids{
		id_ptHist,
	 	id_pzHist,
		id_numpcHist,
		id_numSPCHist,
		id_ptSPCHist,
                id_numHGNPCHist,
            numTH1Hist};
       
	enum th2d_ids{
		id_xyHist,
		id_xywideHist,
		id_rphiHist, 
            numTH2Hist};
