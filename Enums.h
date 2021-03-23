//       #include "histset2enums.h" 
// bookeeping enumeration: if we do this we don't need to worry about hist pointer copies and merging
       enum th1d_ids{
		id_ptHist,
	 	id_pzHist,
		id_numpcHist,
		id_numpccutHist,
		id_numSPCHist,
		id_ptSPCHist,
                id_numHGNPCHist,
		id_effPtN,
		id_effPtD,
		id_purityPtN,
		id_purityPtD,
		id_effRN,
		id_effRD,
		id_purityRN,
		id_purityRD,
		id_effXPN,
		id_effXPD,
		id_purityXPN,
		id_purityXPD,
		id_trueGeom,
		id_ngPrompt,
		id_trueConv,
		id_Ng_BP1,
		id_Nc_BP1,
		id_Ng_BP2,
		id_Nc_BP2,
		id_nconvPt,
		id_nconvR,
		id_nconvXP,
		id_nconvPt_fake,
		id_nconvR_fake,
		id_nconvXP_fake,

            numTH1Hist};
       
	enum th2d_ids{
		id_xyHist,
		id_xywideHist,
		id_rphiHist, 
            numTH2Hist};
