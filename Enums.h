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
		id_g0pt,
		id_g0P,
		id_g0eta,
		id_g201pt,
		id_g201P,
		id_g201eta,
		id_g3pt,
		id_g3P,
		id_g3eta,
		id_gpt,
		id_gP,
		id_geta,

		id_Rgeom1_g0,
		id_Rgeom1_g201,
		id_Rgeom1_g3,
		id_Rgeom1_g,

		id_thgeom1_g0,
		id_thgeom1_g201,
		id_thgeom1_g3,
		id_thgeom1_g,

		id_thgeom1a_g0,
		id_thgeom1a_g201,
		id_thgeom1a_g3,
		id_thgeom1a_g,

            numTH1Hist};
       
	enum th2d_ids{
		id_xyHist,
		id_xywideHist,
		id_rphiHist, 
            numTH2Hist};
