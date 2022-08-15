//    #include "mylocaltree.h"     //All the variable incantations needed
   	
	//always make a local copy, if its a value dereference. 
    //if you dont do this scope/dereferencing will get really weird, clunky, and unmanageable
	//have to auto& or myreader will try to register copy of the readerarray pointer

// Here we retain the old variable names as far as possible to minimize required changes

//	auto w_pvndof0_BC = *(s.w_pvndof0_BC);	
	auto w_pvtk_BC = *(s.w_pvtk_BC);

	auto& PC_x = s.Conv_vtx_X;
	auto& PC_y = s.Conv_vtx_Y;
	auto& PC_z = s.Conv_vtx_Z;
	auto& PC_vtx_chi2 = s.Conv_vtx_chi2;
	auto& PC_vtx_ndof = s.Conv_vtx_ndof;

	auto& PC_vtx_sigmaxx = s.Conv_vtx_cov_00;
	auto& PC_vtx_sigmaxy = s.Conv_vtx_cov_01;
	auto& PC_vtx_sigmayy = s.Conv_vtx_cov_11;
	auto& PC_vtx_sigmazz = s.Conv_vtx_cov_22;

	auto numberOfPC = *(s.nConv);
	auto nPV = *(s.nPV);
	auto& PV_X = s.PV_X;
	auto& PV_Y = s.PV_Y;
	auto& PV_Z = s.PV_Z;
	auto& PV_ndof = s.PV_ndof;
	auto& PV_chi2 = s.PV_chi2;

	auto& PC_vTrack0_pt = s.Conv_Tk0_pt;
	auto& PC_vTrack0_phi = s.Conv_Tk0_phi;
	auto& PC_vTrack0_eta = s.Conv_Tk0_eta;
	auto& PC_vTrack1_pt = s.Conv_Tk1_pt;
	auto& PC_vTrack1_phi = s.Conv_Tk1_phi;
	auto& PC_vTrack1_eta = s.Conv_Tk1_eta;

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

// Quality
    auto& Tk0_chi2 = s.Conv_Tk0_chi2;
    auto& Tk1_chi2 = s.Conv_Tk1_chi2;
    auto& Tk0_ndof = s.Conv_Tk0_ndof;
    auto& Tk1_ndof = s.Conv_Tk1_ndof;
    auto& Tk0_sd0  = s.Conv_tracksSigned_d0_Tk0;
    auto& Tk1_sd0  = s.Conv_tracksSigned_d0_Tk1;

// New variables
    auto& PC_mpair = s.Conv_pairInvariantMass;
    auto& PC_dcottheta = s.Conv_pairCotThetaSeparation;
    auto& PC_dmin = s.Conv_distOfMinimumApproach;
    auto& PC_dphi = s.Conv_dPhiTracksAtVtx;
    auto& PC_zPV = s.Conv_zOfPrimaryVertexFromTracks;

//    auto isRealData = *(s.isRealData);
//    bool isRealData = true;

    auto lumiSection = *(s.luminosityBlock);
//    auto& mcpu = s.MC_PUInfo_numberOfInteractions; 
    auto runNumber = *(s.run);
    auto eventNumber = *(s.event);
//    auto nMCPU = *(s.numberOfMC_PUInfo);
    int nMCPU = 0;

    auto& PC_E = s.Conv_refittedPair4Momentum_E;
    auto& PC_M = s.Conv_refittedPair4Momentum_M;
    auto& PC_Px = s.Conv_refittedPair4Momentum_Px;
    auto& PC_Py = s.Conv_refittedPair4Momentum_Py;
    auto& PC_Pz = s.Conv_refittedPair4Momentum_Pz;

//Also needed for asymmetry
    auto& PC_vTrack0_charge = s.Conv_Tk0_charge;
    auto& PC_vTrack1_charge = s.Conv_Tk1_charge;

    auto& PC_vTrack0_quality = s.Conv_Tk0_quality;
    auto& PC_vTrack1_quality = s.Conv_Tk1_quality;

    auto& PC_vTrack0_found = s.Conv_Tk0_found;
    auto& PC_vTrack1_found = s.Conv_Tk1_found;

    auto& PC_vTrack0_lost = s.Conv_Tk0_lost;
    auto& PC_vTrack1_lost = s.Conv_Tk1_lost;

    auto& Conv_vtx_normalizedChi2 = s.Conv_vtx_normalizedChi2;
// New variables - April 2020.
    auto& PC_vTrack0_nBefore = s.Conv_nHitsBeforeVtx_Tk0;
    auto& PC_vTrack1_nBefore = s.Conv_nHitsBeforeVtx_Tk1;
    auto& PC_nSharedHits = s.Conv_nSharedHits;

    auto& SimVtx_processType = s.SimVtx_processType;
    auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
    auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;
    auto& SimTrk_trackId = s.SimTrk_trackId;
    auto& SimTrk_pt = s.SimTrk_pt;
    auto& SimTrk_eta = s.SimTrk_eta;
    auto& SimTrk_phi = s.SimTrk_phi;
    auto& SimTrk_pdgId = s.SimTrk_pdgId;
    auto& SimVtx_x = s.SimVtx_x;
    auto& SimVtx_y = s.SimVtx_y;
    auto& SimVtx_z = s.SimVtx_z;
    auto nSimVtx = *(s.nSimVtx);
    auto nSimTrk = *(s.nSimTrk);
