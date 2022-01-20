//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 22 12:36:23 2020 by ROOT version 6.14/09
// from TTree Events/Events
// found on file: ../OutputFiles/SingleMuon2017_numEvent20000.root
//////////////////////////////////////////////////////////

#ifndef recosim_h
#define recosim_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
//using namespace std;

class recosim : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nConv = {fReader, "nConv"};
   TTreeReaderArray<Float_t> Conv_EoverP = {fReader, "Conv_EoverP"};
   TTreeReaderArray<Float_t> Conv_EoverPrefittedTracks = {fReader, "Conv_EoverPrefittedTracks"};
   TTreeReaderArray<Float_t> Conv_Tk0_chi2 = {fReader, "Conv_Tk0_chi2"};
   TTreeReaderArray<Float_t> Conv_Tk0_eta = {fReader, "Conv_Tk0_eta"};
   TTreeReaderArray<Float_t> Conv_Tk0_ndof = {fReader, "Conv_Tk0_ndof"};
   TTreeReaderArray<Float_t> Conv_Tk0_normalizedChi2 = {fReader, "Conv_Tk0_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_Tk0_phi = {fReader, "Conv_Tk0_phi"};
   TTreeReaderArray<Float_t> Conv_Tk0_pt = {fReader, "Conv_Tk0_pt"};
   TTreeReaderArray<Float_t> Conv_Tk1_chi2 = {fReader, "Conv_Tk1_chi2"};
   TTreeReaderArray<Float_t> Conv_Tk1_eta = {fReader, "Conv_Tk1_eta"};
   TTreeReaderArray<Float_t> Conv_Tk1_ndof = {fReader, "Conv_Tk1_ndof"};
   TTreeReaderArray<Float_t> Conv_Tk1_normalizedChi2 = {fReader, "Conv_Tk1_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_Tk1_phi = {fReader, "Conv_Tk1_phi"};
   TTreeReaderArray<Float_t> Conv_Tk1_pt = {fReader, "Conv_Tk1_pt"};
   TTreeReaderArray<Float_t> Conv_dPhiTracksAtVtx = {fReader, "Conv_dPhiTracksAtVtx"};
   TTreeReaderArray<Float_t> Conv_distOfMinimumApproach = {fReader, "Conv_distOfMinimumApproach"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_Tk0 = {fReader, "Conv_dlClosestHitToVtx_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_Tk1 = {fReader, "Conv_dlClosestHitToVtx_Tk1"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_err_Tk0 = {fReader, "Conv_dlClosestHitToVtx_err_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_err_Tk1 = {fReader, "Conv_dlClosestHitToVtx_err_Tk1"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_sig_Tk0 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_sig_Tk1 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk1"};
   TTreeReaderArray<Float_t> Conv_dxy = {fReader, "Conv_dxy"};
   TTreeReaderArray<Float_t> Conv_dz = {fReader, "Conv_dz"};
   TTreeReaderArray<Float_t> Conv_lxy = {fReader, "Conv_lxy"};
   TTreeReaderArray<Float_t> Conv_lz = {fReader, "Conv_lz"};
   TTreeReaderArray<Float_t> Conv_pairCotThetaSeparation = {fReader, "Conv_pairCotThetaSeparation"};
   TTreeReaderArray<Float_t> Conv_pairInvariantMass = {fReader, "Conv_pairInvariantMass"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Px = {fReader, "Conv_pairMomentum_Px"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Py = {fReader, "Conv_pairMomentum_Py"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Pz = {fReader, "Conv_pairMomentum_Pz"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_E = {fReader, "Conv_refittedPair4Momentum_E"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_M = {fReader, "Conv_refittedPair4Momentum_M"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Px = {fReader, "Conv_refittedPair4Momentum_Px"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Py = {fReader, "Conv_refittedPair4Momentum_Py"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Pz = {fReader, "Conv_refittedPair4Momentum_Pz"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Px = {fReader, "Conv_refittedPairMomentum_Px"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Py = {fReader, "Conv_refittedPairMomentum_Py"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Pz = {fReader, "Conv_refittedPairMomentum_Pz"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_X_Tk0 = {fReader, "Conv_tracksInnerPosition_X_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_X_Tk1 = {fReader, "Conv_tracksInnerPosition_X_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Y_Tk0 = {fReader, "Conv_tracksInnerPosition_Y_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Y_Tk1 = {fReader, "Conv_tracksInnerPosition_Y_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Z_Tk0 = {fReader, "Conv_tracksInnerPosition_Z_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Z_Tk1 = {fReader, "Conv_tracksInnerPosition_Z_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Px_Tk0 = {fReader, "Conv_tracksPin_Px_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Px_Tk1 = {fReader, "Conv_tracksPin_Px_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Py_Tk0 = {fReader, "Conv_tracksPin_Py_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Py_Tk1 = {fReader, "Conv_tracksPin_Py_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Pz_Tk0 = {fReader, "Conv_tracksPin_Pz_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Pz_Tk1 = {fReader, "Conv_tracksPin_Pz_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Px_Tk0 = {fReader, "Conv_tracksPout_Px_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Px_Tk1 = {fReader, "Conv_tracksPout_Px_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Py_Tk0 = {fReader, "Conv_tracksPout_Py_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Py_Tk1 = {fReader, "Conv_tracksPout_Py_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Pz_Tk0 = {fReader, "Conv_tracksPout_Pz_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Pz_Tk1 = {fReader, "Conv_tracksPout_Pz_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksSigned_d0_Tk0 = {fReader, "Conv_tracksSigned_d0_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksSigned_d0_Tk1 = {fReader, "Conv_tracksSigned_d0_Tk1"};
   TTreeReaderArray<Float_t> Conv_vtx_X = {fReader, "Conv_vtx_X"};
   TTreeReaderArray<Float_t> Conv_vtx_Y = {fReader, "Conv_vtx_Y"};
   TTreeReaderArray<Float_t> Conv_vtx_Z = {fReader, "Conv_vtx_Z"};
   TTreeReaderArray<Float_t> Conv_vtx_chi2 = {fReader, "Conv_vtx_chi2"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_00 = {fReader, "Conv_vtx_cov_00"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_01 = {fReader, "Conv_vtx_cov_01"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_02 = {fReader, "Conv_vtx_cov_02"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_11 = {fReader, "Conv_vtx_cov_11"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_12 = {fReader, "Conv_vtx_cov_12"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_22 = {fReader, "Conv_vtx_cov_22"};
   TTreeReaderArray<Float_t> Conv_vtx_ndof = {fReader, "Conv_vtx_ndof"};
   TTreeReaderArray<Float_t> Conv_vtx_normalizedChi2 = {fReader, "Conv_vtx_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_zOfPrimaryVertexFromTracks = {fReader, "Conv_zOfPrimaryVertexFromTracks"};

// New variables
   TTreeReaderArray<Float_t> Conv_Tk0_dPtRel = {fReader, "Conv_Tk0_dPtRel"};
   TTreeReaderArray<Float_t> Conv_Tk0_dR = {fReader, "Conv_Tk0_dR"};
   TTreeReaderArray<Float_t> Conv_Tk1_dPtRel = {fReader, "Conv_Tk1_dPtRel"};
   TTreeReaderArray<Float_t> Conv_Tk1_dR = {fReader, "Conv_Tk1_dR"};
   TTreeReaderArray<Float_t> Conv_vtxdl = {fReader, "Conv_vtxdl"};

   TTreeReaderArray<Int_t> Conv_Tk0_algo = {fReader, "Conv_Tk0_algo"};
   TTreeReaderArray<Int_t> Conv_Tk0_charge = {fReader, "Conv_Tk0_charge"};
   TTreeReaderArray<Int_t> Conv_Tk0_found = {fReader, "Conv_Tk0_found"};
   TTreeReaderArray<Int_t> Conv_Tk0_lost = {fReader, "Conv_Tk0_lost"};
   TTreeReaderArray<Int_t> Conv_Tk0_quality = {fReader, "Conv_Tk0_quality"};
   TTreeReaderArray<Int_t> Conv_Tk1_algo = {fReader, "Conv_Tk1_algo"};
   TTreeReaderArray<Int_t> Conv_Tk1_charge = {fReader, "Conv_Tk1_charge"};
   TTreeReaderArray<Int_t> Conv_Tk1_found = {fReader, "Conv_Tk1_found"};
   TTreeReaderArray<Int_t> Conv_Tk1_lost = {fReader, "Conv_Tk1_lost"};
   TTreeReaderArray<Int_t> Conv_Tk1_quality = {fReader, "Conv_Tk1_quality"};
   TTreeReaderArray<Int_t> Conv_algo = {fReader, "Conv_algo"};
   TTreeReaderArray<Int_t> Conv_nHitsBeforeVtx_Tk0 = {fReader, "Conv_nHitsBeforeVtx_Tk0"};
   TTreeReaderArray<Int_t> Conv_nHitsBeforeVtx_Tk1 = {fReader, "Conv_nHitsBeforeVtx_Tk1"};
   TTreeReaderArray<Int_t> Conv_nSharedHits = {fReader, "Conv_nSharedHits"};
   TTreeReaderArray<Int_t> Conv_nTracks = {fReader, "Conv_nTracks"};
// New
   TTreeReaderArray<Int_t> Conv_Tk0_Idx = {fReader, "Conv_Tk0_Idx"};
   TTreeReaderArray<Int_t> Conv_Tk1_Idx = {fReader, "Conv_Tk1_Idx"};
   TTreeReaderArray<Int_t> Conv_convVtxIdx = {fReader, "Conv_convVtxIdx"};

   TTreeReaderArray<Bool_t> Conv_isConverted = {fReader, "Conv_isConverted"};
   TTreeReaderValue<UInt_t> nPV = {fReader, "nPV"};
   TTreeReaderArray<Float_t> PV_X = {fReader, "PV_X"};
   TTreeReaderArray<Float_t> PV_Y = {fReader, "PV_Y"};
   TTreeReaderArray<Float_t> PV_Z = {fReader, "PV_Z"};
   TTreeReaderArray<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderArray<Float_t> PV_cov_00 = {fReader, "PV_cov_00"};
   TTreeReaderArray<Float_t> PV_cov_01 = {fReader, "PV_cov_01"};
   TTreeReaderArray<Float_t> PV_cov_02 = {fReader, "PV_cov_02"};
   TTreeReaderArray<Float_t> PV_cov_11 = {fReader, "PV_cov_11"};
   TTreeReaderArray<Float_t> PV_cov_12 = {fReader, "PV_cov_12"};
   TTreeReaderArray<Float_t> PV_cov_22 = {fReader, "PV_cov_22"};
   TTreeReaderArray<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderArray<Float_t> PV_normalizedChi2 = {fReader, "PV_normalizedChi2"};

// New
   TTreeReaderValue<UInt_t> nSimTrk = {fReader, "nSimTrk"};
   TTreeReaderArray<Float_t> SimTrk_charge = {fReader, "SimTrk_charge"};
   TTreeReaderArray<Float_t> SimTrk_eta = {fReader, "SimTrk_eta"};
   TTreeReaderArray<Float_t> SimTrk_phi = {fReader, "SimTrk_phi"};
   TTreeReaderArray<Float_t> SimTrk_pt = {fReader, "SimTrk_pt"};
   TTreeReaderArray<Int_t> SimTrk_pdgId = {fReader, "SimTrk_pdgId"};
   TTreeReaderArray<Int_t> SimTrk_simvtx_Idx = {fReader, "SimTrk_simvtx_Idx"};
   TTreeReaderArray<Int_t> SimTrk_trackId = {fReader, "SimTrk_trackId"};
   TTreeReaderValue<UInt_t> nSimVtx = {fReader, "nSimVtx"};
   TTreeReaderArray<Float_t> SimVtx_tof = {fReader, "SimVtx_tof"};
   TTreeReaderArray<Float_t> SimVtx_x = {fReader, "SimVtx_x"};
   TTreeReaderArray<Float_t> SimVtx_y = {fReader, "SimVtx_y"};
   TTreeReaderArray<Float_t> SimVtx_z = {fReader, "SimVtx_z"};
   TTreeReaderArray<Int_t> SimVtx_processType = {fReader, "SimVtx_processType"};
   TTreeReaderArray<Int_t> SimVtx_simtrk_parent_tid = {fReader, "SimVtx_simtrk_parent_tid"};
   
   //TTreeReaderValue<double> w_pvndof0_BC = {fReader, "w_pvndof0_BC"};
   TTreeReaderValue<double> w_pvtk_BC = {fReader, "w_pvtk_BC"};  
 
   recosim(TTree * /*tree*/ =0) { }
   virtual ~recosim() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

// VERY IMPORTANT TO COMMENT THIS OUT WHEN USED COMPILED ...
//   ClassDef(convsel,0);

};

#endif

#ifdef recosim_cxx
void recosim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t recosim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef recosim_cxx
