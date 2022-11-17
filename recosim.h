//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 16 16:15:46 2021 by ROOT version 6.14/09
// from TTree tree/Photon Conversion Tree
// found on file: /home/t3-ku/janguian/jsingera/PC/MinBias3_2018B_UL/MinimumBias3/MinBias3_2018B_UL/210813_153458/0000/pc_100.root
//////////////////////////////////////////////////////////

#ifndef recosim_h
#define recosim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class recosim : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> run = {fReader, "run"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   //TTreeReaderValue<Int_t> lumiBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<Int_t> LN = {fReader, "LN"};
   TTreeReaderValue<UInt_t> evt_timeStamp = {fReader, "evt_timeStamp"};
   TTreeReaderValue<Int_t> orbit = {fReader, "orbit"};
   TTreeReaderValue<Int_t> bunchCrossing = {fReader, "bunchCrossing"};
   TTreeReaderValue<Int_t> LS = {fReader, "LS"};
   TTreeReaderValue<Int_t> LS_timeStamp = {fReader, "LS_timeStamp"};
   TTreeReaderValue<Double_t> ilumi_del = {fReader, "ilumi_del"};
   TTreeReaderValue<Double_t> ilumi_rec = {fReader, "ilumi_rec"};
   TTreeReaderValue<Int_t> nConv = {fReader, "nConv"};
   TTreeReaderValue<std::vector<bool>> Conv_isConverted = {fReader, "Conv_isConverted"};
   TTreeReaderArray<unsigned int> Conv_nTracks = {fReader, "Conv_nTracks"};
   TTreeReaderArray<double> Conv_pairInvariantMass = {fReader, "Conv_pairInvariantMass"};
   TTreeReaderArray<double> Conv_pairCotThetaSeparation = {fReader, "Conv_pairCotThetaSeparation"};
   TTreeReaderArray<double> Conv_pairMomentum_Px = {fReader, "Conv_pairMomentum_Px"};
   TTreeReaderArray<double> Conv_pairMomentum_Py = {fReader, "Conv_pairMomentum_Py"};
   TTreeReaderArray<double> Conv_pairMomentum_Pz = {fReader, "Conv_pairMomentum_Pz"};
   TTreeReaderArray<double> Conv_refittedPair4Momentum_Px = {fReader, "Conv_refittedPair4Momentum_Px"};
   TTreeReaderArray<double> Conv_refittedPair4Momentum_Py = {fReader, "Conv_refittedPair4Momentum_Py"};
   TTreeReaderArray<double> Conv_refittedPair4Momentum_Pz = {fReader, "Conv_refittedPair4Momentum_Pz"};
   TTreeReaderArray<double> Conv_refittedPair4Momentum_E = {fReader, "Conv_refittedPair4Momentum_E"};
   TTreeReaderArray<double> Conv_refittedPair4Momentum_M = {fReader, "Conv_refittedPair4Momentum_M"};
   TTreeReaderArray<double> Conv_refittedPairMomentum_Px = {fReader, "Conv_refittedPairMomentum_Px"};
   TTreeReaderArray<double> Conv_refittedPairMomentum_Py = {fReader, "Conv_refittedPairMomentum_Py"};
   TTreeReaderArray<double> Conv_refittedPairMomentum_Pz = {fReader, "Conv_refittedPairMomentum_Pz"};
   TTreeReaderArray<double> Conv_EoverP = {fReader, "Conv_EoverP"};
   TTreeReaderArray<double> Conv_EoverPrefittedTracks = {fReader, "Conv_EoverPrefittedTracks"};
   TTreeReaderArray<double> Conv_distOfMinimumApproach = {fReader, "Conv_distOfMinimumApproach"};
   TTreeReaderArray<double> Conv_dPhiTracksAtVtx = {fReader, "Conv_dPhiTracksAtVtx"};
   TTreeReaderArray<double> Conv_dxy = {fReader, "Conv_dxy"};
   TTreeReaderArray<double> Conv_dz = {fReader, "Conv_dz"};
   TTreeReaderArray<double> Conv_lxy = {fReader, "Conv_lxy"};
   TTreeReaderArray<double> Conv_lz = {fReader, "Conv_lz"};
   TTreeReaderArray<double> Conv_zOfPrimaryVertexFromTracks = {fReader, "Conv_zOfPrimaryVertexFromTracks"};
   TTreeReaderArray<double> Conv_tracksSigned_d0_Tk0 = {fReader, "Conv_tracksSigned_d0_Tk0"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_X_Tk0 = {fReader, "Conv_tracksInnerPosition_X_Tk0"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_Y_Tk0 = {fReader, "Conv_tracksInnerPosition_Y_Tk0"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_Z_Tk0 = {fReader, "Conv_tracksInnerPosition_Z_Tk0"};
   TTreeReaderArray<double> Conv_tracksPout_Px_Tk0 = {fReader, "Conv_tracksPout_Px_Tk0"};
   TTreeReaderArray<double> Conv_tracksPout_Py_Tk0 = {fReader, "Conv_tracksPout_Py_Tk0"};
   TTreeReaderArray<double> Conv_tracksPout_Pz_Tk0 = {fReader, "Conv_tracksPout_Pz_Tk0"};
   TTreeReaderArray<double> Conv_tracksPin_Px_Tk0 = {fReader, "Conv_tracksPin_Px_Tk0"};
   TTreeReaderArray<double> Conv_tracksPin_Py_Tk0 = {fReader, "Conv_tracksPin_Py_Tk0"};
   TTreeReaderArray<double> Conv_tracksPin_Pz_Tk0 = {fReader, "Conv_tracksPin_Pz_Tk0"};
   TTreeReaderArray<double> Conv_nHitsBeforeVtx_Tk0 = {fReader, "Conv_nHitsBeforeVtx_Tk0"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_Tk0 = {fReader, "Conv_dlClosestHitToVtx_Tk0"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_err_Tk0 = {fReader, "Conv_dlClosestHitToVtx_err_Tk0"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_sig_Tk0 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk0"};
   TTreeReaderArray<double> Conv_tracksSigned_d0_Tk1 = {fReader, "Conv_tracksSigned_d0_Tk1"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_X_Tk1 = {fReader, "Conv_tracksInnerPosition_X_Tk1"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_Y_Tk1 = {fReader, "Conv_tracksInnerPosition_Y_Tk1"};
   TTreeReaderArray<double> Conv_tracksInnerPosition_Z_Tk1 = {fReader, "Conv_tracksInnerPosition_Z_Tk1"};
   TTreeReaderArray<double> Conv_tracksPout_Px_Tk1 = {fReader, "Conv_tracksPout_Px_Tk1"};
   TTreeReaderArray<double> Conv_tracksPout_Py_Tk1 = {fReader, "Conv_tracksPout_Py_Tk1"};
   TTreeReaderArray<double> Conv_tracksPout_Pz_Tk1 = {fReader, "Conv_tracksPout_Pz_Tk1"};
   TTreeReaderArray<double> Conv_tracksPin_Px_Tk1 = {fReader, "Conv_tracksPin_Px_Tk1"};
   TTreeReaderArray<double> Conv_tracksPin_Py_Tk1 = {fReader, "Conv_tracksPin_Py_Tk1"};
   TTreeReaderArray<double> Conv_tracksPin_Pz_Tk1 = {fReader, "Conv_tracksPin_Pz_Tk1"};
   TTreeReaderArray<double> Conv_nHitsBeforeVtx_Tk1 = {fReader, "Conv_nHitsBeforeVtx_Tk1"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_Tk1 = {fReader, "Conv_dlClosestHitToVtx_Tk1"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_err_Tk1 = {fReader, "Conv_dlClosestHitToVtx_err_Tk1"};
   TTreeReaderArray<double> Conv_dlClosestHitToVtx_sig_Tk1 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk1"};
   TTreeReaderArray<double> Conv_nSharedHits = {fReader, "Conv_nSharedHits"};
   TTreeReaderArray<double> Conv_algo = {fReader, "Conv_algo"};
   TTreeReaderArray<double> Conv_vtx_X = {fReader, "Conv_vtx_X"};
   TTreeReaderArray<double> Conv_vtx_Y = {fReader, "Conv_vtx_Y"};
   TTreeReaderArray<double> Conv_vtx_Z = {fReader, "Conv_vtx_Z"};
   TTreeReaderArray<double> Conv_vtx_cov_00 = {fReader, "Conv_vtx_cov_00"};
   TTreeReaderArray<double> Conv_vtx_cov_01 = {fReader, "Conv_vtx_cov_01"};
   TTreeReaderArray<double> Conv_vtx_cov_02 = {fReader, "Conv_vtx_cov_02"};
   TTreeReaderArray<double> Conv_vtx_cov_11 = {fReader, "Conv_vtx_cov_11"};
   TTreeReaderArray<double> Conv_vtx_cov_12 = {fReader, "Conv_vtx_cov_12"};
   TTreeReaderArray<double> Conv_vtx_cov_22 = {fReader, "Conv_vtx_cov_22"};
   TTreeReaderArray<double> Conv_vtx_chi2 = {fReader, "Conv_vtx_chi2"};
   TTreeReaderArray<double> Conv_vtx_normalizedChi2 = {fReader, "Conv_vtx_normalizedChi2"};
   TTreeReaderArray<double> Conv_vtx_ndof = {fReader, "Conv_vtx_ndof"};
   TTreeReaderArray<double> Conv_Tk0_pt = {fReader, "Conv_Tk0_pt"};
   TTreeReaderArray<double> Conv_Tk0_eta = {fReader, "Conv_Tk0_eta"};
   TTreeReaderArray<double> Conv_Tk0_phi = {fReader, "Conv_Tk0_phi"};
   TTreeReaderArray<double> Conv_Tk0_charge = {fReader, "Conv_Tk0_charge"};
   TTreeReaderArray<double> Conv_Tk0_chi2 = {fReader, "Conv_Tk0_chi2"};
   TTreeReaderArray<double> Conv_Tk0_normalizedChi2 = {fReader, "Conv_Tk0_normalizedChi2"};
   TTreeReaderArray<double> Conv_Tk0_ndof = {fReader, "Conv_Tk0_ndof"};
   TTreeReaderArray<double> Conv_Tk0_found = {fReader, "Conv_Tk0_found"};
   TTreeReaderArray<double> Conv_Tk0_lost = {fReader, "Conv_Tk0_lost"};
   TTreeReaderArray<double> Conv_Tk0_quality = {fReader, "Conv_Tk0_quality"};
   TTreeReaderArray<double> Conv_Tk0_algo = {fReader, "Conv_Tk0_algo"};
   TTreeReaderArray<double> Conv_Tk1_pt = {fReader, "Conv_Tk1_pt"};
   TTreeReaderArray<double> Conv_Tk1_eta = {fReader, "Conv_Tk1_eta"};
   TTreeReaderArray<double> Conv_Tk1_phi = {fReader, "Conv_Tk1_phi"};
   TTreeReaderArray<double> Conv_Tk1_charge = {fReader, "Conv_Tk1_charge"};
   TTreeReaderArray<double> Conv_Tk1_chi2 = {fReader, "Conv_Tk1_chi2"};
   TTreeReaderArray<double> Conv_Tk1_normalizedChi2 = {fReader, "Conv_Tk1_normalizedChi2"};
   TTreeReaderArray<double> Conv_Tk1_ndof = {fReader, "Conv_Tk1_ndof"};
   TTreeReaderArray<double> Conv_Tk1_found = {fReader, "Conv_Tk1_found"};
   TTreeReaderArray<double> Conv_Tk1_lost = {fReader, "Conv_Tk1_lost"};
   TTreeReaderArray<double> Conv_Tk1_quality = {fReader, "Conv_Tk1_quality"};
   TTreeReaderArray<double> Conv_Tk1_algo = {fReader, "Conv_Tk1_algo"};
   TTreeReaderValue<Int_t> nPV = {fReader, "nPV"};
   TTreeReaderValue<Int_t> nPV_cut = {fReader, "nPV_cut"};
   TTreeReaderValue<std::vector<bool>> nPV_mask = {fReader, "nPV_mask"};
   TTreeReaderArray<double> PV_X = {fReader, "PV_X"};
   TTreeReaderArray<double> PV_Y = {fReader, "PV_Y"};
   TTreeReaderArray<double> PV_Z = {fReader, "PV_Z"};
   TTreeReaderArray<double> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderArray<double> PV_normalizedChi2 = {fReader, "PV_normalizedChi2"};
   TTreeReaderArray<double> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderArray<double> PV_cov_00 = {fReader, "PV_cov_00"};
   TTreeReaderArray<double> PV_cov_01 = {fReader, "PV_cov_01"};
   TTreeReaderArray<double> PV_cov_02 = {fReader, "PV_cov_02"};
   TTreeReaderArray<double> PV_cov_11 = {fReader, "PV_cov_11"};
   TTreeReaderArray<double> PV_cov_12 = {fReader, "PV_cov_12"};
   TTreeReaderArray<double> PV_cov_22 = {fReader, "PV_cov_22"};


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

   //ClassDef(recosim,0);

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
