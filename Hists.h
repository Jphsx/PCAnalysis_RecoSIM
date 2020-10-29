// Now in "histset2init.h"

void histset::init(){
//init TH1D
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHist) = new MyTH1D("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_numpcHist) = new MyTH1D("numpcHist", "Number of PC;;Entries per bin", 20,-0.5, 19.5);

    TH1Manager.at(id_numpccutHist) = new MyTH1D("numpccutHist","Number of PC After Selection Cuts;;Entries per bin",20,-0.5,19.5);
    TH1Manager.at(id_ptCutHist) = new MyTH1D("ptCutHist", "p_{T} Distribution after Selection;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzCutHist) = new MyTH1D("pzCutHist", "p_{Z} Distribution after Selection;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);

    TH1Manager.at(id_numSPCHist) = new MyTH1D("numspcHist", "Number of Sim PC;;Entries per bin",100,-0.5,99.5);
    TH1Manager.at(id_ptSPCHist) = new MyTH1D("ptSPCHist", "Sim. #gamma Conv. p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}",100,0.0,5.0);

    TH1Manager.at(id_numHGNPCHist) = new MyTH1D("numHGNPCHist", "Number of PC after disambiguation;;Entries per bin",100,-0.5,99.5);
    TH1Manager.at(id_ptHCutHist) = new MyTH1D("ptHCutHist", "p_{T} Distribution after Selection + HGN;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHCutHist) = new MyTH1D("pzHCutHist", "p_{Z} Distribution after Selection + HGN;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);

    TH1Manager.at(id_xplusSPCHist) = new MyTH1D("xplusSPCHist", "Photon pT < 0 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_minTkPtSPCHist) = new MyTH1D("minTkPtSPCHist", "Minimum P_{T} of e/p; p_{T} [GeV]; Entries per bin",100,0.0,5.0);
    //TH1Manager.at(id_xplusSBGHist) = new MyTH1D("xplusSBGHist", "Positron pT Fraction in Sim. Bkg; Positron PT Fraction; Entries per 0.01 bin", 100,0.0,1.0);

	//the 0 and 1 appended mean 0=before any cuts , 1= after cuts and HGN
    TH1Manager.at(id_mindLc0) = new MyTH1D("mindLc0", "Distance from Reco. PC to nearest Sim. PC;dL [cm];Entries per bin", 100,0.0,5.0);
    TH1Manager.at(id_mindLb0) = new MyTH1D("mindLb0", "Distance from Reco. PC to nearest Sim. Bkg.; dL [cm]; Entries per bin",100,0.0,5.0);
    TH1Manager.at(id_mindLx0) = new MyTH1D("mindLx0", "Distance from Reco. PC to nearest Sim. Vtx.; dL [cm]; Entries per bin",100,0.0,5.0);
    TH1Manager.at(id_mindL2x0) = new MyTH1D("mindL2x0","Distance from Reco. PC to next nearest Sim. Vtx.; dL [cm]; Entries per bin", 100,0.0,5.0);
	
    TH1Manager.at(id_mindLc1) = new MyTH1D("mindLc1", "Distance from Reco. PC to nearest Sim. PC;dL [cm];Entries per bin", 100,0.0,5.0);
    TH1Manager.at(id_mindLb1) = new MyTH1D("mindLb1", "Distance from Reco. PC to nearest Sim. Bkg.; dL [cm]; Entries per bin",100,0.0,5.0);
    TH1Manager.at(id_mindLx1) = new MyTH1D("mindLx1", "Distance from Reco. PC to nearest Sim. Vtx.; dL [cm]; Entries per bin",100,0.0,5.0);
    TH1Manager.at(id_mindL2x1) = new MyTH1D("mindL2x1","Distance from Reco. PC to next nearest Sim. Vtx.; dL [cm]; Entries per bin", 100,0.0,5.0);


    TH1Manager.at(id_comp0) = new MyTH1D("comp0", "Reconstructed Composition; Vertex Type; Fractional Composition",3,-0.5,2.5);
    TH1Manager.at(id_comp1) =  new MyTH1D("comp1", "Reconstructed Composition; Vertex Type; Fractional Composition",3,-0.5,2.5);




    //numerators and denominators for eff 1 -- raw denominator(no longer true it may have a cutmask)
    TH1Manager.at(id_ptnum1) = new MyTH1D("ptnum1","P_{T} Reco. PC matched to Sim. PC;p_{T} [GeV];Entries per bin",50,0,5.0);
    TH1Manager.at(id_ptnumB1) = new MyTH1D("ptnumB1","P_{T} Reco. PC matched to Sim BG;p_{T} [GeV];Entries per bin",50,0,5.0);
    TH1Manager.at(id_ptnumU1) = new MyTH1D("ptnumU1","P_{T} Reco. PC Unmatched to Sim;p_{T} [GeV];Entries per bin",50,0,5.0);
    TH1Manager.at(id_xpnum1) = new MyTH1D("xpnum1","Reco. PC matched to Sim. PC Positron pT Fraction;Positron pT Fraction; Entries per bin",20,0.0,1.0);
    TH1Manager.at(id_xpnumB1) = new MyTH1D("xpnumB1","Reco. PC matched to Sim. BG;p_{T} [GeV];Entries per bin",20,0.0,1.0);
    TH1Manager.at(id_xpnumU1) = new MyTH1D("xpnumU1","Reco. PC Unmatched to Sim.;p_{T} [GeV];Entries per bin",20,0.0,1.0);
    TH1Manager.at(id_minTknum1) = new MyTH1D("minTknum1", "Reco. PC matched to Sim. PC pair min. p_{T};p_{T} [GeV];Entries per bin",40,0.0,2.0);
    TH1Manager.at(id_minTknumB1) = new MyTH1D("minTknumB1", "Reco. PC matched to Sim. BG pair min. p_{T};p_{T} [GeV];Entries per bin",40,0.0,2.0);
    TH1Manager.at(id_minTknumU1) = new MyTH1D("minTknumU1", "Reco. PC Unmatched to Sim.;p_{T} [GeV]; Entries per bin",40,0.0,2.0);
    TH1Manager.at(id_Rnum1) = new MyTH1D("Rnum1", "Reco. PC matched to Sim. PC Radius;R [cm];Entries per bin", 100,0.,10.);
    TH1Manager.at(id_RnumB1) = new MyTH1D("RnumB1","Reco. PC matched to Sim. BG Radius;R [cm];Entries per bin", 100,0.,10.);
    TH1Manager.at(id_RnumU1) = new MyTH1D("RnumU1","Reco. PC Unmatched to Sim.;R [cm];Entries per bin",100,0.,10.);
    TH1Manager.at(id_Rwidenum1) = new MyTH1D("Rwidenum1","Reco. PC matched to Sim. PC Radius; R [cm];Entries per bin",250,0.,25.);
    TH1Manager.at(id_RwidenumB1) = new MyTH1D("RwidenumB1","Reco. PC matched to Sim. BG Radius; R [cm];Entries per bin",250,0.,25.);
    TH1Manager.at(id_RwidenumU1) = new MyTH1D("RwidenumU1","Reco. PC Unmatched to Sim. Radius; R [cm];Entries per bin",250,0.,20.);
    
    TH1Manager.at(id_ptden1) = new MyTH1D("ptden1","Sim. PC;p_{T} [GeV];Entries per bin",50,0,5.0);
    TH1Manager.at(id_xpden1) = new MyTH1D("xpden1","Sim. PC Positron pT Fraction;Positron pT Fraction; Entries per bin",20,0.0,1.0);
    TH1Manager.at(id_minTkden1) = new MyTH1D("minTkden1", "Sim. PC pair min. p_{T};p_{T} [GeV];Entries per bin",40,0.0,2.0);
    TH1Manager.at(id_Rden1) = new MyTH1D("Rden1", "Sim. PC Radius;R [cm];Entries per bin",100,0.,10);
    TH1Manager.at(id_Rwideden1) = new MyTH1D("Rwideden1","Sim. PC Radius;R [cm]; Entries per bin",250,0.,25.);

// init TH2D
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
}//end histogram init

