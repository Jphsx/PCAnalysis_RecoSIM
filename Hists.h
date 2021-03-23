// Now in "histset2init.h"

void histset::init(){
//init TH1D
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHist) = new MyTH1D("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_numpcHist) = new MyTH1D("numpcHist", "Number of PC;;Entries per bin", 20,-0.5, 19.5);

    TH1Manager.at(id_numpccutHist) = new MyTH1D("numpccutHist","Number of PC After Selection Cuts;;Entries per bin",20,-0.5,19.5);

    TH1Manager.at(id_numSPCHist) = new MyTH1D("numspcHist", "Number of Sim PC;;Entries per bin",100,-0.5,99.5);
    TH1Manager.at(id_ptSPCHist) = new MyTH1D("ptSPCHist", "Sim. #gamma Conv. p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}",100,0.0,5.0);

    TH1Manager.at(id_numHGNPCHist) = new MyTH1D("numHGNPCHist", "Number of PC after disambiguation;;Entries per bin",100,-0.5,99.5);


     TH1Manager.at(id_effPtN ) = new MyTH1D("effPtN","Efficiency;p_{T} [GeV];#varepsilon", 6, 0.0, 6.0);
     TH1Manager.at(id_effPtD ) = new MyTH1D("effPtD","Efficiency;p_{T} [GeV];#varepsilon", 6, 0.0, 6.0);
     TH1Manager.at(id_purityPtN ) = new MyTH1D("purityPtN","Purity;p_{T} [GeV];Purity", 6, 0.0, 6.0);
     TH1Manager.at(id_purityPtD ) = new MyTH1D("purityPtD","Purity;p_{T} [GeV];Purity", 6, 0.0, 6.0);
     
     Float_t Rbins[] = { 1,5,9,13,18,20 };
     Int_t  binnum = 5;
 
     TH1Manager.at(id_effRN ) = new MyTH1D("effRN","Efficiency;R [cm];#varepsilon",binnum,Rbins);
     TH1Manager.at(id_effRD ) = new MyTH1D("effRD","Efficiency;R [cm];#varepsilon",binnum,Rbins);
     TH1Manager.at(id_purityRN ) = new MyTH1D("purityRN","Purity;R [cm];Purity",binnum,Rbins);
     TH1Manager.at(id_purityRD ) = new MyTH1D("purityRD","Purity;R [cm];Purity",binnum,Rbins);
     TH1Manager.at(id_effXPN ) = new MyTH1D("effXPN","Efficiency;Positron p_{T} fraction X_{+};#varepsilon",10,0.,1.);
     TH1Manager.at(id_effXPD ) = new MyTH1D("effXPD","Efficiency;Positron p_{T} fraction X_{+};#varepsilon",10,0.,1.);
     TH1Manager.at(id_purityXPN ) = new MyTH1D("purityXPN","Purity;Positron p_{T} fraction X_{+};Purity",10,0.,1.);
     TH1Manager.at(id_purityXPD ) = new MyTH1D("purityXPD","Purity;Positron p_{T} fraction X_{+};Purity",10,0.,1.);

    //true geom counts
    TH1Manager.at(id_trueGeom) = new MyTH1D("trueGeom","N true photons",40,0.0,20.0);
    TH1Manager.at(id_ngPrompt) = new MyTH1D("ngPrompt","N true prompt",1,0.0,1.0);
    TH1Manager.at(id_trueConv) = new MyTH1D("trueConv","N true convs", 40,0.0,20.0);

    //study of BPIX1
    TH1Manager.at(id_Ng_BP1 ) = new MyTH1D("Ng_BP1", "N prompt photons",1,0.0,1.0);
    TH1Manager.at(id_Nc_BP1 ) = new MyTH1D("Nc_BP1", "N conv. BPIX1",1,0.0,1.0);
    //study of BPIX2
    TH1Manager.at(id_Ng_BP2 ) = new MyTH1D("Ng_BP2", "N prompt photons",1,0.0,1.0);
    TH1Manager.at(id_Nc_BP2 ) = new MyTH1D("Nc_BP2", "N conv. BPIX2",1,0.0,1.0);

    TH1Manager.at(id_nconvPt ) = new MyTH1D("nconvPt","N Conversions; p_{T} (GeV); N_{conv}",6,0.0,6.0);
    TH1Manager.at(id_nconvR ) = new MyTH1D("nconvR","N Conversions; R (cm); N_{conv}",40,0.0,20.0);
    TH1Manager.at(id_nconvXP ) = new MyTH1D("nconvXP","N Conversion; Positron p_{T} fraction X_{+}; N_{conv}",10,0.,1.);

    TH1Manager.at(id_nconvPt_fake) = new MyTH1D("nconvPt_fake", "N Fake Conversions;p_{T} (GeV); N_{conv}",6,0.0,6.0);
    TH1Manager.at(id_nconvR_fake) = new MyTH1D("nconvR_fake","N Fake Conversions;R (cm);N_{conv}",40,0.0,20.0);
    TH1Manager.at(id_nconvXP_fake) = new MyTH1D("nconvXP_fake","N Fake Conversions;Positron p_{T} fraction X_{+}",10,0.,1.);

	
	
// init TH2D
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
}//end histogram init

