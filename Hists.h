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
    
   // TH1Manager.at(id_geomCorrHist) = new MyTH1D("geomCorrHist","r with fgeom;r cm; x/x_{0}",
    //TH1Manager.at(id_rHist) = new MyTH1D("rHist","r;r cm; nconv",
    
    TH1Manager.at(id_g0pt) = new MyTH1D("g0pt", "pType 0 Photon p_{T};p_{T} GeV;N_{#gamma}",100,0.0,5.0);
    TH1Manager.at(id_g0eta) = new MyTH1D("g0eta", "pType 0 Photon #eta;#eta GeV;N_{#gamma}",100,-2.4,2.4);
    TH1Manager.at(id_g0P) = new MyTH1D("g0P","pType 0 Photon P;P [GeV];N_{#gamma}",100,0.0,10.0);

    TH1Manager.at(id_g201pt) = new MyTH1D("g201pt", "pType 201 Photon p_{T};p_{T} GeV;N_{#gamma}",100,0.0,5.0);
    TH1Manager.at(id_g201eta) = new MyTH1D("g201eta", "pType 201 Photon #eta;#eta GeV;N_{#gamma}",100,-2.4,2.4);
    TH1Manager.at(id_g201P) = new MyTH1D("g201P","pType 201 Photon P;P [GeV];N_{#gamma}",100,0.0,10.0);


    TH1Manager.at(id_g3pt) = new MyTH1D("g3pt", "pType 3 Photon p_{T};p_{T} GeV;N_{#gamma}",100,0.0,5.0);
    TH1Manager.at(id_g3eta) = new MyTH1D("g3eta", "pType 3 Photon #eta;#eta GeV;N_{#gamma}",100,-2.4,2.4);
    TH1Manager.at(id_g3P) = new MyTH1D("g3P","pType 3 Photon P;P [GeV];N_{#gamma}",100,0.0,10.0);

    TH1Manager.at(id_gpt) = new MyTH1D("gpt", "Inclusive pType Photon p_{T};p_{T} GeV;N_{#gamma}",100,0.0,5.0);
    TH1Manager.at(id_geta) = new MyTH1D("geta", "Inclusive pType Photon #eta;#eta GeV;N_{#gamma}",100,-2.4,2.4);
    TH1Manager.at(id_gP) = new MyTH1D("gP","Inclusive pType Photon P;P [GeV];N_{#gamma}",100,0.0,10.0);


    //Rgeom 1, we fix theta to be fully transverse
    TH1Manager.at(id_Rgeom1_g0) = new MyTH1D("geom1_g0", "pType 0 Radial Flux with #theta^{*} = #phi * = #pi/2;R [cm];N_{#gamma}",100,0.5,100.5);
    TH1Manager.at(id_Rgeom1_g201) = new MyTH1D("geom1_g201", "pType 201 Radial Flux with #theta^{*} = #phi * = #pi/2;R [cm];N_{#gamma}",100,0.5,100.5);
    TH1Manager.at(id_Rgeom1_g3) = new MyTH1D("geom1_g3", "pType 3 Radial Flux with #theta^{*} = #phi * = #pi/2 ;R [cm];N_{#gamma}",100,0.5,100.5);
    TH1Manager.at(id_Rgeom1_g) = new MyTH1D("geom1_g", "Inclusive pType Radial Flux with #theta^{*} = #phi * = #pi/2;R [cm];N_{#gamma}",100,0.5,100.5);

    //thGeom1, we fix R to 2 and scan theta from 0.5 to 2.5 (plotted in eta)
    TH1Manager.at(id_thgeom1_g0) = new MyTH1D("thgeom1_g0", "pType 0 Angular Flux with R=50;#eta *;N_{#gamma}",500,-3,3);
    TH1Manager.at(id_thgeom1_g201) = new MyTH1D("thgeom1_g201", "pType 201 Angular Flux with R=50;#eta *;N_{#gamma}",500,-3,3);
    TH1Manager.at(id_thgeom1_g3) = new MyTH1D("thgeom1_g3", "pType 3 Angular Flux with R=50;#eta *;N_{#gamma}",500,-3,3);
    TH1Manager.at(id_thgeom1_g) = new MyTH1D("thgeom1_g", "Inclusive pType Angular Flux with R=50;#eta *;N_{#gamma}",500,-3,3);

    TH1Manager.at(id_thgeom1a_g0) = new MyTH1D("thgeom1a_g0", "pType 0 Angular Flux with R=50;#cos#theta *;N_{#gamma}",200,-1,1);
    TH1Manager.at(id_thgeom1a_g201) = new MyTH1D("thgeom1a_g201", "pType 201 Angular Flux with R=50;#cos#theta *;N_{#gamma}",200,-1,1);
    TH1Manager.at(id_thgeom1a_g3) = new MyTH1D("thgeom1a_g3", "pType 3 Angular Flux with R=50;#cos#theta *;N_{#gamma}",200,-1,1);
    TH1Manager.at(id_thgeom1a_g) = new MyTH1D("thgeom1a_g", "Inclusive pType Angular Flux with R=50;#cos#theta *;N_{#gamma}",200,-1,1);

    
// init TH2D
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
}//end histogram init

