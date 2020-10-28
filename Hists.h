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
	
// init TH2D
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
}//end histogram init

