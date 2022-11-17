// Now in "histset2init.h"

void histset::init(){
//init TH1D
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};NPC", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHist) = new MyTH1D("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_numpcHist) = new MyTH1D("numpcHist", "Number of PC;;Entries per bin", 20,-0.5, 19.5);
    TH1Manager.at(id_etaHist) = new MyTH1D("etaHist","PC #eta",60,-3,3);

    TH1Manager.at(id_numpccutHist) = new MyTH1D("numpccutHist","Number of PC After Selection Cuts;;Entries per bin",20,-0.5,19.5);
    TH1Manager.at(id_ptCutHist) = new MyTH1D("ptCutHist", "p_{T} Distribution after Selection;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzCutHist) = new MyTH1D("pzCutHist", "p_{Z} Distribution after Selection;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_numHGNPCHist) = new MyTH1D("numHGNPCHist", "Number of PC after disambiguation;;Entries per bin",100,-0.5,99.5);
    TH1Manager.at(id_ptHCutHist) = new MyTH1D("ptHCutHist", "p_{T} Distribution after Selection + HGN;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHCutHist) = new MyTH1D("pzHCutHist", "p_{Z} Distribution after Selection + HGN;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_etaCutHist) = new MyTH1D("etaCutHist", "#eta after selection",60,-3,3);
    TH1Manager.at(id_etaHGNHist) = new MyTH1D("etaHGNHist", "#eta after selection +HGN",60,-3,3);

    TH1Manager.at(id_r25RHist) = new MyTH1D("r25RHist","r dist (Raw)",250,0,25);
    TH1Manager.at(id_r25CHist) = new MyTH1D("r25CHist","r dist (Cut)",250,0,25);
    TH1Manager.at(id_r25HHist) = new MyTH1D("r25Hist","r dist (HGN)",250,0,25);
  
    TH1Manager.at(id_rho25RHist) = new MyTH1D("rho25RHist","rho dist (Raw)",250,0,25);
    TH1Manager.at(id_rho25CHist) = new MyTH1D("rho25CHist","rho dist (Cut)",250,0,25);
    TH1Manager.at(id_rho25HHist) = new MyTH1D("rho25Hist","rho dist (HGN)",250,0,25);

    TH1Manager.at(id_rps25RHist) = new MyTH1D("rps25RHist","rps dist (Raw)",250,0,25);
    TH1Manager.at(id_rps25CHist) = new MyTH1D("rps25CHist","rps dist (Cut)",250,0,25);
    TH1Manager.at(id_rps25HHist) = new MyTH1D("rps25Hist","rps dist (HGN)",250,0,25);

    TH1Manager.at(id_PVndof) = new MyTH1D("PVndof","Primary Vertex n.d.o.f, NPV=1",101,-0.5,100.5);
    TH1Manager.at(id_npv) = new MyTH1D("npv","N Primary Vertex",11,-0.5,10.5);
    TH1Manager.at(id_Rtheta) = new MyTH1D("Rtheta","raw theta",50,-2*3.14159,2*3.14159);
    TH1Manager.at(id_Rphi) = new MyTH1D("Rphi","raw phi",50,-2*3.14159,2*3.14159);
    TH1Manager.at(id_Htheta) = new MyTH1D("Htheta","HGN phi",50,-2*3.14159,2*3.14159);
    TH1Manager.at(id_Hphi) = new MyTH1D("Hphi","HGN theta",50,-2*3.14159,2*3.14159);
  
    TH1Manager.at(id_sumndof) = new MyTH1D("sumndof","sum of all pv ndofs",100,0.5,100.5);//start at 1 0 bin is big
    TH1Manager.at(id_ndofall) = new MyTH1D("ndofall","all ndofs, not npv cut",100,0.5,100.5);// also cut at 4 in code to compare to MC
    TH1Manager.at(id_ndofpv1) = new MyTH1D("ndofpv1","ndof for npv=1",100,0.5,100.5);
    TH1Manager.at(id_pvz) = new MyTH1D("pvz","z of pv's",40,-20,20);
    TH1Manager.at(id_pvr) = new MyTH1D("pvr","r of pv's",30,0,0.3);
    TH1Manager.at(id_pvz_cut) = new MyTH1D("pvzcut","z of pv's after tksum cut", 40, -20,20);
    TH1Manager.at(id_pvr_cut) = new MyTH1D("pvrcut","r of pv's after tksum cut", 30, 0,0.3);
    TH1Manager.at(id_pvdz) = new MyTH1D("pvdz","dz from pv0 and secondary pv",100,-10,10);
    TH1Manager.at(id_pvtrk) = new MyTH1D("pvtrk", "sum of track weights from ndof, all PV", 50,0.5,50.5);
    TH1Manager.at(id_pvtrkpv1) = new MyTH1D("pvtrkpv1", "sum of track weights from ndof, PV=1",50,0.5,50.5);
  
    TH1Manager.at(id_sumtksum_rlo) = new MyTH1D("sumtksum_rlo","tk sum with r<=bpix2",50,0.5,50.5);
    TH1Manager.at(id_sumtksum_rhi) = new MyTH1D("sumtksum_rhi","tk sum with r>bpix2",50,0.5,50.5);
 

     TH1Manager.at(id_rps25_b2p5) = new MyTH1D("rps25_b2p5","r ps b2p5",100,0,25);
     TH1Manager.at(id_rho25_b2p5) = new MyTH1D("rho25_b2p5","rho bp b2p5",100,0,25);
     TH1Manager.at(id_r25_b2p5) = new MyTH1D("r25_b2p5","r b2p5",100,0,25);
	const Int_t NBINS = 5;
        Double_t edges3[NBINS + 1] = {1.0, 5.0, 9.0, 13.5, 18.0, 25.0};
     TH1Manager.at(id_r25HHist_5bin) = new MyTH1D("r25hist_5bin","r custombin",NBINS,edges3);

    const Int_t ptNBINS = 5;
    Double_t ptedges[NBINS + 1] = {0.4, 1.0, 1.5, 2., 3., 5.};
    TH2Manager.at(id_rpt) = new MyTH2D("rpt","custom bin r and pt",NBINS,edges3,ptNBINS,ptedges);
 
// init TH2D
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
    TH2Manager.at(id_xyHistHGN) = new MyTH2D("xyHistHGN", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.); 
    TH2Manager.at(id_xywideHGN) = new MyTH2D("xywideHGN", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);

    TH2Manager.at(id_ndof_pcReta) = new MyTH2D("ndof_pcReta","ndof and pc raw eta",51,-0.5,50.5,60,-3,3);
    TH2Manager.at(id_ndof_pcRpt) = new MyTH2D("ndof_pcRpt","ndof and pc raw pt",51,-0.5,50.5,50,0,20); 
    TH2Manager.at(id_ndof_pcHeta) = new MyTH2D("ndof_pcHeta","ndof and hgn pc eta",51,-0.5,50.5,60,-3,3);
    TH2Manager.at(id_ndof_pcHpt) = new MyTH2D("ndof_pcHpt","ndof and hgn pc pt", 51,-0.5,50.5,50,0,20);
}//end histogram init

