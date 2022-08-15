// Now in "histset2init.h"

void histset::init(){
//init TH1D
     TH1::SetDefaultSumw2();
    int customNbins = 350;//changed from 100 (2p5 mm bin) for truth testing and to resolve BP
    int customEbins = 40;
    //int customNbins2 = 250;
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 20.0);
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

    TH1Manager.at(id_rerrHGNHist) = new MyTH1D("rerrHGNHist", "radial resolution", 30,0,0.3);

    TH1Manager.at(id_r25RHist) = new MyTH1D("r25RHist","r dist (Raw)",250,0,25);
    TH1Manager.at(id_r25CHist) = new MyTH1D("r25CHist","r dist (Cut)",250,0,25);
    TH1Manager.at(id_r25HHist) = new MyTH1D("r25Hist","r dist (HGN)",250,0,25);
    TH1Manager.at(id_r25Hist_b2p5) = new MyTH1D("r25Hist_b2p5","r dist(2.5mm)",customNbins,0,25);
    TH1Manager.at(id_r25Hist_b2p5_nowt) = new MyTH1D("r25Hist_b2p5_nowt","r dist(2.5mm) no wt",customNbins,0,25);

    
 //   TH1Manager.at(id_r25Hist_Truth_b2p5) = new MyTH1D("r25Hist_truth","r dist sim truth", customNbins,0,25);

    TH1Manager.at(id_PVndof) = new MyTH1D("PVndof","Primary Vertex n.d.o.f, NPV=1",101,-0.5,100.5);
    TH1Manager.at(id_npv) = new MyTH1D("npv","N Primary Vertex",11,-0.5,10.5);
	
    TH1Manager.at(id_pvz) = new MyTH1D("pvz","z of pv's",100,-10,10);
    TH1Manager.at(id_pvdz) = new MyTH1D("pvdz","dz from pv0 and secondary pv",100,-10,10);
    TH1Manager.at(id_pvtrksum) = new MyTH1D("pvtrksum", "sum of ( track weight sum from ndof), all PV",50,0.5,50.5);
    TH1Manager.at(id_sumpvtrksum) = new MyTH1D("sumpvtrksum","sum of (sum track weight), all PV",50,0.5,50.5);
   
    //exploring excess of pc in bpix 1/2
    TH1Manager.at(id_sumtksum_rlo) = new MyTH1D("sumtksum_rlo","sum of tk sum with r<=bpix2",50,0.5,50.5);
    TH1Manager.at(id_sumtksum_rhi) = new MyTH1D("sumtksum_rhi","sum of tk sum with r>bpix2",50,0.5,50.5);
 
        //efficiency hists
        const Int_t NBINS = 5;
   	Double_t edges[NBINS + 1] = {1.0, 5.0, 9.0, 13.5, 18.0, 25.0}; //custom radial bins 
    TH1Manager.at(id_eRden) = new MyTH1D("eRden","eff R denominator",NBINS,edges);
    TH1Manager.at(id_eRden_nowt) = new MyTH1D("eRden_nowt","eff R denominator unweighted",NBINS,edges);
    TH1Manager.at(id_eRden_b2p5) = new MyTH1D("eRden_b2p5","eff R denominator b2p5",customNbins,0,25);
    TH1Manager.at(id_eRden_b2p5_nowt) = new MyTH1D("eRden_b2p5_nowt","eff R denominator b2p5",customNbins,0,25);
   // TH1Manager.at(id_eRden_true_b2p5) = new MyTH1D("eRden_true_b2p5","eff R denom truth",customNbins,0,25);
    TH1Manager.at(id_ePtden) = new MyTH1D("ePtden","eff Pt denominator",10,0,5);
    TH1Manager.at(id_etksumden) = new MyTH1D("etksumden","eff ntk denominator",12,8,20);
    TH1Manager.at(id_debug1) = new MyTH1D("debug1","debug1",NBINS,edges);
    TH1Manager.at(id_debug2) = new MyTH1D("debug2","debug2",100,0,25);

    TH1Manager.at(id_r25coarse) = new MyTH1D("r25coarse","r dist (coarse)",NBINS,edges);
   
    TH1Manager.at(id_matchdR) = new MyTH1D("matchdR","dL dist of matches",20,0,2);
    
	Double_t edges2[NBINS + 1] = {1.0, 5.0, 9.0, 13.5, 18.0, 25.0};
    TH1Manager.at(id_eRnum) = new MyTH1D("eRnum","eff R numerator",NBINS,edges2);
    TH1Manager.at(id_eRnum_nowt) = new MyTH1D("eRnum_nowt","eff R numerator unweighted",NBINS,edges2);

    TH1Manager.at(id_eRnum_b2p5) = new MyTH1D("eRnum_b2p5","eff R numerator b2p5",customNbins,0,25);
      TH1Manager.at(id_eRnum_b2p5_nowt) = new MyTH1D("eRnum_b2p5_nowt","eff R numerator b2p5 unweighted",customNbins,0,25);

 //   TH1Manager.at(id_eRnum_true_b2p5) = new MyTH1D("eRnum_true_b2p5","eff R numerator truth", customNbins,0,25);
    TH1Manager.at(id_ePtnum) = new MyTH1D("ePtnum","eff Pt numerator",10,0,5);
    TH1Manager.at(id_etksumnum) = new MyTH1D("etksumnum","eff ntk numerator",12,8,20);

	//fakes
 	//const Int_t NBINS = 5;
	Double_t edges3[NBINS + 1] = {1.0, 5.0, 9.0, 13.5, 18.0, 25.0};
    TH1Manager.at(id_eRnumf) = new MyTH1D("eRnumf","eff R numerator",NBINS,edges3);
    TH1Manager.at(id_eRnumf_b2p5) = new MyTH1D("eRnumf_b2p5","eff R numerator b2p5",50,0,25);
    TH1Manager.at(id_ePtnumf) = new MyTH1D("ePtnumf","eff Pt numerator",10,0,5);
    TH1Manager.at(id_etksumnumf) = new MyTH1D("etksumnumf","eff ntk numerator",12,8,20);

	

	//flux plots
	//radial flux
	//TH1Manager.at(id_radflux) = new MyTH1D("radflux","photon r flux",25,0,25);
	TH1Manager.at(id_radflux) = new MyTH1D("radflux_b2p5","photon r flux b2p5",100,0,25);
	TH1Manager.at(id_radflux_nowt) = new MyTH1D("radflux_b2p5_nowt","photon r flux b2p5 nowt",100,0,25);
	TH1Manager.at(id_radfluxcoarse) = new MyTH1D("radfluxcoarse","photon r flux coarse", NBINS, edges);
	//TH1Manager.at(id_radflux) = new MyTH1D("radfluxp1","photon r flux",
	//flux composition
	TH1Manager.at(id_fluxcomp) = new MyTH1D("fluxcomp","all photon composition",6,-0.5,5.5);
	//TH1Manager.at(id_rfluxcomp) = new MyTH1D("rfluxcomp","photon r composition",6,-0.6,5.5);

	TH1Manager.at(id_gflux) = new MyTH1D("gflux", "photon r flux per mm",customNbins,0,25);	
	TH1Manager.at(id_gflux_eta) = new MyTH1D("gflux_eta","photon physics #eta",100,-5,5);
	TH1Manager.at(id_gflux_nowt) = new MyTH1D("gflux_nowt","photon r flux nowt",customNbins,0,25);
	//TH1Manager.at(id_gflux_b2p5) = new MyTH1D("
	TH1Manager.at(id_pc_chi2ndof) = new MyTH1D("pc_chi2ndof","chi2ndof",20,0,2);	
	//TH1Manager.at(id_ginc) = new MyTH1D("ginc", "reco photon incidence"
//	TH1Manager.at(id_truinc) = new MyTH1D("truinc", "true photon incidence (r-z dip)",30,0,3.14159);
		
	
/*
	TH1Manager.at(id_s1_pc) = new MyTH1D("s1_pc", "reconstructed conversions in eta segment 1;#eta;NPC",16,-0.8,0.8);
	TH1Manager.at(id_s2_pc) = new MyTH1D("s2_pc", "reconstructed conversions in eta segment 2;#eta;NPC",16,-0.8,0.8);
	TH1Manager.at(id_s3_pc) = new MyTH1D("s3_pc", "reconstructed conversions in eta segment 3;#eta;NPC",16,-0.8,0.8);
	TH1Manager.at(id_s4_pc) = new MyTH1D("s4_pc", "reconstructed conversions in eta segment 4;#eta;NPC",16,-0.8,0.8);
	TH1Manager.at(id_s5_pc) = new MyTH1D("s5_pc", "reconstructed conversions in eta segment 5;#eta;NPC",16,-0.8,0.8);
	TH1Manager.at(id_s6_pc) = new MyTH1D("s6_pc", "reconstructed conversions in eta segment 6;#eta;NPC",16,-0.8,0.8);

        TH1Manager.at(id_s1_effN) = new MyTH1D("s1_effN","efficiency numerator in eta segment 1;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s2_effN) = new MyTH1D("s2_effN","efficiency numerator in eta segment 2;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s3_effN) = new MyTH1D("s3_effN","efficiency numerator in eta segment 3;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s4_effN) = new MyTH1D("s4_effN","efficiency numerator in eta segment 4;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s5_effN) = new MyTH1D("s5_effN","efficiency numerator in eta segment 5;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s6_effN) = new MyTH1D("s6_effN","efficiency numerator in eta segment 6;#eta",16,-0.8,0.8);
	
	TH1Manager.at(id_s1_effD) = new MyTH1D("s1_effD","efficiency numerator in eta segment 1;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s2_effD) = new MyTH1D("s2_effD","efficiency numerator in eta segment 2;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s3_effD) = new MyTH1D("s3_effD","efficiency numerator in eta segment 3;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s4_effD) = new MyTH1D("s4_effD","efficiency numerator in eta segment 4;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s5_effD) = new MyTH1D("s5_effD","efficiency numerator in eta segment 5;#eta",16,-0.8,0.8);
        TH1Manager.at(id_s6_effD) = new MyTH1D("s6_effD","efficiency numerator in eta segment 6;#eta",16,-0.8,0.8);

	TH1Manager.at(id_s1_ng) = new MyTH1D("s1_ng", "number of photons in eta segment 1;#eta;N#gamma",16,-0.8,0.8);
        TH1Manager.at(id_s2_ng) = new MyTH1D("s2_ng", "number of photons in eta segment 2;#eta;N#gamma",16,-0.8,0.8);
        TH1Manager.at(id_s3_ng) = new MyTH1D("s3_ng", "number of photons in eta segment 3;#eta;N#gamma",16,-0.8,0.8);
        TH1Manager.at(id_s4_ng) = new MyTH1D("s4_ng", "number of photons in eta segment 4;#eta;N#gamma",16,-0.8,0.8);
        TH1Manager.at(id_s5_ng) = new MyTH1D("s5_ng", "number of photons in eta segment 5;#eta;N#gamma",16,-0.8,0.8);
        TH1Manager.at(id_s6_ng) = new MyTH1D("s6_ng", "number of photons in eta segment 6;#eta;N#gamma",16,-0.8,0.8);
*/
 
// init TH2D
//
    const Int_t ptNBINS = 5;
    Double_t ptedges[NBINS + 1] = {0.4, 1.0, 1.5, 2., 3., 5.};

    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
  
    //TH2Manager.at(id_reta_pc) = new MyTH2D("reta_pc", "reconstructed conversions r-eta;#eta;R",16,-0.8,0.8,100,0,25);
    //TH2Manager.at(id_reta_effN) = new MyTH2D("reta_effN", "eff numerator r-eta;#eta;R",16,-0.8,0.8,100,0,25);
    //TH2Manager.at(id_reta_effD) = new MyTH2D("reta_effD", "eff denominator r-eta;#eta;R",16,-0.8,0.8,100,0,25);
    //TH2Manager.at(id_reta_ng) = new MyTH2D("reta_ng","number of photons r-eta;#eta;R",16,-0.8,0.8,100,0,25);
    TH2Manager.at(id_reta_pc) = new MyTH2D("reta_pc", "reconstructed conversions r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);
    TH2Manager.at(id_reta_effN) = new MyTH2D("reta_effN", "eff numerator r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);
    TH2Manager.at(id_reta_effD) = new MyTH2D("reta_effD", "eff denominator r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);
        TH2Manager.at(id_reta_effN_nowt) = new MyTH2D("reta_effN_nowt", "no wteff numerator r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);
    TH2Manager.at(id_reta_effD_nowt) = new MyTH2D("reta_effD_nowt", "no wt eff denominator r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);

    TH2Manager.at(id_reta_ng) = new MyTH2D("reta_ng","number of photons r-eta;#eta;R",16,-0.8,0.8,NBINS,edges3);


 
    TH2Manager.at(id_ndof_pcReta) = new MyTH2D("ndof_pcReta","ndof and pc raw eta",51,-0.5,50.5,60,-3,3);
    TH2Manager.at(id_ndof_pcRpt) = new MyTH2D("ndof_pcRpt","ndof and pc raw pt",51,-0.5,50.5,50,0,20);
    TH2Manager.at(id_ndof_pcHeta) = new MyTH2D("ndof_pcHeta","ndof and hgn pc eta",51,-0.5,50.5,60,-3,3);
    TH2Manager.at(id_ndof_pcHpt) = new MyTH2D("ndof_pcHpt","ndof and hgn pc pt", 51,-0.5,50.5,50,0,20);

    //TH2Manager.at(id_effptr_num) = new MyTH2D("effptr_num", "eff by layer in pt", NBINS, edges3, 10,0,5);
   // TH2Manager.at(id_effptr_den) = new MyTH2D("effptr_den", "eff by layer in pt", NBINS, edges3, 10,0,5);
   // TH2Manager.at(id_effptr_num_nowt) = new MyTH2D("effptr_num_nowt", "eff by layer in pt no wt", NBINS, edges3, 10,0,5);
   // TH2Manager.at(id_effptr_den_nowt) = new MyTH2D("effptr_den_nowt", "eff by layer in pt no wt", NBINS, edges3, 10,0,5);
      TH2Manager.at(id_effptr_num) = new MyTH2D("effptr_num", "eff by layer in pt", NBINS, edges3, ptNBINS,ptedges);
      TH2Manager.at(id_effptr_den) = new MyTH2D("effptr_den", "eff by layer in pt", NBINS, edges3, ptNBINS,ptedges);
      TH2Manager.at(id_effptr_num_nowt) = new MyTH2D("effptr_num_nowt", "eff by layer in pt no wt", NBINS, edges3, ptNBINS,ptedges);
      TH2Manager.at(id_effptr_den_nowt) = new MyTH2D("effptr_den_nowt", "eff by layer in pt no wt", NBINS, edges3, ptNBINS,ptedges);



    TH2Manager.at(id_pt_leadfound) = new MyTH2D("ptleadfound", "leading found hits",11,-0.5,10.5, 10,0,5);
    TH2Manager.at(id_pt_subfound) = new MyTH2D("ptsubfound", "sub leading found hits", 11, -0.5,10.5, 10,0,5);
    TH2Manager.at(id_pt_leadlost) = new MyTH2D("ptleadlost", "leading lost hits", 11, -0.5, 10.5, 10,0,5);
    TH2Manager.at(id_pt_sublost) = new MyTH2D("ptsublost", "sub leading lost hits", 11, -0.5, 10.5, 10, 0,5);
    TH2Manager.at(id_pt_leadqual) = new MyTH2D("ptleadqual", "leading quality", 9, -1.5, 7.5, 10, 0,5);
    TH2Manager.at(id_pt_subqual) = new MyTH2D("ptsubqual", "sub leading quality", 9, -1.5, 7.5, 10, 0, 5);
    TH2Manager.at(id_pt_shared) = new MyTH2D("ptsharedhit", "shared hits",11,-0.5,10.5,10,0,5);
   
    TH2Manager.at(id_pt_rerr) = new MyTH2D("pt_rerr","pc radial error and pt;rerr;pt", 30,0,0.3,10,0,5);
    TH2Manager.at(id_r_rerr) = new MyTH2D("r_err","pc radial error and r;rerr",30,0,0.3,NBINS,edges3);

	//extract 1d from phi slice
    TH2Manager.at(id_r_phi_p6) = new MyTH2D("id_npc_phi","#pi/6 slices in R-#phi;R(cm);#phi(rad)",customNbins,0,25,12,0,2*PI);
    TH2Manager.at(id_effrphiN) = new MyTH2D("id_effrphiN","eff by layer in phi", NBINS, edges3, 12,0,2*PI);
    TH2Manager.at(id_effrphiD) = new MyTH2D("id_effrphiD","eff by layer in phi", NBINS, edges3, 12,0,2*PI);
    TH2Manager.at(id_effrphiN_nowt) = new MyTH2D("id_effrphiN_nowt","eff by layer in phi no wt", NBINS, edges3, 12,0,2*PI);
    TH2Manager.at(id_effrphiD_nowt) = new MyTH2D("id_effrphiD_nowt","eff by layer in phi no wt", NBINS, edges3, 12,0,2*PI);
		
		// id_pc_chi2ndof
    TH2Manager.at(id_nchi2_r) = new MyTH2D("nchi2_r","reco pc chi2/ndof and r",customNbins,0,25,20,0,2);
    TH2Manager.at(id_nchi2_dr) = new MyTH2D("nchi2_dr","reco pc chi2/ndof and radial error",30,0,0.3,20,0,2);

    TH2Manager.at(id_gfluxE) = new MyTH2D("gfluxE","flux per radial bin as a function of true g E;R (cm);E (GeV)", customNbins,0,25,customEbins,0.4,10);
    TH2Manager.at(id_effER_den) = new MyTH2D("effER_den","true conversion as a function of true E; R(cm);E (Gev)", customNbins,0,25,customEbins,0.4,10); 
}//end histogram init

