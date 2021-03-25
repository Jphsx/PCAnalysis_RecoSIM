
#include "tdrstyle.C"


TCanvas* gCanv(std::string cname, TH1D* h, std::string opt){
	TCanvas* c = new TCanvas(cname.c_str(),cname.c_str());
	h->Draw(opt.c_str());
	return c;
}
void teff(std::string cname, TH1D* n , TH1D* d, std::string axis=""){
	TEfficiency* t = new TEfficiency(*n,*d);
	//if(axis != "");
	t->SetTitle(axis.c_str());
	TCanvas* ceffr = new TCanvas(cname.c_str(),cname.c_str());
	t->Draw();
}

void format3(TH1D* h1, TH1D* h2, TH1D* h3 ){
	h1->SetMarkerStyle(8);
	h1->SetMarkerColor(kBlack);
	h1->SetLineColor(kBlack);
	h1->SetLineWidth(2);

	h2->SetMarkerStyle(8);
	h2->SetMarkerColor(kBlue);
	h2->SetLineColor(kBlue);
	h2->SetLineWidth(2);

	h3->SetMarkerStyle(8);
	h3->SetMarkerColor(kRed);
	h3->SetLineColor(kRed);
	h3->SetLineWidth(2);
}
void legend3(TH1D* h1, TH1D* h2, TH1D* h3 ){

	auto legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->SetHeader("","C"); // option "C" allows to center the header
	legend->AddEntry(h1,"Data","lp");
	legend->AddEntry(h2,"MC (Incl.)","lp");
	legend->AddEntry(h3,"MC Fakes","lp");
	legend->Draw();
}
void legend2(TH1D* h1,TH1D* h2){
	 auto legend = new TLegend(0.1,0.7,0.48,0.9);
      //  legend->SetHeader("","C"); // option "C" allows to center the header
        legend->AddEntry(h1,"MC (Incl.)","lp");
        legend->AddEntry(h2,"MC Fakes","lp");
        legend->Draw();
}
void mcdataplot(){
	
	setTDRStyle();
	gROOT->ForceStyle();
	//TFile* fmc = TFile::Open("/home/justin/work/research/DPG/3-13-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/mc_weighted/Outfile.root");
//	TFile* fmc = TFile::Open("/home/justin/work/research/DPG/3-13-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/local_test/Outfile.root");
	TFile* fmc = TFile::Open("/home/justin/work/research/DPG/3-13-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/mc_weighted_small/Outfile.root");
	TFile* fdata = TFile::Open("/home/justin/work/research/DPG/3-13-21/PCAnalysis_RecoSIM/RUN_ON_DATA/pc_ExecutionDirectory/data_test/Outfile.root");
	TFile* fdata2 = TFile::Open("/home/justin/work/research/DPG/3-13-21/PCAnalysis_RecoSIM/RUN_ON_DATA/pc_ExecutionDirectory/data_test2/Outfile.root");

	TFile* fmc1 = TFile::Open("/home/justin/work/research/DPG/3-23-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/save/mc_weighted_L1neventweight.root");
	TFile* fmc2 = TFile::Open("/home/justin/work/research/DPG/3-23-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/save/mc_weighted_bias0Beventweight.root");

	TH1D* trueGeom = (TH1D*)fmc->Get("trueGeom");
	TH1D* trueConv = (TH1D*)fmc->Get("trueConv");
	TCanvas* c1 = new TCanvas("true1");
	TH1D* trueConv_ = (TH1D*)trueConv->Clone();
	trueConv_->Divide(trueGeom);
	trueConv_->Scale(9./7.);
	trueConv_->SetTitle(";R [cm]; x/x0");
	trueConv_->Draw("HIST E");
	
	TH1D* trueGeom_C = (TH1D*)fmc->Get("trueGeom_Coarse");
	TH1D* trueConv_C = (TH1D*)fmc->Get("trueConv_Coarse");
	TCanvas* c2 = new TCanvas("trueC");
	TH1D* trueConv_C_ = (TH1D*) trueConv_C->Clone();
	trueConv_C_->Divide(trueGeom_C);
	trueConv_C_->Scale(9./7.);
	trueConv_C_->Draw("HIST E");


	TH1D* nconvR_match = (TH1D*)fmc->Get("nconvR_match");
	TH1D* nconvR_match_eff = (TH1D*)nconvR_match->Clone();
	nconvR_match_eff->Divide(trueConv);
	double r,n;
	for(int i=1; i<=nconvR_match_eff->GetNbinsX(); i++){
		r = nconvR_match_eff->GetBinContent(i);
		n = trueConv->GetBinContent(i);
		nconvR_match_eff->SetBinError(i, sqrt( (r*(1-r))/n ));
	}
	nconvR_match_eff->Draw();

	TH1D* nconvR_match_C = (TH1D*)fmc->Get("nconvR_match_Coarse");
	nconvR_match_C->Divide(trueConv_C);
	nconvR_match_C->SetLineColor(kRed);
	nconvR_match_C->Draw("SAMES");

	TCanvas* c3 = new TCanvas("x0_","x0_");
	nconvR_match->Divide(trueGeom);
	nconvR_match->Divide(nconvR_match_eff);
	nconvR_match->Scale(9./7.);
	nconvR_match->Draw("HIST E");

/*
	TH1D* nconvR_data = (TH1D*) fdata->Get("nconvR");
	nconvR_data->Draw("HIST E");
	//nconvR_data->Divide(nconvR_match_eff);
	nconvR_data->Divide(trueGeom);
	nconvR_data->Scale(9./7.);
	//nconvR_data->Draw("HIST E");	


	//num hgn pc hist
	TH1D* nhgnpc_data = (TH1D*) fdata->Get("numHGNPCHist");
	TH1D* nhgnpc_mc = (TH1D*) fmc->Get("numHGNPCHist");
	TCanvas* c4 = new TCanvas("pccounts","pccounts" );
	nhgnpc_data->Draw("HIST E");
	nhgnpc_mc->Draw("HIST SAME E");


	TH1D* npc_data1 = (TH1D*) fdata->Get("numpcHist");
	TH1D* npc_mc = (TH1D*) fmc->Get("numpcHist");
	TH1D* npc_data2 = (TH1D*) fdata2->Get("numpcHist");

	npc_data1->Scale(1./npc_data1->Integral());
	npc_mc->Scale(1./npc_mc->Integral());
	npc_data2->Scale(1./npc_data2->Integral());

	npc_data1->SetLineColor(kBlue);
	npc_mc->SetLineColor(kRed);
	npc_data2->SetLineColor(kBlack);
	npc_data1->SetTitle("L1 Data;num PC");
	npc_data2->SetTitle("2018B Data");
	npc_mc->SetTitle("RecoSIM MC");

	TCanvas* c5 = new TCanvas("rawpc","rawpc");
	npc_data1->Draw("HIST E");
	npc_mc->Draw("HIST SAME E");
	npc_data2->Draw("HIST SAME E");
	c5->BuildLegend();
/*
	TCanvas* c6= new TCanvas("nconvr_comp","nconvr_comp");
	TH1D* nconvR_data_1 = (TH1D*) fdata->Get("nconvR");
	TH1D* nconvR_data_2 = (TH1D*) fdata2->Get("nconvR");
	TH1D* nconvR_mc_1 = (TH1D*) fmc1->Get("nconvR");
	TH1D* nconvR_mc_2 = (TH1D*) fmc2->Get("nconvR");
	nconvR_data_1->SetLineColor(kBlue);
	nconvR_mc_1->SetLineColor(kMagenta);
	nconvR_data_2->SetLineColor(kRed);
	nconvR_mc_2->SetLineColor(kOrange);
//	nconvR_data_1->Draw("HIST E");
//	nconvR_mc_1->Draw("HIST E SAME");
	nconvR_data_2->Draw("HIST E ");
	nconvR_mc_2->Draw("HIST E SAME");
	nconvR_data_1->Draw("HIST SAME E");
	 nconvR_mc_1->Draw("HIST E SAME");
*/
	

	TFile* fout = new TFile("hists.root","RECREATE");
//	fout->WriteTObject(trueConv);
//	fout->WriteTObject(trueConv_C);
/*
	TH1D* ngbp1 = (TH1D*)fmc->Get("Ng_BP1");
	TH1D* ngbp2 = (TH1D*)fmc->Get("Ng_BP2");
	TH1D* ncbp1 = (TH1D*)fmc->Get("Nc_BP1");
	TH1D* ncbp2 = (TH1D*)fmc->Get("Nc_BP2");

	ncbp1->Divide(ngbp1);
	ncbp1->Scale(9./7.);
	ncbp2->Divide(ngbp2);
	ncbp2->Scale(9./7.);
	TCanvas* c1 =new TCanvas("bp1","bp1");
	ncbp1->Draw();
	TCanvas* c2 =new TCanvas("bp2","bp2");
	ncbp2->Draw();
*/


/*
	//create fgeom hist
	Float_t Rbins[] = { 1,5,9,13,18,20 };
	Int_t  binnum = 5;
	double pi = 4*atan2(1,1);
	double z1 = -25.;
	double z2 = 25.;

	
	TH1D* fgeom = new TH1D("fgeom","fgeom;R cm",binnum,Rbins);
	TH1D* fgeomf = new TH1D("fgeomf","fgeom fine;R cm",40,0,20);
	for(int i=1; i<=binnum; i++){
		float r1 = Rbins[i-1]; 
		float r2 = Rbins[i];
		fgeom->SetBinContent(i, 2*pi*(z2-z1)*log(r2/r1)  );
	}

	for(int i=1; i<=fgeomf->GetNbinsX();i++){
		double r1 = fgeomf->GetBinLowEdge(i);
		double r2 = fgeomf->GetBinCenter(i)+ fgeomf->GetBinWidth(i)/2.;
		if(r1!=0)
		fgeomf->SetBinContent(i, 2*pi*(z2-z1)*log(r2/r1) );
		if(r1==0)
		fgeomf->SetBinContent(i,0.);
	}
/////////
	for(int i=1; i<=fgeomf->GetNbinsX(); i++){
		double bc = fgeomf->GetBinCenter(i);
		for(int j=0; j<binnum; j++){
			if( bc < Rbins[0] ){
				fgeomf->SetBinContent(i,fgeom->GetBinContent(1));
				break;
			}
			else if(bc >= Rbins[j] && bc<Rbins[j+1] ){
				fgeomf->SetBinContent(i, fgeom->GetBinContent(j+1));
				break;
			}
		}		
	}
*/
/*	TH1D* effRf = new TH1D("effRf","fine eff R;R cm",40,0,20);

	//Efficiency plots
	TH1D* effPtN =(TH1D*) fmc->Get("effPtN");
	TH1D* effPtD =(TH1D*) fmc->Get("effPtD");
	teff("effPt",effPtN,effPtD,";p_{T} [GeV];#varepsilon");

	TH1D* purityPtN = (TH1D*) fmc->Get("purityPtN");
	TH1D* purityPtD = (TH1D*) fmc->Get("purityPtD");
	teff("purPt",purityPtN,purityPtD,";p_{T} [GeV];Purity");
	

	//purityPtN->Divide(purityPtD);
	//purityPtN->Draw();
	//Purity Plots

	TH1D* effRN =(TH1D*) fmc->Get("effRN");
	TH1D* effRD =(TH1D*) fmc->Get("effRD");
	teff("effRRR", effRN, effRD, ";R [cm];#varepsilon");
	effRN->SetTitle(";R [cm];#varepsilon");
	effRN->Divide(effRD);
	double bc;
	double bcn;
	for(int i=1; i<=5; i++){
		bc = effRN->GetBinContent(i);
		bcn = effRD->GetBinContent(i);
		std::cout<<"bc bn"<<bc<<" "<<bcn<<std::endl;
		effRN->SetBinError(i, sqrt( bc*(1-bc)/bcn) );
	}
	effRN->Draw("E1");
//	TH1D* effR_ = (TH1D*) effR->CreateHistogram();

	TH1D* purityRN = (TH1D*) fmc->Get("purityRN");
	TH1D* purityRD = (TH1D*) fmc->Get("purityRD");
	teff("purR", purityRN, purityRD, ";R [cm];Purity");
//	effRN->Divide(effRD);
//	effRN->Draw("HIST E");

//	purityRN->Divide(purityRD);
//	purityRN->Draw();

	//make fine eff R
	for(int i=1; i<=effRf->GetNbinsX(); i++){
                double bc = effRf->GetBinCenter(i);
                for(int j=0; j<binnum; j++){
                        if( bc < Rbins[0] ){
                                effRf->SetBinContent(i,effRN->GetBinContent(1
));
				effRf->SetBinError(i,effRN->GetBinError(1));
                                break;
                        }
                        else if(bc >= Rbins[j] && bc<Rbins[j+1] ){
                                effRf->SetBinContent(i, effRN->GetBinContent(
j+1));
				effRf->SetBinError(i,effRN->GetBinError(j+1));
                                break;
                        }    
                }
        }

	


	TH1D* effXPN =(TH1D*) fmc->Get("effXPN");
	TH1D* effXPD =(TH1D*) fmc->Get("effXPD");
	teff("effxp",effXPN,effXPD);

	TH1D* purityXPN = (TH1D*) fmc->Get("purityXPN");
	TH1D* purityXPD = (TH1D*) fmc->Get("purityXPD");
	teff("purxp",purityXPN, purityXPD);

	//effXPN->Divide(effXPD);
	//effXPN->Draw();	

	//purityXPN->Divide(purityXPD);
	//purityXPN->Draw();


	TH1D* Npt_data = (TH1D*) fdata->Get("nconvPt");
	TH1D* Npt_mc = (TH1D*) fmc->Get("nconvPt");
	TH1D* Npt_mc_fake = (TH1D*) fmc->Get("nconvPt_fake");
	format3(Npt_data, Npt_mc, Npt_mc_fake);
	TCanvas* c1c = new TCanvas("c1c","c1c");
	Npt_data->Draw();
	Npt_mc->Draw("same");
	Npt_mc_fake->Draw("same");
	legend3(Npt_data,Npt_mc, Npt_mc_fake);

	TH1D* NR_data = (TH1D*) fdata->Get("nconvR");
        TH1D* NR_mc = (TH1D*) fmc->Get("nconvR");
	NR_mc->SetTitle(";R [cm];N_{conv}");
        TH1D* NR_mc_fake = (TH1D*) fmc->Get("nconvR_fake");
	format3(NR_data, NR_mc, NR_mc_fake);
	TCanvas* c2c = new TCanvas("c2c","c2c");
        //NR_data->Draw("HIST E");
        NR_mc->Draw("HIST E");
        NR_mc_fake->Draw("same");
	//legend3(NR_data, NR_mc, NR_mc_fake);
	legend2(NR_mc, NR_mc_fake);

	TH1D* NXP_data = (TH1D*) fdata->Get("nconvXP");
        TH1D* NXP_mc = (TH1D*) fmc->Get("nconvXP");
        TH1D* NXP_mc_fake = (TH1D*) fmc->Get("nconvXP_fake");
	format3(NXP_data, NXP_mc, NXP_mc_fake);
        TCanvas* c3c = new TCanvas("c3c","c3c");
        NXP_data->Draw();
        NXP_mc->Draw("same");
        NXP_mc_fake->Draw("same");
	legend3(NXP_data, NXP_mc, NXP_mc_fake);


	TH1D* x0_data = (TH1D*) NR_data->Clone();
	x0_data->SetTitle(";R [cm]; x/X_{0}");
	TH1D* x0_mc = (TH1D*) NR_mc->Clone();
	TH1D* x0_fake = (TH1D*)NR_mc_fake->Clone();
	format3( x0_data, x0_mc, x0_fake);

	fgeomf->Scale(1./fgeomf->Integral());
//	fgeomf->Scale(6159);
	TH1D* ngPrompt = (TH1D*) fmc->Get("ngPrompt");
	double Ng = ngPrompt->GetBinContent(1);
	fgeomf->Scale(Ng);

	TH1D* trueGeom = (TH1D*) fmc->Get("trueGeom");
	TH1D* trueConv = (TH1D*) fmc->Get("trueConv");
	trueConv->Divide(trueGeom);
	trueConv->Scale(9./7.);
	TCanvas* ctruth = new TCanvas("x0true","x0true");
	trueConv->Draw("HIST E");


	TH1D* trueGeomtest = (TH1D*) trueGeom->Clone();
	trueGeomtest->Scale(2);
//	x0_data->Divide(fgeomf);
//	x0_mc->Divide(fgeomf);
//	x0_fake->Divide(fgeomf);
	x0_data->Divide(trueGeom);
	x0_mc->Divide(trueGeomtest);
	x0_fake->Divide(trueGeomtest);
	x0_data->Divide(effRf);
	x0_mc->Divide(effRf);
	x0_fake->Divide(effRf);
	x0_data->Scale(9./7.);
	x0_mc->Scale(9./7.);
	x0_fake->Scale(9./7.);

	//x0_data->Scale(1./x0_data->Integral());
	// x0_mc->Scale(1./x0_mc->Integral());
	//  x0_fake->Scale(1./x0_mc->Integral());

	TCanvas* c4c = new TCanvas("c4c","c4c");
	//x0_data->Draw();
	x0_mc->Draw("HIST");
	x0_fake->Draw("SAME");
//	legend3(x0_data, x0_mc, x0_fake);
	legend2(x0_mc, x0_fake);
	TFile* fout = new TFile("hists.root","RECREATE");
	fout->WriteTObject(gCanv("c1",effPtN,"E1"));
	fout->WriteTObject(gCanv("c2",effRN,"HIST E"));
	fout->WriteTObject(gCanv("c3",effXPN,"E1"));
	fout->WriteTObject(gCanv("c4",fgeom,"HIST"));
	fout->WriteTObject(gCanv("c5",fgeomf,"HIST E"));
	fout->WriteTObject(gCanv("c6",effRf,"HIST E"));
	fout->WriteTObject(gCanv("c7",purityPtN,"E1"));
	fout->WriteTObject(gCanv("c8",purityRN,"E1"));
	fout->WriteTObject(gCanv("c9",purityXPN,"E1"));
	fout->WriteTObject(c1c);
	fout->WriteTObject(c2c);
	fout->WriteTObject(c3c);
	fout->WriteTObject(c4c);
	fout->WriteTObject(ctruth);
	fout->Close();
*/
}
