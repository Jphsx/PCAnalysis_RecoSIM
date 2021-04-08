

void mcdataplot2(){

	Float_t Rbins[] = { 1,5,9,13,18,20 };
        Int_t  binnum = 5;
        double pi = 4*atan2(1,1);


	TFile* f = TFile::Open("/home/justin/work/research/DPG/4-07-21/PCAnalysis_RecoSIM/pc_ExecutionDirectory/local_test/Outfile.root");

	TH1D* effRN = (TH1D*) f->Get("effRN");
	TH1D* effRD = (TH1D*) f->Get("effRD");

	TH1D* effR = (TH1D*) effRN->Clone();
	effR->Divide(effRD);

	TH1D* effRf = (TH1D*) f->Get("nconvR_all");
	effRf = (TH1D*) effRf->Clone();
	for(int i=1; i<=effRf->GetNbinsX(); i++){
                double bc = effRf->GetBinCenter(i);
                for(int j=0; j<binnum; j++){
                        if( bc < Rbins[0] ){
                                effRf->SetBinContent(i,effR->GetBinContent(1
));
                                effRf->SetBinError(i,effR->GetBinError(1));
                                break;
                        }
                        else if(bc >= Rbins[j] && bc<Rbins[j+1] ){
                                effRf->SetBinContent(i, effR->GetBinContent(
j+1));
                                effRf->SetBinError(i,effR->GetBinError(j+1));
                                break;
                        }    
                }
        }

       
	TH1D* nconvR_match = (TH1D*) f->Get("nconvR_all");
	nconvR_match = (TH1D*) nconvR_match->Clone();
	
	TH1D* trueGeom = (TH1D*) f->Get("trueGeom");	

	nconvR_match->Divide(trueGeom);
	nconvR_match->Divide(effRf);
	nconvR_match->Scale(9./7.);
	nconvR_match->Draw("hist e");
}
