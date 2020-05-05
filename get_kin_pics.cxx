int get_kin_pics(TString fname = "", TString outdir = ""){
  TCanvas *c = new TCanvas();
  TLine *line;
  const int Nq2 = 4, Nx = 7, Nz = 9, NpT = 6;
  Double_t Q2[Nq2] = {1.0, 2.2, 3.6, 9.0};
  Double_t x[Nx]  = {0.05, 0.16, 0.22, 0.28, 0.36, 0.46, 0.70};
  Double_t z[Nz]  = {0.25, 0.33, 0.41, 0.49, 0.58, 0.66, 0.74, 0.82, 0.90};
  Double_t pT[NpT] = {0.00, 0.15, 0.30, 0.45, 0.80, 1.20};

  TH2D *h2;
  TH1D *h;
  //Q2 - x lines
  Float_t min_x;
  Float_t max_x;
  Float_t min_y;
  Float_t max_y;
  Int_t Nmax;
  Double_t *arr;
  TFile *infile;
  TTree *t;
  // Getting zx plot
  infile = new TFile(fname);
  /*  infile->GetObject("outdata",t);
  if (!t)
  {
    std::cout<<"file doesn't cointain outdata tree!\n";
    return 1;
  }
  h2 = new TH2D("hxz","x_{B} vs z",150,0.23,0.92,150,0,0.75);
  t->Draw("Xb:Z>>hxz","Q2>1&&W>2&&y<0.8&&(statPart[][0]%4000)/2000>=1&&(statPart[][1]%4000)/2000>=1&&e[][0]>1&&e[][1]>1&&e[][0]/Nu>0.1&&e[][1]/Nu>0.1&&(e[][0]+e[][1])/Nu<0.95");
  */
  ////
  infile->GetObject("kinematics/hxz",h2);
  h2->GetXaxis()->SetRangeUser(0.23,0.92);
  h2->GetYaxis()->SetRangeUser(0.04,0.73);
  h2->Draw("colz");

  min_y = x[0]; max_y = x[Nx-1];
  min_x = z[0]; max_x = z[Nz-1];

  Nmax = Nx;
  arr = x;
  for (int k=0; k<Nmax;k++){
    line = new TLine(min_x,arr[k],max_x,arr[k]);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    line->Draw();
  }
  Nmax = Nz;
  arr = z;
  for (int k=0; k<Nmax;k++){
    line = new TLine(arr[k],min_y,arr[k],max_y);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    line->Draw();
  }

  c->SaveAs(outdir + "/xz_bins.gif");
  c->SaveAs(outdir + "/xz_bins.C");

  std::cout<<"end"<<std::endl;
  infile->Close();
  return 0;
}
