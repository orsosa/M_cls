#define M_survey_cls_cxx
#include "M_survey_cls.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void M_survey_cls::Loop()
{
//   In a ROOT session, you can do:
//      root> .L M_survey_cls_sim.C
//      root> M_survey_cls_sim t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  if (!QUIET)std::cout<<"processing...\n";
  if (!QUIET)std::cout.fill('.');
  if (!QUIET)std::cout<<"# trees to be processed: "<<fChain->GetNtrees()<<std::endl;

  Long64_t nentries = fChain->GetEntries();
  if (options.Contains("nmax:")){
    TString sub=options("nmax:[0-9]*");
    sub=sub(":[0-9]*");
    sub.ReplaceAll(":","0");
    if (sub.IsDigit())
      nentries = sub.Atoll();
  }
  Long64_t nbytes = 0, nb = 0;
  TDirectory *d;
  for (Long64_t jentry = 0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    TH1D *h;
    TProfile *pf;
    DIR["kin"]->GetObject("hNpair_all",h);
    h->Fill(npart);
    //    if (!(npart == 1)) continue;
    if ( !(DIS() && eFID_ec() && eFID_dc() && ePID()) ) continue;
    fillEvHistos();
    
    //for (int k = 0; k<npart; k++){
    for (int k = 0; k<1; k++){ // only first pair
      if ( !(FWD(k) && MMCut(k) && rhoCut(k) && CF(k) && minEnergy(k) && piFID_ec(k) && pipFID_dc(k) && pimFID_dc(k) && pi0PID(k) && DZCut(k)) ) continue; // kin limits.
      //if ( !(FWD(k) && CF(k) && minEnergy(k) && piFID_ec(k) && pipFID_dc(k) && pimFID_dc(k) && pi0PID(k)) ) continue; // kin limits.
      //     
      
      fillPartHistos(k);
      
      TString hnameM = "hM";
      TString hnameMx = "hMx";
      TString hnamedz = "hdz";
      TString hmpi0name = "hmpi0";
      TString infix = "";
      TString suff = "";
      TString ttlsuf = "";
      Bool_t binFound;
      for (auto& x: bedg){
	binFound = kFALSE;
	TString bn = x.first;
	brv[bn]->GetNdata();
	infix += bn.Data();
	for (int n=0;n<x.second.size()-1;n++){
	  if (x.second[n] < brv[bn]->EvalInstance(k) && brv[bn]->EvalInstance(k)< x.second[n+1]){
	    suff += Form("%d",n);
	    ttlsuf +=  Form("(%.2f<%s<%.2f)",x.second[n], bn.Data(), x.second[n+1]);
	    binFound = kTRUE;
	    break;
	  }
	}
	if (!binFound) break;
      }
      if (!binFound) break;
      
      hnameM = TString(DIR[hnameM]->GetName()) + "/" + hnameM + "_" + infix +"_"+suff;
      fillHist(hnameM,M[k]);
      
      if (options.Contains("pi0")){
	hmpi0name = TString(DIR[hmpi0name]->GetName()) + "/" + hmpi0name + "_" + infix +"_"+suff;
	fillMpi0Hist(hmpi0name,k);
      }
      hnameMx = TString(DIR[hnameMx]->GetName()) + "/" + hnameMx + "_" + infix +"_"+suff;
      fillHist(hnameMx,sqrt(Mx2[k]));
      hnamedz = TString(DIR[hnamedz]->GetName()) + "/" + hnamedz + "_" + infix +"_"+suff;
      fillHist(hnamedz,GetDZ(k));
    }

    if (!QUIET) std::cout<<std::setw(15)<<float(jentry+1)/nentries*100<<" %"<<"\r";
    if (!QUIET) std::cout.flush();

  }
  ofile->Write("",TObject::kOverwrite);
  ofile->Close();
  bm->Show("main");

}

Bool_t M_survey_cls::ePID()
{
  return  (vze<20) &&(-2.5<e_chi2pid&&e_chi2pid<2.5);
  //return  (vze<15);
}

Bool_t M_survey_cls::DIS()
{
  return (Q2>1)&&(W>2)&&(y<0.8);
}

Bool_t M_survey_cls::eFID_ec()
{
  return (20<e_pcal_lv)&&(20<e_pcal_lw);
}

Bool_t M_survey_cls::piFID_ec(int k)
{
  Bool_t fid1 = true;
  if (!options.Contains("piFID_ec")) return true;
  if (options.Contains("pi0")){
    fid1 = (20<det_pcal_lv[k][1])&&(20<det_pcal_lw[k][1])&&(20<det_pcal_lv[k][2])&&(20<det_pcal_lw[k][2]);;
  }
  else{
    fid1 = (20<det_pcal_lv[k][1])&&(20<det_pcal_lw[k][1]);
  }
  
  return (20<det_pcal_lv[k][0])&&(20<det_pcal_lw[k][0])&&fid1;
}

Bool_t M_survey_cls::eFID_dc()
{
  return true;
}

Bool_t M_survey_cls::pipFID_dc(int k)
{
  return true;
}

Bool_t M_survey_cls::pimFID_dc(int k)
{
  return true;
}

Bool_t M_survey_cls::MMCut(int k)
{
  Bool_t ret;
  if (options.Contains("noMxCut"))
    ret = true;
  else if (options.Contains("rho-"))
    if (options.Contains("rgb"))
      ret = sqrt(Mx2[k])>1.33; 
    else 
      //ret = sqrt(Mx2[k])>1.68; // rga and rgaob, Delta--
      ret = sqrt(Mx2[k])>1.33; // rga and rgaob, for consistency.

  else if (options.Contains("rgaob"))
    ret = sqrt(Mx2[k])>1.16;
  else if (options.Contains("rgb"))
    ret = sqrt(Mx2[k])>1.13;

  else
    ret = sqrt(Mx2[k])>1.10;//rga and sim
  if (options.Contains("excl")) ret=!ret;
  return ret;

}

Bool_t M_survey_cls::rhoCut(int k)
{
  Bool_t ret = false; 
  if (!options.Contains("rhoCut"))
    ret = true;
  
  else if (options.Contains("rgb"))
    ret = 0.624<M[k]&&M[k]<0.902;

  else if (options.Contains("rgaob"))
    ret = 0.634<M[k]&&M[k]<0.884;
  else
    ret = 0.628<M[k]&&M[k]<0.897; //rga

  if (options.Contains("rhoCutComp"))
    ret = !ret; 

  return ret;
}


Bool_t M_survey_cls::FWD(int k)
{
  bool fwd1 = true;
  if (options.Contains("pi0")){
    fwd1 = ( ((int)det_statPart[k][1]%4000)/2000 >= 1 )&&( ((int)det_statPart[k][2]%4000)/2000 >= 1 );
  }
  else{
    fwd1 = ( ((int)det_statPart[k][1]%4000)/2000 >= 1 );
  }
  
  return ( ((int)det_statPart[k][0]%4000)/2000 >= 1 ) && fwd1;
}

Bool_t M_survey_cls::CF(int k)
{
  Float_t xF_1 = -1111;
  Float_t E_1 = -1111;
  if (options.Contains("pi0")){
    xF_1 = getxF(pdata_e[k][1] + pdata_e[k][2], pdata_px[k][1] + pdata_px[k][2], pdata_py[k][1] + pdata_py[k][2], pdata_pz[k][1] + pdata_pz[k][2]);
    E_1 = (pdata_e[k][1] + pdata_e[k][2]);
  }
  else{
    xF_1 = xFm1[k];
    E_1 = pdata_e[k][1];
  }
  
  return (xFm0[k]>0) && (xF_1>0) && (pdata_e[k][0]/Nu>0.1) && (E_1/Nu>0.1) && ((pdata_e[k][0] + E_1)/Nu<0.95);
}

Bool_t M_survey_cls::minEnergy(int k)
{
  Float_t xF_1 = -1111;
  Float_t E_1 = -1111;
  if (options.Contains("pi0")){
    E_1 = (pdata_e[k][1] + pdata_e[k][2]);
  }
  else{
    E_1 = pdata_e[k][1];
  }

  return (pdata_e[k][0]>1.0) && (E_1>1.0);
}


Float_t M_survey_cls::getALU(TString hpname, TString hnname, TString pv, TString tv, Float_t &val, Float_t &err){
  TH1D *hp,*hn, *hs, *hm, *hALU;
  TFitResultPtr res;
  TF1 *ff = new TF1("ff"," [A]*sin(x*TMath::DegToRad())",0,360);

  gDirectory->GetObject(hpname,hp);
  gDirectory->GetObject(hnname,hn);
  hp->Sumw2();
  hn->Sumw2();
  hs = (TH1D *)hp->Clone("hs_" + pv);
  hm = (TH1D *)hp->Clone("hm_" + pv);
  hALU = (TH1D * )hp->Clone("hALU_" + pv);
  hALU -> SetTitle("A_{LU}^{sin(" + tv + ")} = (N^{+} - N^{-}) / (N^{+} + N^{-})");
  hs->Add(hp,hn,1,1);
  hm->Add(hp,hn,1,-1);
  hALU->Divide(hm,hs);
  ff->SetParameter(0,0.01);
  if (!QUIET)std::cout<<tv<<std::endl;
  res = hALU->Fit(ff,"Rs+0");
  if ( (0.4<(ff->GetChisquare()/ff->GetNDF())&&(ff->GetChisquare()/ff->GetNDF()<2.5)) ){
    val = ff->GetParameter(0);
    err = ff->GetParError(0);
  }
  else{
    val = -111;
    err = 0;
  }

  return ff->GetChisquare()/ff->GetNDF();
}

Float_t M_survey_cls::getALU2D(TString bn){
  TString psname, nsname;
  TH2D *hp,*hn;
  TH1D *halup,*halun,*halus;

  ////  pltv ////
  psname = "hsin"+pltv + "_" + bn + "_p";
  nsname = "hsin"+pltv + "_" + bn + "_n";
  gDirectory->GetObject(psname,hp);
  gDirectory->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = hp->ProfileX()->ProjectionX("hpALU_"+ pltv + "_" + bn);
  halun = hn->ProfileX()->ProjectionX("hnALU_"+ pltv + "_" + bn);
  halus = (TH1D *)halup->Clone("hsALU_"+ pltv + "_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(" + ttlv + ")>^{+}");
  halun->SetTitle("-2<sin(" + ttlv + ")>^{-}");
  halus->SetTitle("(<sin(" + ttlv + ")>^{+} - <sin(" + ttlv+ ")>^{-})");
  configHisto(halup,bn,"A_{LU}(<sin(" + ttlv + ")>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(" + ttlv + ")>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(" + ttlv + ")>^{+/-})",kMagenta+1,kFullStar);
    
  ///// end pltv /////
  //// phiH ////
  psname = "hsinphiH_" + bn + "_p";
  nsname = "hsinphiH_" + bn + "_n";
  gDirectory->GetObject(psname,hp);
  gDirectory->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = (TH1D*)hp->ProfileX("hpALU_phiH_" + bn);
  halun = (TH1D*)hn->ProfileX("hnALU_phiH_" + bn);
  halus = (TH1D*)halup->Clone("hsALU_phiH_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(#phi_{H})>^{+}");
  halun->SetTitle("-2<sin(#phi_{H})>^{-}");
  halus->SetTitle("(<sin(#phi_{H})>^{+} - <sin(#phi_{H})>^{-})");
  configHisto(halup,bn,"A_{LU}(<sin(#phi_{H})>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(#phi_{H})>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(#phi_{H})>^{+/-})",kMagenta+1,kFullStar);

  //// end phiH ///
  ////  phiH 0 ////
  psname = "hsinphiH_0_" + bn + "_p";
  nsname = "hsinphiH_0_" + bn + "_n";
  gDirectory->GetObject(psname,hp);
  gDirectory->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = hp->ProfileX()->ProjectionX("hpALU_phiH_0_" + bn);
  halun = hn->ProfileX()->ProjectionX("hnALU_phiH_0_" + bn);
  halus = (TH1D*)halup->Clone("hsALU_phiH_0_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(#phi_{H})>^{+}");
  halun->SetTitle("-2<sin(#phi_{H})>^{-}");
  halus->SetTitle("(<sin(#phi_{H})>^{+} - <sin(#phi_{H})>^{-})");

  configHisto(halup,bn,"A_{LU}(<sin(#phi_{H})>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(#phi_{H})>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(#phi_{H})>^{+/-})",kMagenta+1,kFullStar);
  ///// end phiH 0 /////

  ////  phiH 1 ////
  psname = "hsinphiH_1_" + bn + "_p";
  nsname = "hsinphiH_1_" + bn + "_n";
  gDirectory->GetObject(psname,hp);
  gDirectory->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = hp->ProfileX()->ProjectionX("hpALU_phiH_1_" + bn);
  halun = hn->ProfileX()->ProjectionX("hnALU_phiH_1_" + bn);
  halus = (TH1D*)halup->Clone("hsALU_phiH_1_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(#phi_{H})>^{+}");
  halun->SetTitle("-2<sin(#phi_{H})>^{-}");
  halus->SetTitle("(<sin(#phi_{H})>^{+} - <sin(#phi_{H})>^{-})");

  configHisto(halup,bn,"A_{LU}(<sin(#phi_{H})>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(#phi_{H})>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(#phi_{H})>^{+/-})",kMagenta+1,kFullStar);
  
  ///// end phiH 1 /////

  return 0;
}

Int_t M_survey_cls::fitMpi0Hist(TString hname){
  TH1D *h;
  gDirectory->GetObject(hname,h);
  if (h->GetEntries()==0) return 1;
  TF1 *ff = new TF1 (hname +"_ff","gaus+pol2(3)",MPI0_LL,MPI0_HL);
  ff->SetParName(0,"C"); ff->SetParName(1,"m_{#gamma#gamma}"); ff->SetParName(2,"#sigma_{#gamma#gamma}"); ff->SetParName(3,"a_{0}"); ff->SetParName(4,"a_{1}"); ff->SetParName(5,"a_{2}");
  ff->SetParLimits(0,0.,1e6);
  ff->SetParLimits(1,0.12,0.14);
  ff->SetParLimits(2,0.007,0.03);
  ff->SetParLimits(3,0.,1e4);
  ff->SetParLimits(5,-1e5,0);
  ff->SetParameters(800,0.135,0.013,392,409,-5447);
  ff->SetLineColor(kRed);
  h->Fit(ff,"RNB");
  h->GetListOfFunctions()->Add(ff);
  
  return 0;
}

Int_t M_survey_cls::fillMpi0Hist(TString hname, Int_t k){
  TH1D *h;
  ofile->GetObject(hname,h);
  Float_t mpi0 = sqrt(2*(pdata_e[k][1]*pdata_e[k][2] - pdata_px[k][1]*pdata_px[k][2] - pdata_py[k][1]*pdata_py[k][2] - pdata_pz[k][1]*pdata_pz[k][2]));

  h->Fill(mpi0);
  
  return 0;
}

Int_t M_survey_cls::fillHist(TString hname, Float_t value){
  TH1D *h;
  ofile->GetObject(hname,h);
  h->Fill(value);
  return 0;
}

Int_t M_survey_cls::fillHist2D(TString hname, Float_t x, Float_t y){
  TH2D *h;
  if (*helicity == -1) { // positive helicity, it is flipped!
    ofile->GetObject(hname + "_p",h);
    h->Fill(x,y);
  }
  else if(*helicity == 1){ // negative helicity, it is flipped!
    ofile->GetObject(hname + "_n",h);
    h->Fill(x,-y);
  }
  else
    return 1;
  
  return 0;
}

Int_t M_survey_cls::configHisto(TH1D *h, TString xtitle, TString ytitle, Color_t c, EMarkerStyle ms){

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetMarkerStyle(ms);
  h->SetMarkerSize(1.3);
  
  return 0;
  
}

Int_t M_survey_cls::LoadElecFIDPar(){
  std::ifstream fpar("/home/orsosa/Dropbox/INFN_work/Phys_ana/PID/DC_elec_par.txt");

  char junk[100];
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);

  for (int k=0;k<NSECTORS;k++)
  {
    fpar>>pl0_e[k]>>pl1_e[k]>>pl2_e[k]>>pl3_e[k]>>pr0_e[k]>>pr1_e[k]>>pr2_e[k]>>pr3_e[k];
  }
  fpar.close();
  return 0;
}


Int_t M_survey_cls::setStyle(){

  myStyle  = new TStyle("orsosaStyle","My Root Styles");
  myStyle->SetHistMinimumZero(0);
  myStyle->SetPalette(1,0);
  //myStyle->SetPalette(kBlueYellow);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(0);

  myStyle->SetStatColor(0);
  myStyle->SetOptStat("e");

  myStyle->SetLabelSize(0.05,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.05,"xyz"); // size of axis title font
  myStyle->SetTitleOffset(0.75,"xyz"); // axis title offset 
  myStyle->SetTitleFont(22,"xyz"); // font option
  myStyle->SetTitleFont(22,"a"); // pad font option
  myStyle->SetLabelFont(22,"xyz");


  myStyle->SetLabelSize(0.02,"z"); // size of axis value font
  myStyle->SetLabelOffset(-0.03,"z"); // size of axis value font
  myStyle->SetTickLength(0.002,"z"); // size of axis value font

  
  // Stat and legend fonts
  myStyle->SetStatFont(22); 
  myStyle->SetLegendFont(22); 
  // hiostogram style
  myStyle->SetHistLineWidth(2);
  myStyle->SetCanvasDefH(768);
  myStyle->SetCanvasDefW(1024);

  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadLeftMargin(0.1);
  myStyle->SetPadRightMargin(0.075);

  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);

  myStyle->SetFrameBorderMode(0);
  
  myStyle->SetGridStyle(3);
  myStyle->SetGridWidth(2);
  myStyle->SetPadGridX(kTRUE);
  myStyle->SetPadGridY(kTRUE);

  
  gROOT->SetStyle("orsosaStyle"); //uncomment to set this style
  
  return 0;
}

Bool_t M_survey_cls::pi0PID(Int_t k){
  if (!options.Contains("pi0")) return kTRUE;
  Float_t minth_brem = 10;
  Float_t minE = 0.5;
  Bool_t ret = true;
  ret = ret && acos((pdata_px[k][1]*Pex + pdata_py[k][1]*Pey + pdata_pz[k][1]*Pez)/pdata_e[k][1]/(Nu/y-Nu))*TMath::RadToDeg()>minth_brem;
  ret = ret && acos((pdata_px[k][2]*Pex + pdata_py[k][2]*Pey + pdata_pz[k][2]*Pez)/pdata_e[k][2]/(Nu/y-Nu))*TMath::RadToDeg()>minth_brem; 
  ret = ret && (pdata_e[k][1]>minE)&&(pdata_e[k][2]>minE);
  Float_t mpi0 = sqrt(2*(pdata_e[k][1]*pdata_e[k][2] - pdata_px[k][1]*pdata_px[k][2] - pdata_py[k][1]*pdata_py[k][2] - pdata_pz[k][1]*pdata_pz[k][2]));
  ret = ret && det_pcal_lv[k][1]>0&&det_pcal_lv[k][1]>0;
  ret = ret && det_pcal_lw[k][2]>0&&det_pcal_lw[k][2]>0;
  ret = ret && 10<pdata_phiHs[k][1]&&pdata_phiHs[k][1]<350;
  ret = ret && 10<pdata_phiHs[k][2]&&pdata_phiHs[k][2]<350;
  TH1D *h;
  gDirectory->GetObject("hmpi0",h);
  if (ret) h->Fill(mpi0);
  ret = ret && MPI0_LL<mpi0&&mpi0<MPI0_HL;
  gDirectory->GetObject("hmpi0_cut",h);
  if (ret) h->Fill(mpi0);
  if (options.Contains("noPi0Cut")) ret = kTRUE;
  return ret;
}


Float_t M_survey_cls::getxF(Float_t E, Float_t Px, Float_t Py, Float_t Pz){
  Float_t kMprt = 0.93827;
  Float_t kEbeam = Nu/y;
  Float_t P2 = Px*Px + Py*Py + Pz*Pz;
  Float_t cospq = ((kEbeam-Pez)*Pz - Pex*Px - Pey*Py)/( sqrt((Q2 + Nu*Nu)*P2) );
  Float_t Pt2 = P2*(1-cospq*cospq);
  Float_t Pl2 = P2*cospq*cospq;
  Float_t Pl = sqrt(P2)*cospq;
  ////// LORENTZ BOOST //////////
  Float_t b=TMath::Sqrt(Q2 + Nu*Nu)/(Nu + kMprt);
  Float_t g=(Nu + kMprt)/W;

  Float_t PlCM = g*(Pl - b*E);
      
  Float_t xFm = 2*PlCM/W;
  return xFm;
}

Double_t M_survey_cls::get_phiR(int k){
  Double_t phi = phiR[k];
  Float_t kEbeam = Nu/y;
  if (options.Contains("pi0")){
    TVector3 Ph(Phx[k],Phy[k],Phz[k]);
    TVector3 P0(pdata_px[k][0],pdata_py[k][0],pdata_pz[k][0]);
    TVector3 P1(pdata_px[k][1]+pdata_px[k][2],pdata_py[k][1]+pdata_py[k][2],pdata_pz[k][1]+pdata_pz[k][2]);
    TVector3 q(-Pex,-Pey,kEbeam-Pez);
    TVector3 k_in(0,0,kEbeam);

    TVector3 Ph_u = Ph.Unit();
    TVector3 R = P0 - P1;
    R=R*0.5;
    TVector3 RT = R-(R*Ph_u)*Ph_u;

    Float_t qxkRT_sign = q.Cross(k_in)*RT;
    qxkRT_sign /= TMath::Abs(qxkRT_sign);

    // Float_t qxkST_sign = q.Cross(k_in)*ST;
    // qxkST_sign /= TMath::Abs(qxkST_sign);

    Float_t qxkqxRT = (q.Cross(k_in))*(q.Cross(RT));
    Float_t qxkqxRT_max = (q.Cross(k_in)).Mag()*(q.Cross(RT)).Mag();
    
    phi =  qxkRT_sign*TMath::ACos(qxkqxRT/qxkqxRT_max)*TMath::RadToDeg();
    phi=phi<0?phi+360:phi;
  }
    
  return phi;
}

Float_t M_survey_cls::GetDZ(int k){
  Double_t DZ = (pdata_e[k][0] - pdata_e[k][1])/Nu;
  if (options.Contains("pi0"))
    DZ = (pdata_e[k][0]-(pdata_e[k][1]+pdata_e[k][2]))/Nu;

  return DZ;
}

Bool_t M_survey_cls::DZCut(int k)
{
  Double_t DZ = (pdata_e[k][0] - pdata_e[k][1])/Nu;
  Bool_t ret = false;
  if (!options.Contains("DZCut"))
    ret = true;
  else 
  {
    if (options.Contains("pi0"))
      DZ = (pdata_e[k][0]-(pdata_e[k][1]+pdata_e[k][2]))/Nu;
    
    ret =  -0.2<DZ&&DZ<0.2; 
  }
  
  return ret;
}


Double_t M_survey_cls::get_phiH1(int k){
  Double_t phi = pdata_phiHs[k][1];
  Double_t kEbeam = Nu/y;
  if (options.Contains("pi0")){
    TVector3 Phv(pdata_px[k][1] + pdata_px[k][2], pdata_py[k][1] + pdata_py[k][2], pdata_pz[k][1] + pdata_pz[k][2]);
    TVector3 q(-Pex,-Pey,kEbeam-Pez);
    TVector3 k_in(0,0,kEbeam);
    Float_t qxkPhv_sign = q.Cross(k_in)*Phv;
    qxkPhv_sign /= TMath::Abs(qxkPhv_sign);
    
    Float_t qxkqxPhv = (q.Cross(k_in))*(q.Cross(Phv));
    Float_t qxkqxPhv_max = (q.Cross(k_in)).Mag()*(q.Cross(Phv)).Mag();
    
    phi =  qxkPhv_sign*TMath::ACos(qxkqxPhv/qxkqxPhv_max)*TMath::RadToDeg();
    phi=phi<0?phi+360:phi;
  }

  return phi;
}
