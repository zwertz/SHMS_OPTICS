#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void make_hist_shms_optics(Int_t nrun=1813,Bool_t CutYtarFlag=kTRUE,Bool_t CutYpFpYFpFlag=kTRUE,Bool_t CutXpFpXFpFlag=kTRUE,Int_t FileID=-2){
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  //  Get info for that optics run
  TString OpticsFile = "DATfiles/list_of_optics_run.dat";
  ifstream file_optics(OpticsFile.Data());
  TString opticsline;
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Int_t NumFoil=0;
  TString temp;
  //
  vector <Double_t> ztar_foil;
  Int_t ndelcut=-1;
  vector<Double_t > delcut;
  vector<Double_t > delwidth;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ndelcut = temp.Atoi();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut-1;nd++) {
        temp.ReadToDelim(file_optics,',');
      	delcut.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delcut.push_back(temp.Atof());
      for (Int_t nw=0;nw<ndelcut-1;nw++) {
	temp.ReadToDelim(file_optics,',');
	delwidth.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delwidth.push_back(temp.Atof());
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
  //
  TString inputroot;
  TString outputhist;
  inputroot = Form("ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_%d.root",OpticsID.Data(),FileID);

  outputhist=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
  cout << " input root = " << inputroot << endl;
  TObjArray HList(0);
  //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << "Ytar Cut file = " << YtarDeltaCutFile << endl;
    for (Int_t nc=0;nc<NumFoil;nc++) {
      fYtarDeltaCut->cd();
      TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
      if (tempcut) {
	Int_t npt = tempcut->GetN();
	cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
	ytar_delta_cut.push_back(tempcut);
      } else {
	cout << " No hYtarDelta_cut = " << nc << endl;
      }
    }
  }
  //
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_yfp_cut;
  vector<vector<vector<Int_t> > > ypfp_yfp_cut_flag;
  ypfp_yfp_cut.resize(NumFoil);
  ypfp_yfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    ypfp_yfp_cut[nf].resize(ndelcut);
    ypfp_yfp_cut_flag[nf].resize(ndelcut);
  }
  if (CutYpFpYFpFlag) {
    outCutFile=Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<11;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
//
  TString xpfp_xfp_outCutFile;
  TFile *xpfp_xfp_fcut;
  vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
  vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
  xpfp_xfp_cut.resize(NumFoil);
  xpfp_xfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    xpfp_xfp_cut[nf].resize(ndelcut);
    xpfp_xfp_cut_flag[nf].resize(ndelcut);
  }
  if (CutXpFpXFpFlag) {
    xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
    cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
    xpfp_xfp_fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<11;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
  TFile *fsimc = new TFile(inputroot); 
  TTree *tsimc = (TTree*) fsimc->Get("T");
  // Define branches
  Double_t  sumnpe;
  tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
  //Double_t  sumhgnpe;
  //tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
  Double_t  etracknorm;
  tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
  Double_t  ytar;
  tsimc->SetBranchAddress("P.gtr.y",&ytar);
  Double_t  xtar;
  tsimc->SetBranchAddress("P.gtr.x",&xtar);
  Double_t  reactx;
  tsimc->SetBranchAddress("P.react.x",&reactx);
  Double_t  reacty;
  tsimc->SetBranchAddress("P.react.y",&reacty);
  Double_t  reactz;
  tsimc->SetBranchAddress("P.react.z",&reactz);
  Double_t  delta;
  tsimc->SetBranchAddress("P.gtr.dp",&delta);
  Double_t  yptar;
  tsimc->SetBranchAddress("P.gtr.ph",&yptar);
  Double_t  xptar;
  tsimc->SetBranchAddress("P.gtr.th",&xptar);
  Double_t  yfp;
  tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
  Double_t  ypfp;
  tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
  Double_t  xfp;
  tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
  Double_t  xpfp;
  tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
  Double_t  ysieve;
  tsimc->SetBranchAddress("P.extcor.ysieve",&ysieve);
  Double_t  xsieve;
  tsimc->SetBranchAddress("P.extcor.xsieve",&xsieve);
  Double_t  xbpm_tar;
   tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpm_tar);
   Double_t  ybpm_tar;
   tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpm_tar);
   Double_t fr_xa;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&fr_xa);
    Double_t fr_ya;
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&fr_ya);
   
  // Define histograms
  TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etotnorm ; Counts",nrun),100,0.,2.);
  HList.Add(hetot);
  TH1F *hngsum = new TH1F("hngsum",Form("Run %d ; NG Npe SUM ; Counts",nrun),100,0.,40.);
  HList.Add(hngsum);
  TH1F *hytar = new TH1F("hytar",Form("Run %d ; Ytar; Counts",nrun),500,-10.,20.);
  HList.Add(hytar);
  TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar; Counts",nrun),500,-35.,25.);
  HList.Add(hztar);
  TH1F *hztarg = new TH1F("hztarg",Form("Run %d ; Ztarg; Counts",nrun),500,-35.,25.);
  HList.Add(hztarg);
  TH1F *hztarCalc = new TH1F("hztarCalc",Form("Run %d ; ZtarCalc; Counts",nrun),500,-35.,25.);
  HList.Add(hztarCalc);
  TH2F *hXptarDelta = new TH2F("hXptarDelta",Form("Run %d ; Xptar ; Delta",nrun),120,-.06,.06,100,-15.,25.);
  HList.Add(hXptarDelta);
  TH2F *hYptarDelta = new TH2F("hYptarDelta",Form("Run %d ; Yptar ; Delta",nrun),120,-.04,.04,100,-15.,25.);
  HList.Add(hYptarDelta);
  TH2F *hYtarDelta = new TH2F("hYtarDelta",Form("Run %d ; Ytar ; Delta",nrun),100,-15.,15.,100,-15.,25.);
  HList.Add(hYtarDelta);
  //
  TH2F *hYpFpYFp_all = new TH2F("hYpFpXFp_all",Form("Run %d ; Ypfp ; Yfp",nrun),100,-.04,.04,100,-40.,40.);
  HList.Add(hYpFpYFp_all);
  TH2F *hYFpXFp_all = new TH2F("hYFpXFp_all",Form("Run %d ; Yfp ; Xfp",nrun),100,-40,40,100,-40.,40.);
  HList.Add(hYFpXFp_all);
  TH2F *hXpFpXFp_all = new TH2F("hXpFpXFp_all",Form("Run %d ; Xpfp ; Xfp",nrun),100,-.1,.1,100,-40.,40.);
  HList.Add(hXpFpXFp_all);
  TH2F *hYtarYptar = new TH2F("hYtarYptar",Form("Run %d ; Yptar ; Ytar",nrun),100,-.05,.05,100,-10.,20.);
  HList.Add(hYtarYptar);
  TH2F *hZtarDelta = new TH2F("hZtarDelta",Form("Run %d ; Ztar ; Delta",nrun),100,-35.,25.,100,-15.,25.);
  HList.Add(hZtarDelta);
  TH1F *hxbpm_tar = new TH1F("hxbpm_tar",Form("Run %d ; Xbpm_tar ; Counts",nrun),100,-2.,2.);
  HList.Add(hxbpm_tar);
  TH1F *hybpm_tar = new TH1F("hybpm_tar",Form("Run %d ; Ybpm_tar ; Counts",nrun),100,-2.,2.);
  HList.Add(hybpm_tar);
  TH2F *hZtarFrXa = new TH2F("hZtarFrXa",Form("Run %d ; fr_xa ; Ztar",nrun),100,-0.35,0.25,100,-11.,11.);
  HList.Add(hZtarFrXa);
  TH2F *hZtarFrYa = new TH2F("hZtarFrYa",Form("Run %d ; fr_ya ; Ztar",nrun),100,-0.35,0.25,100,-11.,11.);
  HList.Add(hZtarFrYa);
  //
  vector <TH2F*> hYsDelta;
  hYsDelta.resize(NumFoil);
  vector <TH2F*> hXsDelta;
  hXsDelta.resize(NumFoil);
  vector <TH2F*> hYpFpYFp;
  hYpFpYFp.resize(NumFoil);
  vector <TH2F*> hXFpYFp;
  vector <TH2F*> hXpFpXFp;
  hXpFpXFp.resize(NumFoil);
  hXFpYFp.resize(NumFoil);
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_YpYfpCut;
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_XpXfpCut;
  vector<vector<vector<TH1F*> > > hXs_DelCut_YpYfpCut;
  vector<vector<TH2F*> > hYsXs_DelCut;
  vector<vector<TH2F*> > hYpFpYFp_DelCut;
  vector<vector<TH2F*> > hXpFpXFp_DelCut;
  cout << " setup DelCut 2d" << endl;
  hYsXs_DelCut.resize(NumFoil);
  hYsXs_DelCut_YpYfpCut.resize(NumFoil);
  hYsXs_DelCut_XpXfpCut.resize(NumFoil);
  hXs_DelCut_YpYfpCut.resize(NumFoil);
  hYpFpYFp_DelCut.resize(NumFoil);
  hXpFpXFp_DelCut.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    hYsXs_DelCut[nf].resize(ndelcut);
    hYsXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYsXs_DelCut_XpXfpCut[nf].resize(ndelcut);
    hXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYpFpYFp_DelCut[nf].resize(ndelcut);
    hXpFpXFp_DelCut[nf].resize(ndelcut);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut_YpYfpCut[nf][nd].resize(11);
      hYsXs_DelCut_XpXfpCut[nf][nd].resize(11);
      hXs_DelCut_YpYfpCut[nf][nd].resize(11);
    }
  }
  cout << " finish setup DelCut 2d" << endl;
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    hYsDelta[nc] = new TH2F(Form("hYsDelta_Foil_%d",nc),Form("Run %d Foil %d; Ys ; Delta",nc,nrun),100,-12,12,50,-15.,25.);
    HList.Add(hYsDelta[nc]);
    hXsDelta[nc] = new TH2F(Form("hXsDelta_Foil_%d",nc),Form("Run %d Foil %d; Xs ; Delta",nc,nrun),100,-15,15,50,-15.,25.);
    HList.Add(hXsDelta[nc]);
    hYpFpYFp[nc] = new TH2F(Form("hYpFpYFp_%d",nc),Form("Run %d Foil %d; Ypfp ; Yfp",nrun,nc),100,-.05,.05,100,-35.,35.);
    HList.Add(hYpFpYFp[nc]);
    hXpFpXFp[nc] = new TH2F(Form("hXpFpXFp_%d",nc),Form("MC Run %d Foil %4.1f; Xpfp ; Xfp",nrun,ztar_foil[nc]),100,-.1,.1,100,-50.,50.);
    HList.Add(hXpFpXFp[nc]);
    hXFpYFp[nc] = new TH2F(Form("hXFpYFp_%d",nc),Form("MC Run %d Foil %4.1f; Yfp ; Xfp",nrun,ztar_foil[nc]),100,-40.,40,100,-50.,50.);
    HList.Add(hXFpYFp[nc]);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut[nc][nd]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Ys ; Xs",nrun,nc,delcut[nd]),50,-12,12,100,-15.,15.);
      HList.Add(hYsXs_DelCut[nc][nd]);
      for  (Int_t ny=0;ny<11;ny++) {
	hYsXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d DelCut %3.1f Ys=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-12,12,100,-15.,15.);
	HList.Add(hYsXs_DelCut_YpYfpCut[nc][nd][ny]);
	hYsXs_DelCut_XpXfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_XFpCut_%d",nc,nd,ny),Form("Run %d Foil %d DelCut %3.1f Xs=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-12,12,100,-15.,15.);
	HList.Add(hYsXs_DelCut_XpXfpCut[nc][nd][ny]);
	hXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH1F(Form("hXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d DelCut %3.1f Ys=%d; Xs",nrun,nc,delcut[nd],ny),100,-15.,15.);
	HList.Add(hXs_DelCut_YpYfpCut[nc][nd][ny]);
      }
      hYpFpYFp_DelCut[nc][nd]  = new TH2F(Form("hYpFpYFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Ypfp ; Yfp",nrun,nc,delcut[nd]),75,-.05,.05,150,-35.,35.);
      HList.Add(hYpFpYFp_DelCut[nc][nd]);
      hXpFpXFp_DelCut[nc][nd]= new TH2F(Form("hXpFpXFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Xpfp ; Xfp",nrun,nc,delcut[nd]),150,-.1,.1,150,-40.,40.);
      HList.Add(hXpFpXFp_DelCut[nc][nd]);
    }
  }	  
  //
  // loop over entries
  Long64_t nentries = tsimc->GetEntries();
  cout << " start loop " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
  tsimc->GetEntry(i);
  if (i%50000==0) cout << " Entry = " << i << endl;
  hxbpm_tar->Fill(xbpm_tar);
  hybpm_tar->Fill(ybpm_tar);
  hZtarFrXa->Fill(fr_xa,reactz);
  hZtarFrYa->Fill(fr_ya,reactz);
    if (sumnpe > 6.) hetot->Fill(etracknorm);
    if (etracknorm>.8) hngsum->Fill(sumnpe);
    if (sumnpe > 6. && delta>-15 && delta<24) {
      if (delta>-10 && delta<24) hytar->Fill(ytar);
      if (delta>-10 && delta<24) hztar->Fill(reactz);
      hXptarDelta->Fill(xptar,delta);
      hYptarDelta->Fill(yptar,delta);
      hYtarDelta->Fill(ytar,delta);
      hYtarYptar->Fill(yptar,ytar);
      hYpFpYFp_all->Fill(ypfp,yfp);
      hXpFpXFp_all->Fill(xpfp,xfp);
      hYFpXFp_all->Fill(yfp,xfp); 
      hYtarYptar->Fill(yptar,ytar);
      hZtarDelta->Fill(reactz,delta);
      for  (UInt_t nc=0;nc<ytar_delta_cut.size();nc++) {
	if (ytar_delta_cut[nc]->IsInside(ytar,delta))	{ 
	  hYsDelta[nc]->Fill(ysieve,delta);
	  hXsDelta[nc]->Fill(xsieve,delta);
	  hYpFpYFp[nc]->Fill(ypfp,yfp);
	  for  (Int_t nd=0;nd<ndelcut;nd++) {
	    //if ( delta >=delcut[nd] && delta <delcut[nd+1]) {
	    if ( delta >=delcut[nd]-delwidth[nd] && delta <delcut[nd]+delwidth[nd]) {
	      hYsXs_DelCut[nc][nd]->Fill(ysieve,xsieve); 
	      hYpFpYFp_DelCut[nc][nd]->Fill(ypfp,yfp);
	      hXpFpXFp_DelCut[nc][nd]->Fill(xpfp,xfp);
	      Int_t f_ny=-1;
	      for  (UInt_t ny=0;ny<11;ny++) {
		if (CutYpFpYFpFlag && ypfp_yfp_cut[nc][nd][ny] && ypfp_yfp_cut[nc][nd][ny]->IsInside(ypfp,yfp)) {
		  hYsXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(ysieve,xsieve);
		  hXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(xsieve);
		  f_ny=ny;
		}
	      }
	      // temporary to clean up
	      if (f_ny !=-1) hXpFpXFp_DelCut[nc][nd]->Fill(xpfp,xfp);
	      hXpFpXFp_DelCut[nc][nd]->Fill(xpfp,xfp);
	      for  (UInt_t nx=0;nx<11;nx++) {
		if (f_ny != -1 && CutXpFpXFpFlag && xpfp_xfp_cut[nc][nd][nx] && xpfp_xfp_cut[nc][nd][nx]->IsInside(xpfp,xfp)) {
		  hYsXs_DelCut_XpXfpCut[nc][nd][nx]->Fill(ysieve,xsieve);
		}
	      }
	    }			     
	  }
	}
      }
    }
  }
  //
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  //
}
