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

void make_fit_ntuple(Int_t nrun=1814,Int_t FileID=-2){
  Double_t yMP = 0.0;
  Double_t xMP = 0.0;
  
  Bool_t CutYtarFlag=kTRUE;
  Bool_t CutYpFpYFpFlag=kTRUE;
  Bool_t CutXpFpXFpFlag=kTRUE;
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
  Int_t ndelcut;
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
  for (Int_t nf=0;nf<NumFoil;nf++) {
    cout << nf << " foil = " << ztar_foil[nf] << endl;
  }
  vector <Double_t> ys_cent;
  vector <Double_t> xs_cent;
  for (Int_t nys=0;nys<11;nys++) {
    Double_t ypos=nys*1.64-1.64*5;
    ys_cent.push_back(ypos);
    Double_t xpos=nys*2.5-2.5*5;
    xs_cent.push_back(xpos);
  }
  //
 //
 //
   TString inputroot;
   TString outputroot;
   //inputroot = Form("ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_%d_newfit_global_zbin_allA1n.root",OpticsID.Data(),FileID);
   inputroot = Form("ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_%d_newfit_lad_shms.root",OpticsID.Data(),FileID);
   outputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
  //
 //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << " Cut file = " << YtarDeltaCutFile << endl;
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
      Int_t npt = tempg->GetN();
      //cout << "hYpFpYFp_cut = " << nf << " " << nd << " " << nc << " npts = " << npt << endl;
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	  ypfp_yfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  ypfp_yfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }
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
 Double_t  sumhgnpe;
 //   tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
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
   Double_t xptarT,ytarT,yptarT,ysieveT,xsieveT,ztarT,ztar,xtarT;
   //
 TFile hroot(outputroot,"recreate");
 TTree *otree = new TTree("TFit","FitTree");
 otree->Branch("ys",&ysieve);
 otree->Branch("ysT",&ysieveT);
 otree->Branch("xs",&xsieve);
 otree->Branch("xsT",&xsieveT);
 otree->Branch("ztarT",&ztarT); 
 otree->Branch("ztar",&ztar);
 otree->Branch("xtarT",&xtarT);
 otree->Branch("xtar",&xtar);
 otree->Branch("xptar",&xptar);
 otree->Branch("yptar",&yptar);
 otree->Branch("ytar",&ytar);
 otree->Branch("xptarT",&xptarT);
 otree->Branch("yptarT",&yptarT);
 otree->Branch("ytarT",&ytarT);
 otree->Branch("delta",&delta);
 otree->Branch("xpfp",&xpfp);
 otree->Branch("ypfp",&ypfp);
 otree->Branch("xfp",&xfp);
 otree->Branch("yfp",&yfp);
 otree->Branch("reactxcalc",&reactx);
 otree->Branch("reactycalc",&reacty);
   //
	Double_t zdis_sieve = 253.;
        Double_t zdis_hbmagexit_coll=40.;
	Double_t hb_delta_yptar_coeff = 0.019+zdis_hbmagexit_coll*0.00052;
	Double_t hb_delta_yptar_coeff2 = 0.00019+zdis_hbmagexit_coll*0.0000052;

// loop over entries
	Long64_t nentries = tsimc->GetEntries();
 cout << " start loop " << nentries << endl;
 CentAngle=CentAngle*3.14159/180.;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (etracknorm>.8 && sumnpe > 6. && delta>-15 && delta<24) {
		  Int_t nf_found=-1, nd_found=-1,ny_found=-1,nx_found=-1;
	         for  (UInt_t nf=0;nf<ytar_delta_cut.size();nf++) {
		   if (ytar_delta_cut[nf]->IsInside(ytar,delta)) nf_found=nf;
		 } 
                 for  (UInt_t nd=0;nd<ndelcut;nd++) {
		   if ( delta >=delcut[nd]-delwidth[nd] && delta <delcut[nd]+delwidth[nd])  nd_found=nd;
		     }
		 if (nf_found!=-1 && nd_found!=-1) {
                for  (UInt_t ny=0;ny<11;ny++) {
		  if (ypfp_yfp_cut[nf_found][nd_found][ny] && ypfp_yfp_cut[nf_found][nd_found][ny]->IsInside(ypfp,yfp)) ny_found=ny;
		}
                for  (UInt_t nx=0;nx<11;nx++) {
		  if (xpfp_xfp_cut[nf_found][nd_found][nx] && xpfp_xfp_cut[nf_found][nd_found][nx]->IsInside(xpfp,xfp)) nx_found=nx;
		}
		 }
		if (nf_found !=-1 && nd_found!=-1 && ny_found!=-1 && nx_found!=-1) {
		  //yptarT = (ys_cent[ny_found]+hb_delta_yptar_coeff*delta-hb_delta_yptar_coeff2*delta*delta+ztar_foil[nf_found]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
		  //ytarT = -ztar_foil[nf_found]*(TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle));
		  //xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
		  yptarT = (ys_cent[ny_found]+hb_delta_yptar_coeff*delta-hb_delta_yptar_coeff2*delta*delta+ztar_foil[nf_found]*TMath::Sin(CentAngle)-reactx*TMath::Cos(CentAngle)+yMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
		  ytarT = -ztar_foil[nf_found]*TMath::Sin(CentAngle)-yptarT*(ztar_foil[nf_found]*TMath::Cos(CentAngle)+reactx*TMath::Sin(CentAngle))+reactx*TMath::Cos(CentAngle)-yMP;
		  xptarT = (xs_cent[nx_found]+reacty+xMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
		  xtarT = -reacty - xMP - xptarT*ztar_foil[nf_found]*TMath::Cos(CentAngle);
		  
		  ysieveT=ys_cent[ny_found];
		  xsieveT=xs_cent[nx_found];
		  ztarT=ztar_foil[nf_found];
		  ztar=reactz;
		  otree->Fill();
		}
		}
	}
	//
	otree->Write();
//
}
