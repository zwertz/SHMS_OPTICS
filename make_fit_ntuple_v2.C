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

void make_fit_ntuple_v2(Int_t nrun=1814,Int_t FileID=-2){
  Double_t yMP = -0.03;//0.0;
  Double_t xMP = -0.126;//0.0;
  Double_t xptarOffset = 0.0;//-0.002221854212;//offset used in replay as xptar = xptar + offset
  
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
 TString OpticsFile = "list_of_optics_run.dat";
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
   if (nrun==16023){inputroot=Form("ROOTfiles/cafe_replay_optics_%s_%d.root",OpticsID.Data(),FileID);}
  else if (nrun==16024){inputroot=Form("ROOTfiles/cafe_replay_optics_comb_6p8deg_%d.root",FileID);}
  else if (nrun==16031){inputroot=Form("ROOTfiles/cafe_replay_optics_%s_%d.root",OpticsID.Data(),FileID);}
  else if (nrun==16029){inputroot=Form("ROOTfiles/cafe_replay_prod_%s_%d.root",OpticsID.Data(),FileID);}
  else if (nrun==16033){inputroot=Form("ROOTfiles/cafe_replay_optics_%s_%d.root",OpticsID.Data(),FileID);}
  else {inputroot=Form("ROOTfiles/cafe_replay_optics_comb_8p3deg_%d.root",FileID);}
   outputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
   TH1F *hxbpm_tar = new TH1F("hxbpm_tar",Form("Run %d ; Xbpm_tar ; Counts",nrun),100,-2.,2.);
   TH1F *hybpm_tar = new TH1F("hybpm_tar",Form("Run %d ; Ybpm_tar ; Counts",nrun),100,-2.,2.);
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
   tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
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
 // tsimc->SetBranchAddress("P.extcor.ysieve",&ysieve);
 Double_t  xsieve;
 //tsimc->SetBranchAddress("P.extcor.xsieve",&xsieve);
 Double_t  xbpm_tar;
   tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&xbpm_tar);
 Double_t  ybpm_tar;
   tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&ybpm_tar);
    Double_t frx;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&frx);
    Double_t fry;
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&fry);
   
   Double_t xptarT,ytarT,yptarT,ysieveT,xsieveT,ztarT,ztar,xtarT;
   //
   TFile hroot(outputroot,"recreate");
   TTree *otree = new TTree("TFit","FitTree");
   otree->Branch("ys",&ysieve);
   otree->Branch("ysT",&ysieveT);
   otree->Branch("xs",&xsieve);
   otree->Branch("xsT",&xsieveT);
   otree->Branch("ztarT",&ztarT);
   otree->Branch("xtar",&xtar);
   otree->Branch("xtarT",&xtarT);
   otree->Branch("ztar",&ztar);
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

 //parse the input matrix elements
   //string coeffsfilename="newfit_global_zbin_30deg.dat";
   //string coeffsfilename="newfit_10790_zbin.dat";
   string coeffsfilename="hsv_fit_global.dat";
  ifstream coeffsfile(coeffsfilename.c_str());
  TString currentline;
  int num_recon_terms=0;

  vector<double> xptarcoeffs;
  vector<double> yptarcoeffs;
  vector<double> ytarcoeffs;
  vector<double> deltacoeffs;
  vector<int> xfpexpon;
  vector<int> xpfpexpon;
  vector<int> yfpexpon;
  vector<int> ypfpexpon;
  vector<int> xtarexpon;
  
   while( currentline.ReadLine(coeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    
    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs.push_back(sc1.Atof());
    ytarcoeffs.push_back(sc2.Atof());
    yptarcoeffs.push_back(sc3.Atof());
    deltacoeffs.push_back(sc4.Atof());
    
    int expontemp[5];

    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

    xfpexpon.push_back(expontemp[0]);
    xpfpexpon.push_back(expontemp[1]);
    yfpexpon.push_back(expontemp[2]);
    ypfpexpon.push_back(expontemp[3]);
    xtarexpon.push_back(expontemp[4]);


    cout << num_recon_terms << " " <<  xptarcoeffs[num_recon_terms] << " " << ytarcoeffs[num_recon_terms] << " " <<  yptarcoeffs[num_recon_terms] << " " << deltacoeffs[num_recon_terms] << " " << xfpexpon[num_recon_terms] << " " << xpfpexpon[num_recon_terms] << " " << yfpexpon[num_recon_terms] << " " << ypfpexpon[num_recon_terms] << " " << xtarexpon[num_recon_terms] << " " << endl;

    num_recon_terms++;   
  }
   
   Double_t zdis_sieve = 253.;
   Double_t zdis_hbmagexit_coll=40.;
   Double_t hb_delta_yptar_coeff = 0.019+zdis_hbmagexit_coll*0.00052;
   Double_t hb_delta_yptar_coeff2 = 0.00019+zdis_hbmagexit_coll*0.0000052;
   
   
   Long64_t nentries = tsimc->GetEntries();
   for (int i = 0; i < nentries; i++) {
     tsimc->GetEntry(i);
     if (i%50000==0) cout << " Entry = " << i << endl;
     if (etracknorm>.8 && sumnpe > 6. && delta>-10 && delta<10) {
       hxbpm_tar->Fill(xbpm_tar);
       hybpm_tar->Fill(ybpm_tar);		  
     }} //
   Double_t xbeam = -hxbpm_tar->GetMean(); // horizontal beam in Hall coordinates
   Double_t ybeam = hybpm_tar->GetMean();
   cout << " xbeam = " << xbeam << " ybeam = " << ybeam << endl;
   cout << " start loop " << nentries << endl;
   
   
   
   // loop over entries
   cout << " start loop " << nentries << endl;
   CentAngle=CentAngle*3.14159/180.;
   for (int i = 0; i < nentries; i++) {
     tsimc->GetEntry(i);

     ////////////////////////////////////////////////////////////////////
      // Calculate corrections & recalculate ,,,track parameters
    Double_t xptar_save=0.,xptar_diff=10000., xtar_new=-reacty;
    Double_t x_tg = -reacty-xMP; // units of cm, beam position in spectrometer coordinate system

    Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
    Double_t etemp;
    for( int icoeff=0; icoeff<num_recon_terms; icoeff++ ){
      etemp= 
	pow( xfp / 100.0, xfpexpon[icoeff] ) * 
	pow( yfp / 100.0, yfpexpon[icoeff] ) * 
	pow( xpfp, xpfpexpon[icoeff] ) * 
	pow( ypfp, ypfpexpon[icoeff] ) * 
	pow( x_tg/100., xtarexpon[icoeff] );
      deltatemp += deltacoeffs[icoeff] * etemp;
      ytartemp += ytarcoeffs[icoeff] * etemp;
      yptartemp += yptarcoeffs[icoeff] * etemp;
      xptartemp += xptarcoeffs[icoeff] *etemp; 
    } // for icoeffold loop

    xptar_save = xptartemp;
    Int_t niter=0;
    while ( xptar_diff > 2 && niter < 5) {
      xtar_new = x_tg - xptartemp*reactz*cos(CentAngle); //units of cm
      ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
      for( int icoeff=0; icoeff<num_recon_terms; icoeff++ ){
	etemp= 
	  pow( xfp / 100.0, xfpexpon[icoeff] ) * 
	  pow( yfp / 100.0, yfpexpon[icoeff] ) * 
	  pow( xpfp, xpfpexpon[icoeff] ) * 
	  pow( ypfp, ypfpexpon[icoeff] ) * 
	  pow( xtar_new/100., xtarexpon[icoeff] );
	deltatemp += deltacoeffs[icoeff] * etemp;
	ytartemp += ytarcoeffs[icoeff] * etemp;
	yptartemp += yptarcoeffs[icoeff] * etemp;
	xptartemp += xptarcoeffs[icoeff] *etemp; 
      } // for icoeffold loop
      xptar_diff = abs(xptartemp-xptar_save)*1000;
      xptar_save = xptartemp;
      niter++;
    }
    ///////////////////////////////////////////////////////////////////
    xsieve = xtar_new + xptartemp*253;
    Double_t delta_per = deltatemp*100.0;
    ysieve = ytartemp*100+yptartemp*253.-(0.019+40.*.01*0.052)*delta_per+(0.00019+40*.01*.00052)*delta_per*delta_per;
    ytartemp *=100;
    Double_t ztarg=(ytartemp-yMP-0*(cos(CentAngle)-yptartemp*sin(CentAngle)))/(-sin(CentAngle)-cos(CentAngle)*yptartemp);//mine
     ////////////////////////////////////////////////////////////////////
     
     if (i%50000==0) cout << " Entry = " << i << endl;
     if (etracknorm>.8 && sumnpe > 6. && delta>-15 && delta<24) {
       Int_t nf_found=-1, nd_found=-1,ny_found=-1,nx_found=-1;
       for  (UInt_t nf=0;nf<ytar_delta_cut.size();nf++) {
	 if (ytar_delta_cut[nf]->IsInside(ytartemp,delta_per)) nf_found=nf;
       } 
       for  (UInt_t nd=0;nd<ndelcut;nd++) {
	 if ( delta_per >=delcut[nd]-delwidth[nd] && delta_per <delcut[nd]+delwidth[nd])  nd_found=nd;
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



	 /*
	 Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
	 Double_t etemp;
	 for( int icoeff=0; icoeff<num_recon_terms; icoeff++ ){
	   etemp= 
	     pow( xfp / 100.0, xfpexpon[icoeff] ) * 
	     pow( yfp / 100.0, yfpexpon[icoeff] ) * 
	     pow( xpfp, xpfpexpon[icoeff] ) * 
	     pow( ypfp, ypfpexpon[icoeff] ) * 
	     pow( xtar/100., xtarexpon[icoeff] );
	   deltatemp += deltacoeffs[icoeff] * etemp;
	   ytartemp += ytarcoeffs[icoeff] * etemp;
	   yptartemp += yptarcoeffs[icoeff] * etemp;
	   xptartemp += xptarcoeffs[icoeff] *etemp; 
	 } // for icoeffold loop
	 
	 reactx = xbeam - (frx+0.01);
	 reacty = -ybeam - (fry+0.01);
	 if (CentAngle*180/3.14 > 25){reacty = -ybeam - (fry+0.005);}
	 
	 xtar = -reacty -(ztar_foil[nf_found]*TMath::Cos(CentAngle)+reactx*TMath::Sin(CentAngle))*xptartemp;
	 ysieve = ytartemp*100.0+yptartemp*zdis_sieve;
	 xsieve = xtar +xptartemp*zdis_sieve;
	 yptar = yptartemp;
	 xptar = xptartemp;
	 ytar = ytartemp*100.0;
	 delta=deltatemp*100.0;
	 */
	 
	 yptarT = (ys_cent[ny_found]+hb_delta_yptar_coeff*delta-hb_delta_yptar_coeff2*delta_per*delta_per+ztar_foil[nf_found]*TMath::Sin(CentAngle)-reactx*TMath::Cos(CentAngle)+yMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
	 ytarT = -ztar_foil[nf_found]*TMath::Sin(CentAngle)-yptarT*(ztar_foil[nf_found]*TMath::Cos(CentAngle)+reactx*TMath::Sin(CentAngle))+reactx*TMath::Cos(CentAngle)-yMP;
	 xptarT = (xs_cent[nx_found]+reacty+xMP)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle)-reactx*TMath::Sin(CentAngle));
	 xtarT = -reacty - xMP - xptarT*ztar_foil[nf_found]*TMath::Cos(CentAngle);
	 
	 ysieveT=ys_cent[ny_found];
	 xsieveT=xs_cent[nx_found];
	 ztarT=ztar_foil[nf_found];
	 ztar=reactz;
	 //ztar = (ytartemp-yMP-reactx*(cos(CentAngle)-yptartemp*sin(CentAngle)))/(-sin(CentAngle)-cos(CentAngle)*yptartemp);
	 ztar = (ytartemp-yMP-0*(cos(CentAngle)-yptartemp*sin(CentAngle)))/(-sin(CentAngle)-cos(CentAngle)*yptartemp);
	 //xptar = xptar - xptarOffset;
	 otree->Fill();
       }
     }
   }
   //
   otree->Write();
   //
}
