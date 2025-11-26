#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
void fit_opt_matrix(TString fname="Optics_1814_1919_shms_newfit3_tree") {
  
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetPadLeftMargin(0.17);
  //
  string newcoeffsfilename="newfit.dat";
  //string oldcoeffsfilename="shms_newfit_april2020.dat";
  //string oldcoeffsfilename="mkj_fit_10776.dat";//shms_newfit_april2020.dat";
  //string oldcoeffsfilename="hsv_fit_10790_stripped.dat";//shms_newfit_april2020.dat";
  string oldcoeffsfilename="hsv_fit_10790_low_stripped.dat";//shms_newfit_april2020.dat";
  int nfit=0,npar,nfit_max=19328,npar_final=0,max_order=6,norder;
  Int_t MaxPerBin=200;
  Int_t MaxZtarPerBin=10000;
  //
  TH1F *hDelta = new TH1F("hDelta","Delta ",20,-10.,30.);
  TH1F *hDeltarecon = new TH1F("hDeltarecon","Delta Recon % ",40,-10,30);
  TH1F *hDeltadiff = new TH1F("hDeltadiff","Delta Diff % ",40,-1,1);
  TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon % ",40,-10,30);
  TH1F *hDeltanewdiff = new TH1F("hDeltanewdiff","Delta New Diff % ",40,-1,1);
  TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
  TH1F *hytarrecon = new TH1F("hytarrecon","ytar recon(cm)",70,-7.,7.);
  TH1F *hytarnew = new TH1F("hytarnew","ytar new(cm)",70,-7.,7.);
  TH1F *hytardiff = new TH1F("hytardiff","ytar diff(cm)",70,-5.,5.);
  TH1F *hytarnewdiff = new TH1F("hytarnewdiff","ytar new diff(cm)",70,-5.,5.);
  TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
  TH1F *hxptarrecon = new TH1F("hxptarrecon","xptar recon",100,-.12,.12);
  TH1F *hxptardiff = new TH1F("hxptardiff","xptar diff (mr)",100,-15,15);
  TH1F *hxptarnew = new TH1F("hxptarnew","xptar new recon",100,-.12,.12);
  TH1F *hxptarnewdiff = new TH1F("hxptarnewdiff","xptar new diff (mr)",100,-15,15);
  TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
  TH1F *hyptarrecon = new TH1F("hyptarrecon","yptar recon ",100,-.1,.1);
  TH1F *hyptardiff = new TH1F("hyptardiff","yptar diff (mr) ",100,-10,10);
  TH1F *hyptarnew = new TH1F("hyptarnew","yptar new recon ",100,-.1,.1);
  TH1F *hyptarnewdiff = new TH1F("hyptarnewdiff","yptar new diff (mr) ",100,-10,10);
  //
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
  
  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_fit;
  vector<double> yptarcoeffs_fit;
  vector<double> ytarcoeffs_fit;
  vector<double> deltacoeffs_fit;
  vector<int> xfpexpon_fit;
  vector<int> xpfpexpon_fit;
  vector<int> yfpexpon_fit;
  vector<int> ypfpexpon_fit;
  vector<int> xtarexpon_fit;

  vector<double> xptarcoeffs_xtar;
  vector<double> yptarcoeffs_xtar;
  vector<double> ytarcoeffs_xtar;
  vector<double> deltacoeffs_xtar;
  vector<int> xfpexpon_xtar;
  vector<int> xpfpexpon_xtar;
  vector<int> yfpexpon_xtar;
  vector<int> ypfpexpon_xtar;
  vector<int> xtarexpon_xtar;

  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  TString currentline;
  int num_recon_terms_old;
  int num_recon_terms_fit;
  int num_recon_terms_xtar;

  num_recon_terms_old = 0;
  num_recon_terms_fit = 0;
  num_recon_terms_xtar = 0;
  // add zero order term to fit
  xptarcoeffs_fit.push_back(0.0);
  ytarcoeffs_fit.push_back(0.0);
  yptarcoeffs_fit.push_back(0.0);
  deltacoeffs_fit.push_back(0.0);
  xfpexpon_fit.push_back(0);
  xpfpexpon_fit.push_back(0);
  yfpexpon_fit.push_back(0);
  ypfpexpon_fit.push_back(0);
  xtarexpon_fit.push_back(0);
  num_recon_terms_fit = 1;
  
  
  //
  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    
    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs_old.push_back(sc1.Atof());
    ytarcoeffs_old.push_back(sc2.Atof());
    yptarcoeffs_old.push_back(sc3.Atof());
    deltacoeffs_old.push_back(sc4.Atof());
    
    int expontemp[5];

    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

  

    xfpexpon_old.push_back(expontemp[0]);
    xpfpexpon_old.push_back(expontemp[1]);
    yfpexpon_old.push_back(expontemp[2]);
    ypfpexpon_old.push_back(expontemp[3]);
    xtarexpon_old.push_back(expontemp[4]);

   
    
    num_recon_terms_old++;
    norder= expontemp[0]+expontemp[1]+expontemp[2]+expontemp[3]+expontemp[4];
    if (expontemp[4]==0) {
    xptarcoeffs_fit.push_back(sc1.Atof());
    ytarcoeffs_fit.push_back(sc2.Atof());
    yptarcoeffs_fit.push_back(sc3.Atof());
    deltacoeffs_fit.push_back(sc4.Atof());
    xfpexpon_fit.push_back(expontemp[0]);
    xpfpexpon_fit.push_back(expontemp[1]);
    yfpexpon_fit.push_back(expontemp[2]);
    ypfpexpon_fit.push_back(expontemp[3]);
    xtarexpon_fit.push_back(expontemp[4]);
    //cout << num_recon_terms_fit << " " <<  xptarcoeffs_fit[num_recon_terms_fit] << " " << ytarcoeffs_fit[num_recon_terms_fit] << " " <<  yptarcoeffs_fit[num_recon_terms_fit] << " " << deltacoeffs_fit[num_recon_terms_fit] << " " << xfpexpon_fit[num_recon_terms_fit] << " " << xpfpexpon_fit[num_recon_terms_fit] << " " << yfpexpon_fit[num_recon_terms_fit] << " " << ypfpexpon_fit[num_recon_terms_fit] << " " << xtarexpon_fit[num_recon_terms_fit] << " " << endl;
    num_recon_terms_fit++;
    } else {
    xptarcoeffs_xtar.push_back(sc1.Atof());
    ytarcoeffs_xtar.push_back(sc2.Atof());
    yptarcoeffs_xtar.push_back(sc3.Atof());
    deltacoeffs_xtar.push_back(sc4.Atof());
    xfpexpon_xtar.push_back(expontemp[0]);
    xpfpexpon_xtar.push_back(expontemp[1]);
    yfpexpon_xtar.push_back(expontemp[2]);
    ypfpexpon_xtar.push_back(expontemp[3]);
    xtarexpon_xtar.push_back(expontemp[4]);
    num_recon_terms_xtar++;
    }
  }

  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  cout << "num recon terms in fit matrix = " << num_recon_terms_fit << endl;
  cout << "num recon terms in xtar matrix = " << num_recon_terms_xtar << endl;
  npar= num_recon_terms_fit ;
  //
  TVectorD b_ytar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
  //
  TString inputroot;
  inputroot = "hist/"+fname+".root";
  cout << " INfile = " << inputroot << endl;
  TFile *fsimc = new TFile(inputroot);
  TTree *FitTree = (TTree*)fsimc->Get("TFit");
  //Declaration of leaves types
  Double_t  ys,xtar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT;
  Double_t xfp,xpfp,yfp,ypfp,ysieveT,ysieve;
  FitTree->SetBranchAddress("ys",&ysieve);
  FitTree->SetBranchAddress("ysT",&ysieveT);
  FitTree->SetBranchAddress("xtar",&xtar);
  FitTree->SetBranchAddress("xptar",&xptar);
  FitTree->SetBranchAddress("yptar",&yptar);
  FitTree->SetBranchAddress("ytar",&ytar);
  FitTree->SetBranchAddress("xptarT",&xptarT);
  FitTree->SetBranchAddress("yptarT",&yptarT);
  FitTree->SetBranchAddress("ytarT",&ytarT);
  FitTree->SetBranchAddress("ztarT",&ztarT);
  FitTree->SetBranchAddress("delta",&delta);
  FitTree->SetBranchAddress("xpfp",&xpfp);
  FitTree->SetBranchAddress("ypfp",&ypfp);
  FitTree->SetBranchAddress("xfp",&xfp);
  FitTree->SetBranchAddress("yfp",&yfp);
  //
  vector<Int_t > Ztar_Cnts;
  vector<vector<vector<Int_t> > > Ztar_Ys_Delta_Cnts;
  const Int_t nfoils=3;
  vector <Double_t> ztar_foil;
  //ztar_foil.push_back(20.);
  ztar_foil.push_back(6.67*2);
  ztar_foil.push_back(0.0);
  ztar_foil.push_back(-20.0);
  //ztar_foil.push_back(-30.0);
  Ztar_Cnts.resize(nfoils);
  Ztar_Ys_Delta_Cnts.resize(nfoils);
  static const Int_t ndelcut=7;
  Double_t delcut[ndelcut]={-11.,-8.,-4.,0,4.,8.,17.5};
  Double_t delwidth[ndelcut]={1.0,2.0,2.0,2.0,2.0,2.0,7.5};
  const Int_t nysieve=11;
  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*1.64-1.64*5;
    ys_cent.push_back(pos);
  }
  for (Int_t nf=0;nf<nfoils;nf++) {
    Ztar_Ys_Delta_Cnts[nf].resize(ndelcut);
    for (Int_t nd=0;nd<ndelcut;nd++) {
      Ztar_Ys_Delta_Cnts[nf][nd].resize(nysieve);
    }
  }
  //
  for (Int_t nf=0;nf<nfoils;nf++) {
    for (Int_t nd=0;nd<ndelcut;nd++) {
      for (Int_t ny=0;ny<nysieve;ny++) {	
	Ztar_Ys_Delta_Cnts[nf][nd][ny]=0;
      }}}
  
  //
  Long64_t nentries = FitTree->GetEntries();
  for (int i = 0; i < nentries; i++) {
    FitTree->GetEntry(i);
    //
    if (TMath::Abs(delta) < 100. ) {
      Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
      Double_t etemp;
      for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar/100., xtarexpon_old[icoeffold] );
	deltatemp += deltacoeffs_old[icoeffold] * etemp;
	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
      } // for icoeffold loop
      hytarrecon->Fill(ytartemp*100.);
      hyptarrecon->Fill(yptartemp);
      hxptarrecon->Fill(xptartemp);
      if ( delta>-15. &&  delta<30. ) {
	//
	Int_t found_nf=-1;
	Int_t found_nd=-1;
	Int_t found_ny=-1;
	Bool_t good_bin=kFALSE;
	for (Int_t nf=0;nf<nfoils;nf++) {
	  if (abs(ztarT-ztar_foil[nf])<1) found_nf=nf;
	}
	for (Int_t nd=0;nd<ndelcut;nd++) {
	  if (delta >=delcut[nd]-delwidth[nd] && delta <delcut[nd]+delwidth[nd]) found_nd=nd;
	}
	for (Int_t ny=0;ny<nysieve;ny++) {	
	  if (abs(ysieveT-ys_cent[ny])<.5) found_ny=ny;
	}
	if (found_nf!=-1 &&found_nd!=-1 && found_ny!=-1)  {
	  good_bin=kTRUE;
	}
	if (good_bin && nfit < nfit_max && Ztar_Ys_Delta_Cnts[found_nf][found_nd][found_ny]< MaxPerBin && Ztar_Cnts[found_nf]< MaxZtarPerBin) {
	  //
          Double_t ytar_xtar = 0.0,yptar_xtar=0.0,xptar_xtar=0.0;
          for( int icoeff_xtar=0; icoeff_xtar<num_recon_terms_xtar; icoeff_xtar++ ){
	    etemp= 
	      pow( xfp / 100.0, xfpexpon_xtar[icoeff_xtar] ) * 
	      pow( yfp / 100.0, yfpexpon_xtar[icoeff_xtar] ) * 
	      pow( xpfp, xpfpexpon_xtar[icoeff_xtar] ) * 
	      pow( ypfp, ypfpexpon_xtar[icoeff_xtar] ) * 
	      pow( xtar/100., xtarexpon_xtar[icoeff_xtar] );
	    ytar_xtar += ytarcoeffs_xtar[icoeff_xtar] * etemp;
	    yptar_xtar += yptarcoeffs_xtar[icoeff_xtar] * etemp;
	    xptar_xtar += xptarcoeffs_xtar[icoeff_xtar] *etemp; 
	  }
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
	    etemp= 
	      pow( xfp / 100.0, xfpexpon_fit[icoeff_fit] ) * 
	      pow( yfp / 100.0, yfpexpon_fit[icoeff_fit] ) * 
	      pow( xpfp, xpfpexpon_fit[icoeff_fit] ) * 
	      pow( ypfp, ypfpexpon_fit[icoeff_fit] ) * 
	      pow( xtar/100., xtarexpon_fit[icoeff_fit] );
	    if (nfit < nfit_max ) {
              lambda[icoeff_fit][nfit] = etemp;
	      b_xptar[icoeff_fit] += (xptarT-xptar_xtar) * etemp;
	      b_yptar[icoeff_fit] += (yptarT-yptar_xtar) * etemp;
	      b_ytar[icoeff_fit] += (ytarT-ytar_xtar*100) /100.0 * etemp;
	    }
	  } // for icoeff_fit loop
	  hytar->Fill(ytar);
	  hyptar->Fill(yptar);
	  hxptar->Fill(xptar);
	  hytardiff->Fill(ytar-ytarT);
	  hyptardiff->Fill(1000.*(yptar-yptarT));
	  hxptardiff->Fill(1000.*(xptar-xptarT));
	  Ztar_Cnts[found_nf]++;
	  Ztar_Ys_Delta_Cnts[found_nf][found_nd][found_ny]++;
	  nfit++;
	  xfptrue.push_back( xfp );
	  yfptrue.push_back( yfp );
	  xpfptrue.push_back( xpfp );
	  ypfptrue.push_back( ypfp );
	  xtartrue.push_back( xtar );
	  xptartrue.push_back( xptarT );
	  ytartrue.push_back( ytarT  );
	  yptartrue.push_back( yptarT  );
	}
      }
    }
  }
  //
  for (Int_t nf=0;nf<nfoils;nf++) cout << " counts foil " << nf << " : " << Ztar_Cnts[nf] << endl;
  if (nfit < nfit_max) {
    cout << " nfit < nfit_max set nfit_max = " << nfit << endl;
    return;
  }
  //
  for (Int_t nf=0;nf<nfoils;nf++) {
    cout << " ztar = " << ztar_foil[nf] << endl;
    for (Int_t nd=0;nd<ndelcut;nd++) {
      cout << " Ndelta = " << delcut[nd] << endl;       
      for (Int_t ny=0;ny<nysieve;ny++) {	
	cout <<  Ztar_Ys_Delta_Cnts[nf][nd][ny] << " " ;
      }
      cout << endl;
    }}
  
  //
  //
  cout << " number to fit = " << nfit << " max = " << nfit_max << endl;
  for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
  }
  for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 5000 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
      }
    }
  }
  //
  TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_ytar );
  cout << "ytar solution ok = " << ok << endl;
  //b_ytar.Print();
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  //b_yptar.Print();
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  //b_xptar.Print();
  // calculate target quantities with new fit parameter
  for( int ifit=0; ifit<nfit; ifit++){
    Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
    Double_t etemp;
    for( int ipar=0; ipar<npar; ipar++){
      etemp=lambda[ipar][ifit];
      ytarnew += b_ytar[ipar] * etemp;
      yptarnew += b_yptar[ipar] * etemp;
      xptarnew += b_xptar[ipar] *etemp;        
    }
    Double_t ytar_xtar = 0.0,yptar_xtar=0.0,xptar_xtar=0.0;
    for( int icoeff_xtar=0; icoeff_xtar<num_recon_terms_xtar; icoeff_xtar++ ){
      etemp= 
	pow( xfptrue.at(ifit) / 100.0, xfpexpon_xtar[icoeff_xtar] ) * 
	pow( yfptrue.at(ifit) / 100.0, yfpexpon_xtar[icoeff_xtar] ) * 
	pow( xpfptrue.at(ifit), xpfpexpon_xtar[icoeff_xtar] ) * 
	pow( ypfptrue.at(ifit), ypfpexpon_xtar[icoeff_xtar] ) * 
	pow( xtartrue.at(ifit)/100., xtarexpon_xtar[icoeff_xtar] );
      ytar_xtar += ytarcoeffs_xtar[icoeff_xtar] * etemp;
      yptar_xtar += yptarcoeffs_xtar[icoeff_xtar] * etemp;
      xptar_xtar += xptarcoeffs_xtar[icoeff_xtar] *etemp; 
    }
    hytarnew->Fill((ytarnew+ytar_xtar)*100.);
    hyptarnew->Fill(yptarnew+yptar_xtar);
    hxptarnew->Fill(xptarnew);
    hytarnewdiff->Fill((ytarnew+ytar_xtar)*100.-ytartrue.at(ifit));
    hyptarnewdiff->Fill(1000*(yptarnew+yptar_xtar-yptartrue.at(ifit)));
    hxptarnewdiff->Fill(1000*(xptarnew+xptar_xtar-xptartrue.at(ifit)));
  }
  // write out coeff
  char coeffstring[100];
  Double_t tt;
  cout << "writing new coeffs file" << endl;
  newcoeffsfile << "! new fit to "+fname << endl;
  newcoeffsfile << " ---------------------------------------------" << endl;
  for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
    newcoeffsfile << " ";
    //      tt=xptarcoeffs_fit[icoeff_fit];
    tt=b_xptar[icoeff_fit] ;
    sprintf( coeffstring, "%16.9g", tt );
    newcoeffsfile << coeffstring; 
    //      newcoeffsfile << " ";
    sprintf( coeffstring, "%16.9g", b_ytar[icoeff_fit] );
    newcoeffsfile << coeffstring;
    sprintf( coeffstring, "%16.9g", b_yptar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    sprintf( coeffstring, "%16.9g", deltacoeffs_fit[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    newcoeffsfile << " ";
    newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_fit[icoeff_fit]; 
    newcoeffsfile << endl;
    
  }
  //
  for( int icoeff_fit=0; icoeff_fit<num_recon_terms_xtar; icoeff_fit++ ){
    newcoeffsfile << " ";
    //      tt=xptarcoeffs_fit[icoeff_fit];
    tt=xptarcoeffs_xtar[icoeff_fit] ;
    sprintf( coeffstring, "%16.9g", tt );
    newcoeffsfile << coeffstring; 
    //      newcoeffsfile << " ";
    sprintf( coeffstring, "%16.9g", ytarcoeffs_xtar[icoeff_fit] );
    newcoeffsfile << coeffstring;
    sprintf( coeffstring, "%16.9g", yptarcoeffs_xtar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    sprintf( coeffstring, "%16.9g", deltacoeffs_xtar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    newcoeffsfile << " ";
    newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_xtar[icoeff_fit]; 
    newcoeffsfile << endl;
    
  }
  //
  newcoeffsfile << " ---------------------------------------------" << endl;
  
  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
  //
  TCanvas *cdiff = new TCanvas("cdiff","Old matrix Diff target",800,800);
  cdiff->Divide(2,2);
  cdiff->cd(1);
  hytardiff->Draw();
  hytardiff->Fit("gaus");
  TF1 *fitcydiff=hytardiff->GetFunction("gaus");
  cdiff->cd(2);
  hyptardiff->Draw();
  hyptardiff->Fit("gaus");
  TF1 *fitcypdiff=hyptardiff->GetFunction("gaus");
  cdiff->cd(3);
  hxptardiff->Draw();
  hxptardiff->Fit("gaus");
  TF1 *fitcxpdiff=hxptardiff->GetFunction("gaus");
  cdiff->cd(4);
  //  hDeltadiff->Draw();
  //hDeltadiff->Fit("gaus");
  //TF1 *fitcdeldiff=hDeltadiff->GetFunction("gaus");
  //
  TCanvas *cnewdiff = new TCanvas("cnewdiff","Newfit diff target",800,800);
  cnewdiff->Divide(2,2);
  cnewdiff->cd(1);
  hytarnewdiff->Draw();
  hytarnewdiff->Fit("gaus");
  TF1 *fitcynewdiff=hytarnewdiff->GetFunction("gaus");
  cnewdiff->cd(2);
  hyptarnewdiff->Draw();
  hyptarnewdiff->Fit("gaus");
  TF1 *fitcypnewdiff=hyptarnewdiff->GetFunction("gaus");
  cnewdiff->cd(3);
  hxptarnewdiff->Draw();
  hxptarnewdiff->Fit("gaus");
  TF1 *fitcxpnewdiff=hxptarnewdiff->GetFunction("gaus");
  cnewdiff->cd(4);
  //   hDeltanewdiff->Draw();
  //hDeltanewdiff->Fit("gaus");
  //TF1 *fitcdelnewdiff=hDeltanewdiff->GetFunction("gaus");
  //
}
