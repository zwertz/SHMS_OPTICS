void plotOpticsCompare(const char *nrun1,const char *nrun2, double foilPOS){

  //read the input file
  //TFile *f = new TFile(Form("/net/cdaq/cdaql3data/cdaq/hallc-online/ROOTfiles/coin_replay_production_%i_-1.root",nrun));
  //TFile *f1 = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_global_zbin_allA1n.root",nrun1));
  TFile *f1 = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_lad_shms.root",nrun2));
  TFile *f2 = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_lad_shms_v2.root",nrun2));
  //TFile *f2 = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_lad_shms_v3.root",nrun2));
  //TFile *f2 = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_lad_shms_v4.root",nrun2));
  TTree *tt1 = (TTree*)f1->Get("T");
  TTree *tt2 = (TTree*)f2->Get("T");

  //TString out_string1("newfit_global_zbin_allA1n");
  TString out_string1("newfit_lad_shms");
  TString out_string2("lad_shms_v2");
  
 //foil cut
  double foilmin = foilPOS-2.5;
  double foilmax = foilPOS+2.5;
  cout << foilmin << " is lower foilcut" << endl;
  cout << foilmax << " is upper foilcut" << endl;

  TString foil_pos(Form("Foil position: %0.1f",foilPOS));
  TString all_foil("All Foils");
  
  //here's the cut
  TCut cut = Form("P.gtr.dp<20&&P.gtr.dp>-15&&P.cal.etracknorm>0.8&&P.ngcer.npeSum>5&&P.react.z>%f&&P.react.z<%f",foilmin,foilmax);
  TCut cut_nom = "P.gtr.dp<20&&P.gtr.dp>-15&&P.cal.etracknorm>0.8&&P.ngcer.npeSum>5";
 TCut cutCentral = "abs(P.gtr.x+P.gtr.th*253.0)<1&&abs((-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph))<1";

  //make the output file
  TCanvas *canvas = new TCanvas("canvas","plots",1200,800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);

  TCanvas *canvas2 = new TCanvas("canvas2","plots2",1200,800);
  canvas2->SetFillColor(0);
  canvas2->SetBorderMode(0);
  canvas2->SetBorderSize(0);
  canvas2->SetFrameFillColor(0);
  canvas2->SetFrameBorderMode(0);

  std::string pdf_file_name= Form("output_plots_%s_%s_%s_%s_%1.f_sieveComparison.pdf",nrun1,nrun2,out_string1.Data(),out_string2.Data(),foilPOS);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("output_plots_%s_%s_%s_%s_%1.f_sieveComparison.root",nrun1,nrun2,out_string1.Data(),out_string2.Data(),foilPOS),"RECREATE");

  TString sieve_string1(Form("ysieve vs xsieve run1 %s;ySieve [cm];xSieve[cm]",all_foil.Data()));
  TString sieve_string2(Form("ysieve vs xsieve run1 %s;ySieve [cm];xSieve[cm]",foil_pos.Data()));

  TString sieve_string3(Form("ysieve vs xsieve run2 %s;ySieve [cm];xSieve[cm]",all_foil.Data()));
  TString sieve_string4(Form("ysieve vs xsieve run2 %s;ySieve [cm];xSieve[cm]",foil_pos.Data()));
  //make plots
  TH1F *h_z1_all = new TH1F("h_z1_all",";P.react.z [cm]",100,-14,14);
  TH1F *h_xsieve1_all = new TH1F("h_xsieve1_all",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve1_all = new TH1F("h_ysieve1_all",";ySieve[cm]",250,-9.0,9.0);
  TH2F *h2_sieve1_all = new TH2F("h2_sieve1_all",sieve_string1.Data(),250,-9.0,9.0,250,-14.0,14.0);
  TH1F *h_Vy1_all = new TH1F("h_Vy1_all",";yTar [cm];",100,-5,5);
  
  
  TH1F *h_z1 = new TH1F("h_z1",";P.react.z [cm]",100,-14,14);
  TH1F *h_xsieve1 = new TH1F("h_xsieve1",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve1 = new TH1F("h_ysieve1",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_xsieve_overlay1 = new TH1F("h_xsieve_overlay1",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve_overlay1 = new TH1F("h_ysieve_overlay1",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_Vy1 = new TH1F("h_Vy1",";yTar [cm];",100,-5,5);
  TH2F *h2_ypVy1 = new TH2F("h2_ypVy1","yTar vs ypTar run1;yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp1 = new TH2F("h2_yfpVxfp1","xfp vs yfp run1;xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz1 = new TH2F("h2_dpVz1","delta vs z run1;zVertex [cm];delta",100,-14,14,100,-15,20);
  TH2F *h2_ypVz1 = new TH2F("h2_ypVz1","z vs ypTar run1;zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve1 = new TH2F("h2_sieve1","ysieve vs xsieve run1;ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_sieve_overlay1 = new TH2F("h2_sieve_overlay1",sieve_string2.Data(),250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_xpVd1 = new TH2F("h2_xpVd1","delta vs xpfp run1;delta;xpfp",100,-10,20,100,-0.15,0.15);
  
  //plots with central hole only
  TH1F *h_z_c1 = new TH1F("h_z_c1","central sieve hole-run1;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy_c1 = new TH2F("h2_ypVy_c1","central sieve hole-run1;yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp_c1 = new TH2F("h2_yfpVxfp_c1","central sieve hole-run1;yfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz_c1 = new TH2F("h2_dpVz_c1","central sieve hole-run1;zVertex [cm];delta",100,-14,14,100,-10,20);
  TH2F *h2_ypVz_c1 = new TH2F("h2_ypVz_c1","central sieve hole-run1;zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve_c1 = new TH2F("h2_sieve_c1","central sieve hole-run1;ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
  TH2F *h2_xpVd_c1 = new TH2F("h2_xpVd_c1","central sieve hole-run1;delta;xpfp",100,-15,20,100,-0.15,0.15);

  //plots for other run
  TH1F *h_z2_all = new TH1F("h_z2_all","run2;P.react.z [cm]",100,-14,14);
  TH1F *h_xsieve2_all = new TH1F("h_xsieve2_all",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve2_all = new TH1F("h_ysieve2_all",";ySieve[cm]",250,-9.0,9.0);
  TH2F *h2_sieve2_all = new TH2F("h2_sieve2_all",sieve_string3.Data(),250,-9.0,9.0,250,-14.0,14.0);
  TH1F *h_Vy2_all = new TH1F("h_Vy2_all",";yTar [cm];",100,-5,5);
  
  TH1F *h_z2 = new TH1F("h_z2",";P.react.z [cm]",100,-14,14);
  TH1F *h_xsieve2 = new TH1F("h_xsieve2",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve2 = new TH1F("h_ysieve2",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_xsieve_overlay2 = new TH1F("h_xsieve_overlay2",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve_overlay2 = new TH1F("h_ysieve_overlay2",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_Vy2 = new TH1F("h_Vy2",";yTar [cm];",100,-5,5);
  TH2F *h2_ypVy2 = new TH2F("h2_ypVy2","yTar vs ypTar (run2);yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp2 = new TH2F("h2_yfpVxfp2","xfp vs yfp (run2);xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz2 = new TH2F("h2_dpVz2","z vs delta (run2);zVertex [cm];delta",100,-14,14,100,-15,20);
  TH2F *h2_ypVz2 = new TH2F("h2_ypVz2","z vs ypTar (run2);zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve2 = new TH2F("h2_sieve2","ysieve vs xsieve (run2);ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_sieve_overlay2 = new TH2F("h2_sieve_overlay2",sieve_string4.Data(),250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_xpVd2 = new TH2F("h2_xpVd2","delta vs xpfp (run2);delta;xpfp",100,-10,20,100,-0.15,0.15);
  
  //plots with central hole only
  TH1F *h_z_c2 = new TH1F("h_z_c2","central sieve hole-sk;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy_c2 = new TH2F("h2_ypVy_c2","central sieve hole-sk;yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp_c2 = new TH2F("h2_yfpVxfp_c2","central sieve hole-sk;yfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz_c2 = new TH2F("h2_dpVz_c2","central sieve hole-sk;zVertex [cm];delta",100,-14,14,100,-10,20);
  TH2F *h2_ypVz_c2 = new TH2F("h2_ypVz_c2","central sieve hole-sk;zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve_c2 = new TH2F("h2_sieve_c2","central sieve hole-sk;ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
  TH2F *h2_xpVd_c2 = new TH2F("h2_xpVd_c2","central sieve hole-sk;delta;xpfp",100,-15,20,100,-0.15,0.15);
  
  
  TH2F *h2_ypVzSlice[8];


  //plot this stuff
  tt1->Draw("P.react.z>>h_z1_all",cut_nom);
  tt1->Draw("P.extcor.xsieve>>h_xsieve1_all",cut_nom);
  tt1->Draw("P.extcor.ysieve>>h_ysieve1_all",cut_nom);
  tt1->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve1_all",cut_nom);
  tt1->Draw("P.gtr.y>>h_Vy1_all",cut_nom);
  
  tt1->Draw("P.gtr.ph:P.react.z>>h2_ypVz1",cut);
  tt1->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp1",cut);
  tt1->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy1",cut);
  tt1->Draw("P.react.z>>h_z1",cut);
  tt1->Draw("P.extcor.xsieve>>h_xsieve1",cut);
  tt1->Draw("P.extcor.ysieve>>h_ysieve1",cut);
  tt1->Draw("P.extcor.xsieve>>h_xsieve_overlay1",cut);
  tt1->Draw("P.extcor.ysieve>>h_ysieve_overlay1",cut);
  tt1->Draw("P.gtr.y>>h_Vy1",cut);

  tt1->Draw("P.gtr.dp:P.react.z>>h2_dpVz1",cut);
//tt1->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve1",cut);
  tt1->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve1",cut);
  tt1->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_overlay1",cut);
  tt1->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd1",cut);

  tt1->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c1",cut && cutCentral);
  tt1->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c1",cut && cutCentral);
  tt1->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c1",cut && cutCentral);
  tt1->Draw("P.react.z>>h_z_c1",cut && cutCentral);
  tt1->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c1",cut && cutCentral);
//tt1->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c1",cut && cutCentral);
  tt1->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_c1",cut && cutCentral);
  tt1->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c1",cut && cutCentral);

  tt2->Draw("P.react.z>>h_z2_all",cut_nom);
  tt2->Draw("P.extcor.xsieve>>h_xsieve2_all",cut_nom);
  tt2->Draw("P.extcor.ysieve>>h_ysieve2_all",cut_nom);
  tt2->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve2_all",cut_nom);
  tt2->Draw("P.gtr.y>>h_Vy2_all",cut_nom);
  
  tt2->Draw("P.gtr.ph:P.react.z>>h2_ypVz2",cut);
  tt2->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp2",cut);
  tt2->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy2",cut);
  tt2->Draw("P.react.z>>h_z2",cut);
  tt2->Draw("P.extcor.xsieve>>h_xsieve2",cut);
  tt2->Draw("P.extcor.ysieve>>h_ysieve2",cut);
  tt2->Draw("P.extcor.xsieve>>h_xsieve_overlay2",cut);
  tt2->Draw("P.extcor.ysieve>>h_ysieve_overlay2",cut);
  tt2->Draw("P.gtr.y>>h_Vy2",cut);

  tt2->Draw("P.gtr.dp:P.react.z>>h2_dpVz2",cut);
  //tt2->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve2",cut);
  tt2->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve2",cut);
  tt2->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_overlay2",cut);
  tt2->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd2",cut);

  tt2->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c2",cut && cutCentral);
  tt2->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c2",cut && cutCentral);
  tt2->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c2",cut && cutCentral);
  tt2->Draw("P.react.z>>h_z_c2",cut && cutCentral);
  tt2->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c2",cut && cutCentral);
//tt2->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c2",cut && cutCentral);
  tt2->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_c2",cut && cutCentral);
  tt2->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c2",cut && cutCentral);

  TH1F *h_z_res = new TH1F("h_z_res","Residual P.react.z",100,-14,14);
  h_z_res->Add(h_z1_all,h_z2_all,-1);

  TH1F *h_xsieve_res = new TH1F("h_xsieve_res","Residual xsieve foil",250, -14,14);
  h_xsieve_res->Add(h_xsieve1,h_xsieve2,-1);

  TH1F *h_ysieve_res = new TH1F("h_ysieve_res","Residual ysieve foil",250, -9,9);
  h_ysieve_res->Add(h_ysieve1,h_ysieve2,-1);

  TH1F *h_xsieve_all_res = new TH1F("h_xsieve_all_res","Residual xsieve all foils",250, -14,14);
  h_xsieve_all_res->Add(h_xsieve1_all,h_xsieve2_all,-1);

  TH1F *h_ysieve_all_res = new TH1F("h_ysieve_all_res","Residual ysieve all foils",250, -9,9);
  h_ysieve_all_res->Add(h_ysieve1_all,h_ysieve2_all,-1);
  
  for (int ii=0; ii<8; ii++){
    h2_ypVzSlice[ii] = new TH2F(Form("h2_ypVzSlice_%i",ii),Form("xfp=%i cm +/- 0.5cm;zVertex [cm]; ypTar",ii+1),100,-15,15,100,-0.05,0.05);
    TCut slice = Form("abs(P.dc.x_fp - (%i+1))<0.5",ii);
    tt1->Draw(Form("P.gtr.ph:P.react.z>>h2_ypVzSlice_%i",ii),cut && slice);
  }

  //Create lines where sieve holes should be
  //For 2D Sieve
  TLine* ys_line[11];
 TText* ys_text[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*1.64-1.64*5;
    ys_line[nys]= new TLine(pos,-14.,pos,14.);
    ys_text[nys]= new TText(pos,-17,Form("%d",nys));
     ys_text[nys]->SetTextColor(2);
     ys_line[nys]->SetLineColor(2);
     ys_line[nys]->SetLineWidth(1);
  }
 TLine* xs_line[11];
 TText* xs_text[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*2.5-2.5*5;
    xs_line[nys]= new TLine(-9.,pos,9.,pos);
    xs_text[nys]= new TText(-10.5,pos,Form("%d",nys));
     xs_text[nys]->SetTextColor(2);
     xs_line[nys]->SetLineColor(2);
     xs_line[nys]->SetLineWidth(1);
  }

  //For 1D sieve plots
  TLine* ys_line_2[11];
  TText* ys_text_2[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*1.64-1.64*5;
    ys_line_2[nys]= new TLine(pos,0.,pos,h_ysieve_overlay1->GetMaximum());
    ys_text_2[nys]= new TText(pos,-30,Form("%d",nys));
     ys_text_2[nys]->SetTextColor(2);
     ys_line_2[nys]->SetLineColor(6);
     ys_line_2[nys]->SetLineWidth(1);
  }
 TLine* xs_line_2[11];
 TText* xs_text_2[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*2.5-2.5*5;
    xs_line_2[nys]= new TLine(pos,0.,pos,h_xsieve_overlay1->GetMaximum());
    xs_text_2[nys]= new TText(pos,-30,Form("%d",nys));
     xs_text_2[nys]->SetTextColor(2);
     xs_line_2[nys]->SetLineColor(6);
     xs_line_2[nys]->SetLineWidth(1);
  }

  //For 1D sieve plots
  TLine* ys_line_3[11];
  TText* ys_text_3[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*1.64-1.64*5;
    ys_line_3[nys]= new TLine(pos,0.,pos,h_ysieve1_all->GetMaximum());
    ys_text_3[nys]= new TText(pos,-55,Form("%d",nys));
     ys_text_3[nys]->SetTextColor(2);
     ys_line_3[nys]->SetLineColor(6);
     ys_line_3[nys]->SetLineWidth(1);
  }
 TLine* xs_line_3[11];
 TText* xs_text_3[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*2.5-2.5*5;
    xs_line_3[nys]= new TLine(pos,0.,pos,h_xsieve1_all->GetMaximum());
    xs_text_3[nys]= new TText(pos,-55,Form("%d",nys));
     xs_text_3[nys]->SetTextColor(2);
     xs_line_3[nys]->SetLineColor(6);
     xs_line_3[nys]->SetLineWidth(1);
  }

  h_z1->SetLineColor(2);
  h_z2->SetLineColor(4);
  
  TLegend *leg1 = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg1->AddEntry(h_z1,"before opt","l");
  leg1->AddEntry(h_z2,"after opt","l");
  
  //save plots
  canvas->Update();
  canvas2->Update();
  gPad->Clear();

  canvas2->Divide(1,1);
  canvas2->cd(1);
  h_z1->SetLineWidth(2);
  h_z1->SetTitle(foil_pos);
  //h_z1->SetLineColor(2);
  h_z1->GetYaxis()->SetRangeUser(0.0,h_z2->GetMaximum()+200);
  h_z1->Draw();
  h_z2->SetLineWidth(2);
  //h_z2->SetLineColor(4);
  h_z2->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_z1_all->SetLineWidth(2);
  h_z1_all->SetTitle(all_foil);
  h_z1_all->SetLineColor(2);
  h_z1_all->GetYaxis()->SetRangeUser(0.0,h_z2_all->GetMaximum()+200);
  h_z1_all->Draw();
  h_z2_all->SetLineWidth(2);
  h_z2_all->SetLineColor(4);
  h_z2_all->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_z_res->SetLineWidth(2);
  h_z_res->SetLineColor(1);
  h_z_res->Draw("");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_Vy1->SetLineWidth(2);
  h_Vy1->SetLineColor(2);
  h_Vy1->GetYaxis()->SetRangeUser(0.0,h_Vy2->GetMaximum()+200);
  h_Vy1->Draw();
  h_Vy2->SetLineWidth(2);
  h_Vy2->SetLineColor(4);
  h_Vy2->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_Vy1_all->SetLineWidth(2);
  h_Vy1_all->SetLineColor(2);
  h_Vy1_all->GetYaxis()->SetRangeUser(0.0,h_Vy2_all->GetMaximum()+200);
  h_Vy1_all->Draw();
  h_Vy2_all->SetLineWidth(2);
  h_Vy2_all->SetLineColor(4);
  h_Vy2_all->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());
  
  canvas2->cd(1);
  h_xsieve1->SetLineWidth(2);
  h_xsieve1->SetTitle(foil_pos);
  h_xsieve1->SetLineColor(2);
  h_xsieve1->GetYaxis()->SetRangeUser(0.0,h_xsieve2->GetMaximum()+25);
  h_xsieve1->Draw();
  h_xsieve2->SetLineWidth(2);
  h_xsieve2->SetLineColor(4);
  h_xsieve2->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_xsieve_overlay1->SetLineWidth(2);
  h_xsieve_overlay1->SetTitle(foil_pos);
  h_xsieve_overlay1->SetLineColor(2);
  h_xsieve_overlay1->GetYaxis()->SetRangeUser(0.0,h_xsieve_overlay2->GetMaximum()+25);
  h_xsieve_overlay1->Draw();
  h_xsieve_overlay2->SetLineWidth(2);
  h_xsieve_overlay2->SetLineColor(4);
  h_xsieve_overlay2->Draw("same");
  leg1->Draw("same");
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line_2[nys]->Draw();
    xs_text_2[nys]->Draw();
  }
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_ysieve1->SetLineWidth(2);
  h_ysieve1->SetTitle(foil_pos);
  h_ysieve1->SetLineColor(2);
  h_ysieve1->GetYaxis()->SetRangeUser(0.0,h_ysieve2->GetMaximum()+25);
  h_ysieve1->Draw();
  h_ysieve2->SetLineWidth(2);
  h_ysieve2->SetLineColor(4);
  h_ysieve2->Draw("same");
  leg1->Draw("same");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_ysieve_overlay1->SetLineWidth(2);
  h_ysieve_overlay1->SetTitle(foil_pos);
  h_ysieve_overlay1->SetLineColor(2);
  h_ysieve_overlay1->GetYaxis()->SetRangeUser(0.0,h_ysieve_overlay2->GetMaximum()+25);
  h_ysieve_overlay1->Draw();
  h_ysieve_overlay2->SetLineWidth(2);
  h_ysieve_overlay2->SetLineColor(4);
  h_ysieve_overlay2->Draw("same");
  leg1->Draw("same");
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    ys_line_2[nys]->Draw();
    ys_text_2[nys]->Draw();
  }
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

   canvas2->cd(1);
  h_xsieve1_all->SetLineWidth(2);
  h_xsieve1_all->SetTitle(all_foil);
  h_xsieve1_all->SetLineColor(2);
  h_xsieve1_all->GetYaxis()->SetRangeUser(0.0,h_xsieve2_all->GetMaximum()+25);
  h_xsieve1_all->Draw();
  h_xsieve2_all->SetLineWidth(2);
  h_xsieve2_all->SetLineColor(4);
  h_xsieve2_all->Draw("same");
  leg1->Draw("same");
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line_3[nys]->Draw();
    xs_text_3[nys]->Draw();
  }
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

   canvas2->cd(1);
  h_ysieve1_all->SetLineWidth(2);
  h_ysieve1_all->SetTitle(all_foil);
  h_ysieve1_all->SetLineColor(2);
  h_ysieve1_all->GetYaxis()->SetRangeUser(0.0,h_ysieve2_all->GetMaximum()+50);
  h_ysieve1_all->Draw();
  h_ysieve2_all->SetLineWidth(2);
  h_ysieve2_all->SetLineColor(4);
  h_ysieve2_all->Draw("same");
  leg1->Draw("same");
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    ys_line_3[nys]->Draw();
    ys_text_3[nys]->Draw();
  }
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_xsieve_res->SetLineWidth(2);
  h_xsieve_res->SetLineColor(1);
  h_xsieve_res->Draw("");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_ysieve_res->SetLineWidth(2);
  h_ysieve_res->SetLineColor(1);
  h_ysieve_res->Draw("");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_xsieve_all_res->SetLineWidth(2);
  h_xsieve_all_res->SetLineColor(1);
  h_xsieve_all_res->Draw("");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());

  canvas2->cd(1);
  h_ysieve_all_res->SetLineWidth(2);
  h_ysieve_all_res->SetLineColor(1);
  h_ysieve_all_res->Draw("");
  gPad->SetGrid();
  canvas2->Print((pdf_file_name +"(").c_str());
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  h2_ypVz1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_ypVz2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_dpVz1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_dpVz2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_ypVy1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_ypVy2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_yfpVxfp1->Draw("colz");
  gPad->SetGrid();  
  canvas->cd(2);
  h2_yfpVxfp2->Draw("colz"); 
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_sieve1->Draw("colz"); 
  gPad->SetGrid();
  canvas->cd(2);
  h2_sieve2->Draw("colz");
  gPad->SetGrid();  
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_sieve_overlay1->Draw("colz");
  gPad->SetGrid();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line[nys]->Draw();
    xs_text[nys]->Draw();
    ys_line[nys]->Draw();
    ys_text[nys]->Draw();
  }
  canvas->cd(2);
  h2_sieve_overlay2->Draw("colz");
  gPad->SetGrid();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line[nys]->Draw();
    xs_text[nys]->Draw();
    ys_line[nys]->Draw();
    ys_text[nys]->Draw();
  }
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_sieve1_all->Draw("colz");
  gPad->SetGrid();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line[nys]->Draw();
    xs_text[nys]->Draw();
    ys_line[nys]->Draw();
    ys_text[nys]->Draw();
  }
  canvas->cd(2);
  h2_sieve2_all->Draw("colz");
  gPad->SetGrid();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line[nys]->Draw();
    xs_text[nys]->Draw();
    ys_line[nys]->Draw();
    ys_text[nys]->Draw();
  }
  canvas->Print((pdf_file_name +"(").c_str());
  
  canvas->cd(1);
  h2_xpVd1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_xpVd2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h_z_c1->SetLineWidth(2);
  h_z_c1->Draw();
  gPad->SetGrid();
  canvas->cd(2);
  h_z_c2->SetLineWidth(2);
  h_z_c2->Draw();
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_ypVz_c1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_ypVz_c2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_dpVz_c1->Draw("colz");
  canvas->cd(2);
  h2_dpVz_c2->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_ypVy_c1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_ypVy_c2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_yfpVxfp_c1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_yfpVxfp_c2->Draw("colz");  
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_sieve_c1->Draw("colz"); 
  gPad->SetGrid();
  canvas->cd(2);
  h2_sieve_c2->Draw("colz");  
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  canvas->cd(1);
  h2_xpVd_c1->Draw("colz");
  gPad->SetGrid();
  canvas->cd(2);
  h2_xpVd_c2->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

  //  for (int jj=0;jj<8; jj++){
  //h2_ypVzSlice[jj]->Draw("colz");
  // canvas->Print((pdf_file_name +"(").c_str());
  // }

  canvas->cd(1);
  h2_xpVd_c1->Draw("colz");
  gPad->SetGrid();
  canvas->Print((pdf_file_name +")").c_str());

  fout->Write();
  fout->Close();

}
