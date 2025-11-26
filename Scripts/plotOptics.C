void plotOptics(const char *nruns, double foilPOS){

  //read the input file
  //TFile *f = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1.root",nruns));
  //TFile *f = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_global_zbin_allA1n.root",nruns));
  TFile *f = new TFile(Form("/lustre24/expphy/volatile/hallc/c-lad/ewertz/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_Optics_%s_0_0_-1_newfit_lad_shms.root",nruns));
  
  TTree *tt = (TTree*)f->Get("T");

  //foil cut
  double foilmin = foilPOS-2.5;
  double foilmax = foilPOS+2.5;
  cout << foilmin << " is lower foilcut" << endl;
  cout << foilmax << " is upper foilcut" << endl;
  
  //here's the cut
  TCut cut = Form("P.gtr.dp<20&&P.gtr.dp>-15&&P.cal.etracknorm>0.8&&P.ngcer.npeSum>5&&P.react.z>%f&&P.react.z<%f",foilmin,foilmax);
  TCut cutCentral = "abs(P.gtr.x+P.gtr.th*253.0)<1&&abs((-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph))<1";

  //make the output file
  TCanvas *canvas = new TCanvas("canvas","plots",1200,800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  std::string pdf_file_name= Form("output_plots_%s_sieve0.pdf",nruns);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  canvas->SetGridx();
  canvas->SetGridy();
  TFile *fout = new TFile(Form("output_plots_%s_sieve0.root",nruns),"RECREATE");


  //make plots
  TH1F *h_z = new TH1F("h_z",";P.react.z [cm]",100,-14,14);
  TH1F *h_xsieve = new TH1F("h_xsieve",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve = new TH1F("h_ysieve",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_xsieve_overlay = new TH1F("h_xsieve_overlay",";xSieve[cm]",250,-14.0,14.0); 
  TH1F *h_ysieve_overlay = new TH1F("h_ysieve_overlay",";ySieve[cm]",250,-9.0,9.0);
  TH1F *h_Vy = new TH1F("h_Vy",";yTar [cm];",100,-5,5);
  TH2F *h2_ypVy = new TH2F("h2_ypVy",";yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp = new TH2F("h2_yfpVxfp",";xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz = new TH2F("h2_dpVz",";zVertex [cm];delta",100,-14,14,100,-15,20);
  TH2F *h2_ypVz = new TH2F("h2_ypVz",";zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve = new TH2F("h2_sieve",";ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_sieve_overlay = new TH2F("h2_sieve_overlay",";ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_xpVd = new TH2F("h2_xpVd",";delta;xpfp",100,-10,20,100,-0.15,0.15);
  
  //plots with central hole only
  TH1F *h_z_c = new TH1F("h_z_c","central sieve hole;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy_c = new TH2F("h2_ypVy_c","central sieve hole;yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp_c = new TH2F("h2_yfpVxfp_c","central sieve hole;yfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz_c = new TH2F("h2_dpVz_c","central sieve hole;zVertex [cm];delta",100,-14,14,100,-10,20);
  TH2F *h2_ypVz_c = new TH2F("h2_ypVz_c","central sieve hole;zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve_c = new TH2F("h2_sieve_c","central sieve hole;ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_sieve_c_overlay = new TH2F("h2_sieve_c_overlay","central sieve hole;ySieve [cm];xSieve[cm]",250,-9.0,9.0,250,-14.0,14.0);
  TH2F *h2_xpVd_c = new TH2F("h2_xpVd_c","central sieve hole;delta;xpfp",100,-15,20,100,-0.15,0.15);

  TH2F *h2_ypVzSlice[8];


  //plot this stuff
  tt->Draw("P.gtr.ph:P.react.z>>h2_ypVz",cut);
  tt->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp",cut);
  tt->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy",cut);
  tt->Draw("P.react.z>>h_z",cut);
  tt->Draw("P.extcor.xsieve>>h_xsieve",cut);
  tt->Draw("P.extcor.ysieve>>h_ysieve",cut);
  tt->Draw("P.extcor.xsieve>>h_xsieve_overlay",cut);
  tt->Draw("P.extcor.ysieve>>h_ysieve_overlay",cut);
  tt->Draw("P.gtr.y>>h_Vy",cut);
  //tt->Draw("P.gtr.x+P.gtr.th*253.0>>h_xsieve",cut);
  //tt->Draw("(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h_ysieve",cut);
  tt->Draw("P.gtr.dp:P.react.z>>h2_dpVz",cut);
  //tt->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve",cut);
  //tt->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_overlay",cut);
  tt->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve",cut);
  tt->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_overlay",cut);
  tt->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd",cut);

  tt->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c",cut && cutCentral);
  tt->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c",cut && cutCentral);
  tt->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c",cut && cutCentral);
  tt->Draw("P.react.z>>h_z_c",cut && cutCentral);
  tt->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c",cut && cutCentral);
  //tt->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c",cut && cutCentral);
  //tt->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c_overlay",cut && cutCentral);
  tt->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_c",cut && cutCentral);
  tt->Draw("P.extcor.xsieve:P.extcor.ysieve>>h2_sieve_c_overlay",cut && cutCentral);
  tt->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c",cut && cutCentral);

  for (int ii=0; ii<8; ii++){
    h2_ypVzSlice[ii] = new TH2F(Form("h2_ypVzSlice_%i",ii),Form("xfp=%i cm +/- 0.5cm;zVertex [cm]; ypTar",ii+1),100,-15,15,100,-0.05,0.05);
    TCut slice = Form("abs(P.dc.x_fp - (%i+1))<0.5",ii);
    tt->Draw(Form("P.gtr.ph:P.react.z>>h2_ypVzSlice_%i",ii),cut && slice);
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
    ys_line_2[nys]= new TLine(pos,0.,pos,h_ysieve_overlay->GetMaximum());
    ys_text_2[nys]= new TText(pos,-30,Form("%d",nys));
     ys_text_2[nys]->SetTextColor(2);
     ys_line_2[nys]->SetLineColor(2);
     ys_line_2[nys]->SetLineWidth(1);
  }
 TLine* xs_line_2[11];
 TText* xs_text_2[11];
  for (Int_t nys=0;nys<11;nys++) {
    Double_t pos=nys*2.5-2.5*5;
    xs_line_2[nys]= new TLine(pos,0.,pos,h_xsieve_overlay->GetMaximum());
    xs_text_2[nys]= new TText(pos,-30,Form("%d",nys));
     xs_text_2[nys]->SetTextColor(2);
     xs_line_2[nys]->SetLineColor(2);
     xs_line_2[nys]->SetLineWidth(1);
  }
  
  //save plots
  canvas->Update();

  h_z->SetLineWidth(2);
  h_z->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h_Vy->SetLineWidth(2);
  h_Vy->Draw();
  canvas->Print((pdf_file_name +"(").c_str());
  
  h_xsieve->SetLineWidth(2);
  h_xsieve->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h_ysieve->SetLineWidth(2);
  h_ysieve->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h_xsieve_overlay->SetLineWidth(2);
  h_xsieve_overlay->Draw();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    xs_line_2[nys]->Draw();
    xs_text_2[nys]->Draw();
  }
  canvas->Print((pdf_file_name +"(").c_str());

  h_ysieve_overlay->SetLineWidth(2);
  h_ysieve_overlay->Draw();
   //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    ys_line_2[nys]->Draw();
    ys_text_2[nys]->Draw();
  }
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVz->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_dpVz->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVy->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_yfpVxfp->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());

  h2_sieve->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());
  
  h2_sieve_overlay->Draw("colz");
  //Actually draw the lines
  for (Int_t nys=0;nys<11;nys++) {
    ys_line[nys]->Draw();
    ys_text[nys]->Draw();
    xs_line[nys]->Draw();
    xs_text[nys]->Draw();
  }
  canvas->Print((pdf_file_name +"(").c_str());

  h2_xpVd->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h_z_c->SetLineWidth(2);
  h_z_c->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVz_c->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_dpVz_c->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVy_c->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_yfpVxfp_c->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());

  h2_sieve_c->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());

  h2_sieve_c->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());
  
  h2_xpVd_c->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  for (int jj=0;jj<8; jj++){
    h2_ypVzSlice[jj]->Draw("colz");
    canvas->Print((pdf_file_name +"(").c_str());
  }
  h2_xpVd_c->Draw("colz");
  canvas->Print((pdf_file_name +")").c_str());

  fout->Write();
  fout->Close();

}
