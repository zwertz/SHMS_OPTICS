void plotOpticsCompare(int nrun1, int nrun2){

  //read the input file
  // TFile *f = new TFile(Form("/net/cdaq/cdaql3data/cdaq/hallc-online/ROOTfiles/coin_replay_production_%i_-1.root",nrun));
  TFile *f1 = new TFile(Form("/volatile/hallc/spring17/holly/shmsoptics/shms_replay_production_all_%i_-1.root",nrun1));
 TFile *f2 = new TFile(Form("/volatile/hallc/spring17/holly/shmsoptics/shms_replay_production_all_%i_-1.root",nrun2));
  TTree *tt1 = (TTree*)f1->Get("T");
  TTree *tt2 = (TTree*)f2->Get("T");

  //here's the cut
  TCut cut = "P.gtr.dp<20&&P.gtr.dp>-15&&P.cal.etracknorm>0.8&&P.ngcer.npeSum>5";
 TCut cutCentral = "abs(P.gtr.x+P.gtr.th*253.0)<1&&abs((-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph))<1";

  //make the output file
  TCanvas *canvas = new TCanvas("canvas","plots",1100,800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  std::string pdf_file_name= Form("output_plots_%i_%i_sieveComparison.pdf",nrun1,nrun2);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("output_plots_%i_%i_sieveComparison.root",nrun1,nrun2),"RECREATE");

  //make plots
  TH1F *h_z1 = new TH1F("h_z1","nominal;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy1 = new TH2F("h2_ypVy1","yTar vs ypTar (nom);yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp1 = new TH2F("h2_yfpVxfp1","xfp vs yfp (nom);xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz1 = new TH2F("h2_dpVz1","delta vs z (nom);zVertex [cm];delta",100,-14,14,100,-15,20);
  TH2F *h2_ypVz1 = new TH2F("h2_ypVz1","z vs ypTar (nom);zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve1 = new TH2F("h2_sieve1","ysieve vs xsieve (nom);ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
  TH2F *h2_xpVd1 = new TH2F("h2_xpVd1","delta vs xpfp (nom);delta;xpfp",100,-10,20,100,-0.15,0.15);
  
  //plots with central hole only
  TH1F *h_z_c1 = new TH1F("h_z_c1","central sieve hole-nom;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy_c1 = new TH2F("h2_ypVy_c1","central sieve hole-nom;yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp_c1 = new TH2F("h2_yfpVxfp_c1","central sieve hole-nom;yfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz_c1 = new TH2F("h2_dpVz_c1","central sieve hole-nom;zVertex [cm];delta",100,-14,14,100,-10,20);
  TH2F *h2_ypVz_c1 = new TH2F("h2_ypVz_c1","central sieve hole-nom;zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve_c1 = new TH2F("h2_sieve_c1","central sieve hole-nom;ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
  TH2F *h2_xpVd_c1 = new TH2F("h2_xpVd_c1","central sieve hole-nom;delta;xpfp",100,-15,20,100,-0.15,0.15);

  //plots for other run
 TH1F *h_z2 = new TH1F("h_z2","run2;P.react.z [cm]",100,-14,14);
  TH2F *h2_ypVy2 = new TH2F("h2_ypVy2","yTar vs ypTar (run2);yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp2 = new TH2F("h2_yfpVxfp2","xfp vs yfp (run2);xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz2 = new TH2F("h2_dpVz2","z vs delta (run2);zVertex [cm];delta",100,-14,14,100,-15,20);
  TH2F *h2_ypVz2 = new TH2F("h2_ypVz2","z vs ypTar (run2);zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_sieve2 = new TH2F("h2_sieve2","ysieve vs xsieve (run2);ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
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
  tt1->Draw("P.gtr.ph:P.react.z>>h2_ypVz1",cut);
  tt1->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp1",cut);
  tt1->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy1",cut);
  tt1->Draw("P.react.z>>h_z1",cut);
  tt1->Draw("P.gtr.dp:P.react.z>>h2_dpVz1",cut);
  tt1->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve1",cut);
  tt1->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd1",cut);

  tt1->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c1",cut && cutCentral);
  tt1->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c1",cut && cutCentral);
  tt1->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c1",cut && cutCentral);
  tt1->Draw("P.react.z>>h_z_c1",cut && cutCentral);
  tt1->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c1",cut && cutCentral);
  tt1->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c1",cut && cutCentral);
  tt1->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c1",cut && cutCentral);

 tt2->Draw("P.gtr.ph:P.react.z>>h2_ypVz2",cut);
  tt2->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp2",cut);
  tt2->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy2",cut);
  tt2->Draw("P.react.z>>h_z2",cut);
  tt2->Draw("P.gtr.dp:P.react.z>>h2_dpVz2",cut);
  tt2->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve2",cut);
  tt2->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd2",cut);

  tt2->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c2",cut && cutCentral);
  tt2->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c2",cut && cutCentral);
  tt2->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c2",cut && cutCentral);
  tt2->Draw("P.react.z>>h_z_c2",cut && cutCentral);
  tt2->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c2",cut && cutCentral);
  tt2->Draw("P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + 40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c2",cut && cutCentral);
  tt2->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c2",cut && cutCentral);

  for (int ii=0; ii<8; ii++){
    h2_ypVzSlice[ii] = new TH2F(Form("h2_ypVzSlice_%i",ii),Form("xfp=%i cm +/- 0.5cm;zVertex [cm]; ypTar",ii+1),100,-15,15,100,-0.05,0.05);
    TCut slice = Form("abs(P.dc.x_fp - (%i+1))<0.5",ii);
    tt1->Draw(Form("P.gtr.ph:P.react.z>>h2_ypVzSlice_%i",ii),cut && slice);
  }

  //save plots
  canvas->Update();
  gPad->Clear();

  canvas->Divide(2,1);

  canvas->cd(1);
  h_z1->SetLineWidth(2);
  h_z1->Draw();
  gPad->SetGrid();
  canvas->cd(2);
  h_z2->SetLineWidth(2);
  h_z2->Draw();
  gPad->SetGrid();
  canvas->Print((pdf_file_name +"(").c_str());

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
