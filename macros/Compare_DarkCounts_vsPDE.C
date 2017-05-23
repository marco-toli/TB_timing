{
  
    
    const int NPDE = 8;
    float scint_PDE[NPDE];
    scint_PDE[0] = 0.05;    
    scint_PDE[1] = 0.10;
    scint_PDE[2] = 0.15;
    scint_PDE[3] = 0.20;
    scint_PDE[4] = 0.25;
    scint_PDE[5] = 0.30;
//     scint_PDE[5] = 0.40;
    scint_PDE[6] = 0.50;
//     scint_PDE[8] = 0.60;
    scint_PDE[7] = 0.70;
    
    
    const int nNoise = 8;
    float DCounts[nNoise];		//dark counts in 1 ns per mm^2
    DCounts[0] = 0;
    DCounts[1] = 1.;
    DCounts[2] = 2.5;
    DCounts[3] = 5;
    DCounts[4] = 10;
    DCounts[5] = 25;
    DCounts[6] = 50;
    DCounts[7] = 100;
  
//     DCounts[3] = 80;
//     DCounts[4] = 160;
    
    
    const int NPH = 15;
    int nPhotons[NPH];
    nPhotons[0] = 1;
    nPhotons[1] = 2;
    nPhotons[2] = 3;
    nPhotons[3] = 4;
    nPhotons[4] = 5;
    nPhotons[5] = 7;
    nPhotons[6] = 10;
    nPhotons[7] = 15;
    nPhotons[8] = 20;
    nPhotons[9] = 30;
    nPhotons[10] = 40;
    nPhotons[11] = 50;
    nPhotons[12] = 75;
    nPhotons[13] = 100;
    nPhotons[14] = 200;

  
    TGraphErrors * g_CTR_NPH [NPDE][nNoise];
    TGraphErrors * g_CTR_NPH_corr [NPDE][nNoise];            
    
    for (int iPDE = 0; iPDE < NPDE; iPDE++)
    {
    
    
//       TFile * inputFile = new TFile(Form("./graphs_histos/12um_25PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", scint_PDE[iPDE]), "READ");
//       TFile * inputFile = new TFile(Form("./graphs_histos/12um_25PDE_10ns_decay_100RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", scint_PDE[iPDE]), "READ");
//       TFile * inputFile = new TFile(Form("./graphs_histos/GLUE_12um_25PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%.0f.root", scint_PDE[iPDE]), "READ");
//       TFile * inputFile = new TFile(Form("./graphs_histos/temp_%.0f.root", scint_PDE[iPDE]), "READ");
      TFile * inputFile = new TFile(Form("./graphs_histos/PDE_scan/ped_subs/out_cms_tile_3mm_12x12_wide_sipm5x5_PDE_%.0f.root", scint_PDE[iPDE]*100), "READ");

      for (int iNoise = 0; iNoise <nNoise; iNoise++)
      {
          std::cout << "getting graph for noise: " << DCounts[iNoise] << std::endl;
	g_CTR_NPH[iPDE][iNoise] 	= (TGraphErrors*) inputFile->Get(Form("g_CTR_NPH_noise_%d", iNoise));	
	g_CTR_NPH_corr[iPDE][iNoise] = (TGraphErrors*) inputFile->Get(Form("g_CTR_NPH_corr_noise_%d", iNoise));	
	
	/*
	for (int iNPH = 0; iNPH < NPH; iNPH++)
	{
	  if (iPDE == 0 
              &&  iNPH >= (NPH-2)
// 	      || (iPDE==0 && iNoise == (nNoise-2) && iNPH >= (NPH-1))
// 	      || (iPDE==1 && iNoise == (nNoise-1) && iNPH >= (NPH-1))
	     )
	  {
	    g_CTR_NPH[iPDE][iNoise]->RemovePoint(iNPH);
	    g_CTR_NPH_corr[iPDE][iNoise]->RemovePoint(iNPH);
	  }
	}*/
      
      }
    }
    
    //read histos and create new graphs
    
    TGraphErrors * gCTR_vs_PDE[nNoise];
    TGraphErrors * gCTR_vs_Noise[NPDE];
    TGraphErrors * gCTR_vs_PDE_corr[nNoise];
    TGraphErrors * gCTR_vs_Noise_corr[NPDE];
    
    for (int iPDE = 0; iPDE < NPDE; iPDE++)      
    {
      gCTR_vs_Noise[iPDE] = new TGraphErrors ();
      gCTR_vs_Noise_corr[iPDE] = new TGraphErrors ();
    }
    for (int iNoise = 0;  iNoise  < nNoise;   iNoise++)	      
    {
      gCTR_vs_PDE[iNoise] = new TGraphErrors ();
      gCTR_vs_PDE_corr[iNoise] = new TGraphErrors ();
    }
    
    double ctr_length_noise[NPDE][nNoise][NPH];
    double min_ctr_length_noise[NPDE][nNoise];
    
    double corr_ctr_length_noise[NPDE][nNoise][NPH];
    double corr_min_ctr_length_noise[NPDE][nNoise];
    
    double temp_nph [NPH];
    
//  double temp
    
    for (int iPDE = 0; iPDE < NPDE; iPDE++)
    {
      for (int iNoise = 0; iNoise < nNoise; iNoise++)
      {
	float min_ctr 	   = 9999;
	float min_ctr_corr = 9999;
	
	for (int iNPH = 0; iNPH < NPH; iNPH++)
	{
            
	    g_CTR_NPH[iPDE][iNoise]	 ->GetPoint(iNPH, temp_nph[iNPH], ctr_length_noise[iPDE][iNoise][iNPH]);
	    g_CTR_NPH_corr[iPDE][iNoise]->GetPoint(iNPH, temp_nph[iNPH], corr_ctr_length_noise[iPDE][iNoise][iNPH]);
            
            std::cout << "iNPH = " << iNPH << " :: ctr = " <<  ctr_length_noise[iPDE][iNoise][iNPH] << " :: corr_ctr = " << corr_ctr_length_noise[iPDE][iNoise][iNPH] << std::endl;
	  
	    if (ctr_length_noise[iPDE][iNoise][iNPH] < min_ctr) min_ctr = ctr_length_noise[iPDE][iNoise][iNPH];
	    if (corr_ctr_length_noise[iPDE][iNoise][iNPH] < min_ctr_corr) min_ctr_corr = corr_ctr_length_noise[iPDE][iNoise][iNPH];	  
	}
	
	min_ctr_length_noise[iPDE][iNoise] = min_ctr;
	corr_min_ctr_length_noise[iPDE][iNoise] = min_ctr_corr;
	
// 	if (iNoise 
	gCTR_vs_PDE[iNoise]->SetPoint(iPDE, scint_PDE[iPDE]*100, min_ctr_length_noise[iPDE][iNoise]);
	gCTR_vs_Noise[iPDE]->SetPoint(iNoise, DCounts[iNoise]*1e9, min_ctr_length_noise[iPDE][iNoise]);
	
	std::cout << "CTR[" << iPDE << "][" << iNoise << "] = " << min_ctr_length_noise[iPDE][iNoise] << std::endl;
	
	gCTR_vs_PDE_corr[iNoise]->SetPoint(iPDE, scint_PDE[iPDE]*100, corr_min_ctr_length_noise[iPDE][iNoise]);
	gCTR_vs_Noise_corr[iPDE]->SetPoint(iNoise, DCounts[iNoise]*1e9, corr_min_ctr_length_noise[iPDE][iNoise]);
      }
    }
    
    
    
    ///drawing
   int sel_PDE = 4;
   std::cout << "defining latex captions..." << std::endl;
   TLegend * leg;
   
   TLatex t(.1,.91,"SIMULATION: Crystal 12x12 + SiPM 5x5"); 
   t.SetTextFont(43);
   t.SetTextSize(19);
   t.SetNDC(kTRUE);
    
   TLatex t_prel(0.7,.91,"PRELIMINARY"); 
   t_prel.SetTextFont(43);
   t_prel.SetTextSize(19);
   t_prel.SetNDC(kTRUE);
    
   
   TPaveText *pt3 = new TPaveText(0.15, 0.84, 0.4, 0.88, "brNDC");
   pt3->SetShadowColor(0);
   pt3->SetFillColor(0);
   pt3->SetLineColor(0);   
//    pt3->AddText("Cherenkov + Scintillation");
   pt3->AddText(Form("Crystal length: %.0f mm", scint_PDE[sel_PDE]));
   pt3->SetTextFont(42);
   pt3->SetTextAlign(12);
   pt3->SetTextSize(0.03);
   
   TPaveText *pt2 = new TPaveText(0.15, 0.78, 0.4, 0.83, "brNDC");
   pt2->SetShadowColor(0);
   pt2->SetFillColor(0);
   pt2->SetLineColor(0);
//    pt2->AddText(Form("PDE = %.2f ", scint_PDE[]));
   pt2->SetTextFont(42);
   pt2->SetTextAlign(12);
   pt2->SetTextSize(0.03);
      
   int SPTR = 66;
   TPaveText *pt = new TPaveText(0.15, 0.72, 0.5, 0.77, "brNDC");
   pt->SetShadowColor(0);
   pt->SetFillColor(0);
   pt->SetLineColor(0);
   pt->AddText(Form("#sigma_{SPTR} = %d ps", SPTR));
   pt->SetTextFont(42);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.03);
  
   std::cout << "drawing canvas..." << std::endl;
    
    TCanvas * cTrendNoise = new TCanvas ("cTrendNoise", "cTrendNoise", 600, 600);
    gCTR_vs_Noise_corr[0]->Draw("ALPE");
    gCTR_vs_Noise_corr[0]->GetXaxis()->SetTitle("DCR [counts per second]");
    gCTR_vs_Noise_corr[0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
    gCTR_vs_Noise_corr[0]->GetYaxis()->SetRangeUser(0, 140);
    gCTR_vs_Noise_corr[0]->GetXaxis()->SetRangeUser(0, DCounts[nNoise-1]*1.2e9);
    gCTR_vs_Noise_corr[0]->GetYaxis()->SetTitleOffset(1.3);
//     gCTR_vs_Noise[0]->SetMarkerStyle(21);
    for (int iPDE = 0; iPDE < NPDE; iPDE++)
    {
      gCTR_vs_Noise_corr[iPDE]->SetMarkerColor(iPDE+1);
      gCTR_vs_Noise_corr[iPDE]->SetLineColor(iPDE+1);
      gCTR_vs_Noise_corr[iPDE]->SetLineWidth(2);
      
      gCTR_vs_Noise_corr[iPDE]->Draw("same LPE");
    }
    if(NPDE > 4)
    {
        gCTR_vs_Noise_corr[2]->SetMarkerColor(kGreen+1);
        gCTR_vs_Noise_corr[2]->SetLineColor(kGreen+1);
        gCTR_vs_Noise_corr[4]->SetMarkerColor(kYellow+1);
        gCTR_vs_Noise_corr[4]->SetLineColor(kYellow+1);
    }
    gPad->SetGridy();
    
    leg = new TLegend(0.12,0.68,0.45,0.88,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);

    for (int iPDE = 0; iPDE < NPDE; iPDE++)
    {
       leg->AddEntry(gCTR_vs_Noise_corr[iPDE], Form("PDE = %.0f ", scint_PDE[iPDE]*100), "lp");
    }
    leg->Draw();
    t.Draw();
    t_prel.Draw();
    
    
    TCanvas * cTrendPDE = new TCanvas ("cTrendPDE", "cTrendPDE", 600, 600);
    gCTR_vs_PDE_corr[0]->Draw("ALPE");
    gCTR_vs_PDE_corr[0]->GetXaxis()->SetTitle("SiPM PDE [%]");
    gCTR_vs_PDE_corr[0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
    gCTR_vs_PDE_corr[0]->GetYaxis()->SetRangeUser(0, 140);
    
    gCTR_vs_PDE_corr[0]->GetXaxis()->SetRangeUser(0, 100);
    gCTR_vs_PDE_corr[0]->GetYaxis()->SetTitleOffset(1.3);
//     gCTR_vs_PDE[0]->SetMarkerStyle(21);
    for (int iNoise = 0; iNoise < nNoise; iNoise++)
    {
      gCTR_vs_PDE_corr[iNoise]->SetMarkerColor(iNoise+1);
      gCTR_vs_PDE_corr[iNoise]->SetLineColor(iNoise+1);
      gCTR_vs_PDE_corr[iNoise]->SetLineWidth(2);
      
      
      gCTR_vs_PDE_corr[iNoise]->Draw("same LPE");
    }
    if(nNoise > 4)
    {
        gCTR_vs_PDE_corr[2]->SetMarkerColor(kGreen+1);
        gCTR_vs_PDE_corr[2]->SetLineColor(kGreen+1);
        gCTR_vs_PDE_corr[4]->SetMarkerColor(kYellow+1);
        gCTR_vs_PDE_corr[4]->SetLineColor(kYellow+1);
    }
    gPad->SetGridy();
    gPad->SetLogy();
    gCTR_vs_PDE_corr[0]->GetYaxis()->SetRangeUser(10, 200);
    
    leg = new TLegend(0.5,0.6,0.8,0.8,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);

    for (int iNoise = 0; iNoise < nNoise; iNoise++)
    {
       leg->AddEntry(gCTR_vs_PDE_corr[iNoise], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
    }
    leg->Draw();
    t.Draw();
    t_prel.Draw();
    
    
    
    
    
    
   TCanvas * cgCTR = new TCanvas ("cgCTR", "cgCTR", 600, 600);
   g_CTR_NPH[sel_PDE][0]->Draw("ALPE");
   g_CTR_NPH[sel_PDE][0]->GetXaxis()->SetTitle("Threshold [number of photons]");
   g_CTR_NPH[sel_PDE][0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
   g_CTR_NPH[sel_PDE][0]->GetYaxis()->SetTitleOffset(1.4);
   g_CTR_NPH[sel_PDE][0]->GetYaxis()->SetRangeUser(0, 100);
   g_CTR_NPH[sel_PDE][0]->GetXaxis()->SetRangeUser(0, 100);
   g_CTR_NPH[sel_PDE][0]->SetLineWidth(2);
   
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     g_CTR_NPH[sel_PDE][iNoise]->SetLineWidth(2);
     g_CTR_NPH[sel_PDE][iNoise]->SetLineColor(iNoise+1);
     g_CTR_NPH[sel_PDE][iNoise]->SetMarkerColor(iNoise+1);
     g_CTR_NPH[sel_PDE][iNoise]->Draw("LPE same");
   }
   if(nNoise > 4)
   {
        g_CTR_NPH[sel_PDE][2]->SetMarkerColor(kGreen+1);
        g_CTR_NPH[sel_PDE][2]->SetLineColor(kGreen+1);
        g_CTR_NPH[sel_PDE][4]->SetMarkerColor(kYellow+1);
        g_CTR_NPH[sel_PDE][4]->SetLineColor(kYellow+1);
   }
   
   gPad->SetGrid();
   
   leg = new TLegend(0.58,0.68,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(g_CTR_NPH[sel_PDE][iNoise], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   pt->Draw();
//    pt2->Draw();
   pt3->Draw();
   t.Draw();
   t_prel.Draw();
   
   TCanvas * cgCTR_corr = new TCanvas ("cgCTR_corr", "cgCTR_corr", 600, 600);
   g_CTR_NPH_corr[sel_PDE][0]->Draw("ALPE");
   g_CTR_NPH_corr[sel_PDE][0]->GetXaxis()->SetTitle("Threshold [number of photons]");
   g_CTR_NPH_corr[sel_PDE][0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
   g_CTR_NPH_corr[sel_PDE][0]->GetYaxis()->SetTitleOffset(1.4);
   g_CTR_NPH_corr[sel_PDE][0]->GetXaxis()->SetLimits(0, 100);
   g_CTR_NPH_corr[sel_PDE][0]->GetYaxis()->SetRangeUser(0, 100);
   g_CTR_NPH_corr[sel_PDE][0]->SetLineWidth(2);
   
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     g_CTR_NPH_corr[sel_PDE][iNoise]->SetLineWidth(2);
     g_CTR_NPH_corr[sel_PDE][iNoise]->SetLineColor(iNoise+1);
     g_CTR_NPH_corr[sel_PDE][iNoise]->SetMarkerColor(iNoise+1);
     g_CTR_NPH_corr[sel_PDE][iNoise]->Draw("LPE same");
   }
   if(nNoise > 4)
   {
    g_CTR_NPH_corr[sel_PDE][2]->SetMarkerColor(kGreen+1);
    g_CTR_NPH_corr[sel_PDE][2]->SetLineColor(kGreen+1);
    g_CTR_NPH_corr[sel_PDE][4]->SetMarkerColor(kYellow+1);
    g_CTR_NPH_corr[sel_PDE][4]->SetLineColor(kYellow+1);
   }
   
   gPad->SetGrid();
   
   leg = new TLegend(0.58,0.68,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(g_CTR_NPH_corr[sel_PDE][iNoise], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   pt->Draw();
//    pt2->Draw();
   pt3->Draw();
   t.Draw();
   t_prel.Draw();
   
   
    string outdaq;
    outdaq = "./temp_plots/tr_vs_dc.pdf";
    const char * daqfile = outdaq.c_str();
    cTrendNoise->cd();
    cTrendNoise->SaveAs(daqfile);
    
    outdaq = "./temp_plots/tr_vs_pde.pdf";
    daqfile = outdaq.c_str();
    cTrendPDE->cd();
    cTrendPDE->SaveAs(daqfile);
    
    
    
    
    
    
}
