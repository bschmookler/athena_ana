void athena_particles(){

  //-------------------
  //Histograms
  //-------------------
  //Set Style
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  //gStyle->SetLabelOffset(0.01,"X");
  //gStyle->SetLabelOffset(0.01,"Y");
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleYOffset(0.9);

  TH1 *hgen1 = new TH1D("hgen1","Positive Pions",100,-6,6);
  hgen1->SetLineWidth(2);hgen1->SetLineColor(kRed);
  hgen1->GetXaxis()->SetTitle("#eta");hgen1->GetXaxis()->CenterTitle();

  TH1 *hrec1 = new TH1D("hrec1","Positive Pions",100,-6,6);
  hrec1->SetLineWidth(2);hrec1->SetLineColor(kBlue);
  hrec1->GetXaxis()->SetTitle("#eta");hrec1->GetXaxis()->CenterTitle();

  TH1 *hgen2 = new TH1D("hgen2","Negative Pions",100,-6,6);
  hgen2->SetLineWidth(2);hgen2->SetLineColor(kRed);
  hgen2->GetXaxis()->SetTitle("#eta");hgen2->GetXaxis()->CenterTitle();

  TH1 *hrec2 = new TH1D("hrec2","Negative Pions",100,-6,6);
  hrec2->SetLineWidth(2);hrec2->SetLineColor(kBlue);
  hrec2->GetXaxis()->SetTitle("#eta");hrec2->GetXaxis()->CenterTitle();

  TH1 *hgen3 = new TH1D("hgen3","Positive Kaons",100,-6,6);
  hgen3->SetLineWidth(2);hgen3->SetLineColor(kRed);
  hgen3->GetXaxis()->SetTitle("#eta");hgen3->GetXaxis()->CenterTitle();

  TH1 *hrec3 = new TH1D("hrec3","Positive Kaons",100,-6,6);
  hrec3->SetLineWidth(2);hrec3->SetLineColor(kBlue);
  hrec3->GetXaxis()->SetTitle("#eta");hrec3->GetXaxis()->CenterTitle();

  TH1 *hgen4 = new TH1D("hgen4","Negative Kaons",100,-6,6);
  hgen4->SetLineWidth(2);hgen4->SetLineColor(kRed);
  hgen4->GetXaxis()->SetTitle("#eta");hgen4->GetXaxis()->CenterTitle();

  TH1 *hrec4 = new TH1D("hrec4","Negative Kaons",100,-6,6);
  hrec4->SetLineWidth(2);hrec4->SetLineColor(kBlue);
  hrec4->GetXaxis()->SetTitle("#eta");hrec4->GetXaxis()->CenterTitle();

  TH1 *hgen5 = new TH1D("hgen5","Protons",100,-6,6);
  hgen5->SetLineWidth(2);hgen5->SetLineColor(kRed);
  hgen5->GetXaxis()->SetTitle("#eta");hgen5->GetXaxis()->CenterTitle();

  TH1 *hrec5 = new TH1D("hrec5","Protons",100,-6,6);
  hrec5->SetLineWidth(2);hrec5->SetLineColor(kBlue);
  hrec5->GetXaxis()->SetTitle("#eta");hrec5->GetXaxis()->CenterTitle();

  TH1 *hgen6 = new TH1D("hgen6","Negative Muons",100,-6,6);
  hgen6->SetLineWidth(2);hgen6->SetLineColor(kRed);
  hgen6->GetXaxis()->SetTitle("#eta");hgen6->GetXaxis()->CenterTitle();

  TH1 *hrec6 = new TH1D("hrec6","Negative Muons",100,-6,6);
  hrec6->SetLineWidth(2);hrec6->SetLineColor(kBlue);
  hrec6->GetXaxis()->SetTitle("#eta");hrec6->GetXaxis()->CenterTitle();

  TH1 *hgen7 = new TH1D("hgen7","All Electrons",100,-6,6);
  hgen7->SetLineWidth(2);hgen7->SetLineColor(kRed);
  hgen7->GetXaxis()->SetTitle("#eta");hgen7->GetXaxis()->CenterTitle();

  TH1 *hrec7 = new TH1D("hrec7","All Electrons",100,-6,6);
  hrec7->SetLineWidth(2);hrec7->SetLineColor(kBlue);
  hrec7->GetXaxis()->SetTitle("#eta");hrec7->GetXaxis()->CenterTitle();

  TH1 *hgen8 = new TH1D("hgen8","Positrons",100,-6,6);
  hgen8->SetLineWidth(2);hgen8->SetLineColor(kRed);
  hgen8->GetXaxis()->SetTitle("#eta");hgen8->GetXaxis()->CenterTitle();

  TH1 *hrec8 = new TH1D("hrec8","Positrons",100,-6,6);
  hrec8->SetLineWidth(2);hrec8->SetLineColor(kBlue);
  hrec8->GetXaxis()->SetTitle("#eta");hrec8->GetXaxis()->CenterTitle();

  TH1 *hgen9 = new TH1D("hgen9","Neutrons (E > 500 MeV)",100,-6,6);
  hgen9->SetLineWidth(2);hgen9->SetLineColor(kRed);
  hgen9->GetXaxis()->SetTitle("#eta");hgen9->GetXaxis()->CenterTitle();

  TH1 *hrec9 = new TH1D("hrec9","Neutrons (E > 500 MeV)",100,-6,6);
  hrec9->SetLineWidth(2);hrec9->SetLineColor(kBlue);
  hrec9->GetXaxis()->SetTitle("#eta");hrec9->GetXaxis()->CenterTitle();

  TH1 *hgen10 = new TH1D("hgen10","Neutral Kaons (E > 500 MeV)",100,-6,6);
  hgen10->SetLineWidth(2);hgen10->SetLineColor(kRed);
  hgen10->GetXaxis()->SetTitle("#eta");hgen10->GetXaxis()->CenterTitle();

  TH1 *hrec10 = new TH1D("hrec10","Neutral Kaons (E > 500 MeV)",100,-6,6);
  hrec10->SetLineWidth(2);hrec10->SetLineColor(kBlue);
  hrec10->GetXaxis()->SetTitle("#eta");hrec10->GetXaxis()->CenterTitle();

  TH1 *hgen11 = new TH1D("hgen11","Photons (E > 500 MeV)",100,-6,6);
  hgen11->SetLineWidth(2);hgen11->SetLineColor(kRed);
  hgen11->GetXaxis()->SetTitle("#eta");hgen11->GetXaxis()->CenterTitle();

  TH1 *hrec11 = new TH1D("hrec11","Photons (E > 500 MeV)",100,-6,6);
  hrec11->SetLineWidth(2);hrec11->SetLineColor(kBlue);
  hrec11->GetXaxis()->SetTitle("#eta");hrec11->GetXaxis()->CenterTitle();


  //--------------------------------//
  //   Analyse ATHENA Simulation    //
  //--------------------------------//
  TFile *f1;
  TTree *tree;

  //Load ROOT Files
  for(int iFile=1;iFile<=1;iFile++){

    //f1 = new TFile(Form("input/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d.root",iFile));
    f1 = new TFile("input/rec-clusterfix.root");
    tree = (TTree*) f1->Get("events");

    //Leaves
    TLeaf *mc_px,*mc_py,*mc_pz;
    TLeaf *mc_status,*mc_pid;
    TLeaf *mc_mass;
    TLeaf *rec_px,*rec_py,*rec_pz,*rec_E;
    TLeaf *rec_pid;

    mc_px = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.ps.x");
    mc_py = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.ps.y");
    mc_pz = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.ps.z");
    mc_status = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.genStatus");
    mc_pid = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.pdgID");
    mc_mass = (TLeaf*) tree->GetLeaf("mcparticles","mcparticles.mass");
    rec_px = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.x");
    rec_py = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.y");
    rec_pz = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.z");
    rec_E = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.energy");
    rec_pid = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.pid");

    int nevents = tree->GetEntries();
    cout<<"Beginning Analysis of file "<<iFile<<"!"<<endl;
    cout<<"Total number of events = "<<nevents<<endl;
  
    // Loop over all events
    for(int iEvent=0;iEvent<nevents;iEvent++){
      //if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
      tree->GetEntry(iEvent);
  
      //Generated Particles
      for(int imc = 0; imc < mc_px->GetLen(); imc++){
          if(mc_status->GetValue(imc)==1){
              auto theta = atan2(sqrt(pow(mc_px->GetValue(imc),2)+pow(mc_py->GetValue(imc),2)),mc_pz->GetValue(imc));
              auto eta = -log(tan(theta/2.));
              auto energy = sqrt( pow(mc_px->GetValue(imc),2) + pow(mc_py->GetValue(imc),2) + pow(mc_pz->GetValue(imc),2) +
                                  pow(mc_mass->GetValue(imc),2) );

              if(mc_pid->GetValue(imc)==211)
                hgen1->Fill(eta);
              if(mc_pid->GetValue(imc)==-211)
                hgen2->Fill(eta);
              if(mc_pid->GetValue(imc)==321)
                hgen3->Fill(eta);
              if(mc_pid->GetValue(imc)==-321)
                hgen4->Fill(eta);
              if(mc_pid->GetValue(imc)==2212)
                hgen5->Fill(eta);
              if(mc_pid->GetValue(imc)==13)
                hgen6->Fill(eta);
              if(mc_pid->GetValue(imc)==11)
                hgen7->Fill(eta);
              if(mc_pid->GetValue(imc)==-11)
                hgen8->Fill(eta);
              if(mc_pid->GetValue(imc)==2112 && energy>0.5)
                hgen9->Fill(eta);
              if(mc_pid->GetValue(imc)==130 && energy>0.5)
                hgen10->Fill(eta);
              if(mc_pid->GetValue(imc)==22 && energy>0.5)
                hgen11->Fill(eta);

          }
      }

      //Reconstructed Particles
        for(int irec = 0; irec < rec_px->GetLen(); irec++){
            auto theta = atan2(sqrt(pow(rec_px->GetValue(irec),2)+pow(rec_py->GetValue(irec),2)),rec_pz->GetValue(irec));
            auto eta = -log(tan(theta/2.));
            auto energy = rec_E->GetValue(irec);

            if(rec_pid->GetValue(irec)==211)
              hrec1->Fill(eta);
            if(rec_pid->GetValue(irec)==-211)
              hrec2->Fill(eta);
            if(rec_pid->GetValue(irec)==321)
              hrec3->Fill(eta);
            if(rec_pid->GetValue(irec)==-321)
              hrec4->Fill(eta);
            if(rec_pid->GetValue(irec)==2212)
              hrec5->Fill(eta);
            if(rec_pid->GetValue(irec)==13)
              hrec6->Fill(eta);
            if(rec_pid->GetValue(irec)==11)
              hrec7->Fill(eta);
            if(rec_pid->GetValue(irec)==-11)
              hrec8->Fill(eta);
            if(rec_pid->GetValue(irec)==2112 && energy>0.5)
              hrec9->Fill(eta);
            if(rec_pid->GetValue(irec)==130 && energy>0.5)
              hrec10->Fill(eta);
            if(rec_pid->GetValue(irec)==22 && energy>0.5)
              hrec11->Fill(eta);  
                    
        }

    }//Finished Event Loop
  }//Finished File Loop

  //Make Latex
  TPaveText* tex_energy = new TPaveText(0.15,0.8,0.425,0.9,"NDCNB");
  tex_energy->AddText("18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
  tex_energy->AddText("Q^{2} > 10 GeV^{2}");
  //tex_energy->SetTextAlign(1);
	tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(15);

  TPaveText* tex_gen = new TPaveText(0.45,0.825,0.65,0.925,"NDCNB");
  tex_gen->AddText("Generated");
  tex_gen->SetTextColor(kRed);
	tex_gen->SetFillStyle(4000);tex_gen->SetTextFont(63);tex_gen->SetTextSize(15);

  TPaveText* tex_rec = new TPaveText(0.45,0.8,0.65,0.85,"NDCNB");
  tex_rec->AddText("Reconstructed");
  tex_rec->SetTextColor(kBlue);
	tex_rec->SetFillStyle(4000);tex_rec->SetTextFont(63);tex_rec->SetTextSize(15);

  //Make Plots
  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  hgen1->Draw();hrec1->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogy();
  hgen2->Draw();hrec2->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogy();
  hgen3->Draw();hrec3->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c4 = new TCanvas("c4");
  c4->SetLogy();
  hgen4->Draw();hrec4->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogy();
  hgen5->Draw();hrec5->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c6 = new TCanvas("c6");
  c6->SetLogy();
  hgen6->Draw();hrec6->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c7 = new TCanvas("c7");
  c7->SetLogy();
  hgen7->Draw();hrec7->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c8 = new TCanvas("c8");
  c8->SetLogy();
  hgen8->Draw();hrec8->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c9 = new TCanvas("c9");
  c9->SetLogy();
  hgen9->Draw();hrec9->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c10 = new TCanvas("c10");
  c10->SetLogy();
  hgen10->Draw();hrec10->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  TCanvas *c11 = new TCanvas("c11");
  c11->SetLogy();
  hgen11->Draw();hrec11->Draw("same");
  tex_energy->Draw();tex_gen->Draw();tex_rec->Draw();

  //Print to File
  c1->Print("plots/athena_particles.pdf[");
  c1->Print("plots/athena_particles.pdf");
  c2->Print("plots/athena_particles.pdf");
  c3->Print("plots/athena_particles.pdf");
  c4->Print("plots/athena_particles.pdf");
  c5->Print("plots/athena_particles.pdf");
  c6->Print("plots/athena_particles.pdf");
  c7->Print("plots/athena_particles.pdf");
  c8->Print("plots/athena_particles.pdf");
  c9->Print("plots/athena_particles.pdf");
  c10->Print("plots/athena_particles.pdf");
  c11->Print("plots/athena_particles.pdf");
  c11->Print("plots/athena_particles.pdf]");
}