#include "PadMxN.h"

//=====================================//
TLorentzVector apply_boost(TLorentzVector ei,TLorentzVector pi, TLorentzVector part){

  //Step 1: Find the needed boosts and rotations from the incoming lepton and hadron beams 
  //(note, this will give you a perfect boost, in principle you will not know the beam momenta exactly and should use an average)
  
  // Define the Boost to make beams back-to-back
  TLorentzVector cmBoost = (1./ei.E())*ei + (1./pi.E())*pi;

  TLorentzVector boost(-cmBoost.Px(),-cmBoost.Py(),-cmBoost.Pz(),cmBoost.E());
  TVector3 b;
  b = boost.BoostVector();

  TLorentzVector boostBack(0.0,0.0,cmBoost.Pz(),cmBoost.E());
  TVector3 bb;
  bb = boostBack.BoostVector(); // This will boost beams from a center of momentum frame back to (nearly) their original energies

  // Boost and rotate the incoming beams to find the proper rotations TLorentzVector
  pi.Boost(b); // Boost to COM frame
  ei.Boost(b);
  double rotAboutY = -1.0*TMath::ATan2(pi.Px(),pi.Pz()); // Rotate to remove x component of beams
  double rotAboutX = 1.0*TMath::ATan2(pi.Py(),pi.Pz()); // Rotate to remove y component of beams

  //Step 2: Apply boosts and rotations to any particle 4-vector 
  //(here too, choices will have to be made as to what the 4-vector is for reconstructed particles)
  
  // Boost and rotate particle 4-momenta into the headon frame
  part.Boost(b);
  part.RotateY(rotAboutY);
  part.RotateX(rotAboutX);
  part.Boost(bb);

  return part;

}

//=====================================//
void athena_Rec_y(){

  //True y bins
  const int nybins = 5;
  float y_low[nybins] = {0.50,0.20,0.10,0.05,0.01};
  float y_hi[nybins]  = {0.80,0.50,0.20,0.10,0.05};

  //We will apply cuts on the true y
  double y_max = 0.95;double y_min = 1e-3;

  //Resolution limit
  double res_limit = 75;

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

  //y resolution (constant beam particle momenta)
  TH1 *h1a[nybins]; //Electron method
  TH1 *h1b[nybins]; //JB method
  TH1 *h1c[nybins]; //DA method
  TH1 *h1d[nybins]; //Sigma method
  TH1 *h1e[nybins]; //eSigma method

  for(int ibin=0;ibin<nybins;ibin++){
    h1a[ibin] = new TH1D(Form("h1a[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1a[ibin]->SetLineColor(kBlue);h1a[ibin]->SetLineWidth(2);
    h1a[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1a[ibin]->GetXaxis()->CenterTitle();

    h1b[ibin] = new TH1D(Form("h1b[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1b[ibin]->SetLineColor(kBlue);h1b[ibin]->SetLineWidth(2);
    h1b[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1b[ibin]->GetXaxis()->CenterTitle();

    h1c[ibin] = new TH1D(Form("h1c[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1c[ibin]->SetLineColor(kBlue);h1c[ibin]->SetLineWidth(2);
    h1c[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1c[ibin]->GetXaxis()->CenterTitle();

    h1d[ibin] = new TH1D(Form("h1d[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1d[ibin]->SetLineColor(kBlue);h1d[ibin]->SetLineWidth(2);
    h1d[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1d[ibin]->GetXaxis()->CenterTitle();

    h1e[ibin] = new TH1D(Form("h1e[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1e[ibin]->SetLineColor(kBlue);h1e[ibin]->SetLineWidth(2);
    h1e[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1e[ibin]->GetXaxis()->CenterTitle();

  }

  //Profile Histograms
  //We want to resolution, which we get by doing GetBinError() with the "s" option.
  //Otherwise it is the error on the mean.
  //---
  //1<Q2<10
  TProfile *prof1 = new TProfile("prof1","",30,0,1,-res_limit,res_limit,"s"); //Electron % deviation vs. y_true
  TProfile *prof2 = new TProfile("prof2","",30,0,1,-res_limit,res_limit,"s"); //JB % deviation vs. y_true
  TProfile *prof3 = new TProfile("prof3","",30,0,1,-res_limit,res_limit,"s"); //DA % deviation vs. y_true
  TProfile *prof4 = new TProfile("prof4","",30,0,1,-res_limit,res_limit,"s"); //Sigma % deviation vs. y_true
  TProfile *prof5 = new TProfile("prof5","",30,0,1,-res_limit,res_limit,"s"); //eSigma % deviation vs. y_true
  //---
  //10<Q2>100
  TProfile *prof1a = new TProfile("prof1a","",30,0,1,-res_limit,res_limit,"s"); //Electron % deviation vs. y_true
  TProfile *prof2a = new TProfile("prof2a","",30,0,1,-res_limit,res_limit,"s"); //JB % deviation vs. y_true
  TProfile *prof3a = new TProfile("prof3a","",30,0,1,-res_limit,res_limit,"s"); //DA % deviation vs. y_true
  TProfile *prof4a = new TProfile("prof4a","",30,0,1,-res_limit,res_limit,"s"); //Sigma % deviation vs. y_true
  TProfile *prof5a = new TProfile("prof5a","",30,0,1,-res_limit,res_limit,"s"); //eSigma % deviation vs. y_true
  //---
  //100<Q2<1000
  TProfile *prof1b = new TProfile("prof1b","",30,0,1,-res_limit,res_limit,"s"); //Electron % deviation vs. y_true
  TProfile *prof2b = new TProfile("prof2b","",30,0,1,-res_limit,res_limit,"s"); //JB % deviation vs. y_true
  TProfile *prof3b = new TProfile("prof3b","",30,0,1,-res_limit,res_limit,"s"); //DA % deviation vs. y_true
  TProfile *prof4b = new TProfile("prof4b","",30,0,1,-res_limit,res_limit,"s"); //Sigma % deviation vs. y_true
  TProfile *prof5b = new TProfile("prof5b","",30,0,1,-res_limit,res_limit,"s"); //eSigma % deviation vs. y_true

  //Try also with histograms
  double fine_1d_left[30];
  double fine_1d_width = 1./30.;
  TH1 *hbin1[30],*hbin2[30],*hbin3[30],*hbin4[30],*hbin5[30];
  TH1 *hbin1a[30],*hbin2a[30],*hbin3a[30],*hbin4a[30],*hbin5a[30];
  TH1 *hbin1b[30],*hbin2b[30],*hbin3b[30],*hbin4b[30],*hbin5b[30];

  for(int ibin=0;ibin<30;ibin++){

      fine_1d_left[ibin] = ibin/30.;

      hbin1[ibin] = new TH1D(Form("hbin1[%d]",ibin),"",500,-150,150);
      hbin2[ibin] = new TH1D(Form("hbin2[%d]",ibin),"",500,-150,150);
      hbin3[ibin] = new TH1D(Form("hbin3[%d]",ibin),"",500,-150,150);
      hbin4[ibin] = new TH1D(Form("hbin4[%d]",ibin),"",500,-150,150);
      hbin5[ibin] = new TH1D(Form("hbin5[%d]",ibin),"",500,-150,150);

      hbin1a[ibin] = new TH1D(Form("hbin1a[%d]",ibin),"",500,-150,150);
      hbin2a[ibin] = new TH1D(Form("hbin2a[%d]",ibin),"",500,-150,150);
      hbin3a[ibin] = new TH1D(Form("hbin3a[%d]",ibin),"",500,-150,150);
      hbin4a[ibin] = new TH1D(Form("hbin4a[%d]",ibin),"",500,-150,150);
      hbin5a[ibin] = new TH1D(Form("hbin5a[%d]",ibin),"",500,-150,150);
      
      hbin1b[ibin] = new TH1D(Form("hbin1b[%d]",ibin),"",500,-150,150);
      hbin2b[ibin] = new TH1D(Form("hbin2b[%d]",ibin),"",500,-150,150);
      hbin3b[ibin] = new TH1D(Form("hbin3b[%d]",ibin),"",500,-150,150);
      hbin4b[ibin] = new TH1D(Form("hbin4b[%d]",ibin),"",500,-150,150);
      hbin5b[ibin] = new TH1D(Form("hbin5b[%d]",ibin),"",500,-150,150);

  }

  //Two Dimensional profiles and histograms
  //-------
  //Q2 Binning                                                                                                                                        
  double Q2_min = 1E-1;  
  double Q2_max = 1E4;
  const int nbins_Q2 = 25;
  double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/(nbins_Q2);   //Determine bin width                                                              
  double log_Q2_div;
  double Q2_bins[nbins_Q2+1];
  for(int i=0;i<nbins_Q2+1;i++){
    log_Q2_div = log10(Q2_min) + (i*log_bw_Q2);
    Q2_bins[i] = pow(10,log_Q2_div);
  }

  //x Binning                                                                                                                                        
  double x_min = 1E-5;
  double x_max = 1;
  const int nbins_x = 25;
  double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x);   //Determine bin width                                                                  
  double log_x_div;
  double x_bins[nbins_x+1];
  for(int i=0;i<nbins_x+1;i++){
    log_x_div = log10(x_min) + (i*log_bw_x);
    x_bins[i] = pow(10,log_x_div);
  }

  TProfile2D *profh1 = new TProfile2D("profh1","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //Electron Method
  TProfile2D *profh2 = new TProfile2D("profh2","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //JB Method
  TProfile2D *profh3 = new TProfile2D("profh3","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //DA Method
  TProfile2D *profh4 = new TProfile2D("profh4","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //Sigma Method
  TProfile2D *profh5 = new TProfile2D("profh5","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //eSigma Method

  TH2 *hh1 = new TH2D("hh1","",nbins_x,x_bins,nbins_Q2,Q2_bins); //Electron Method
  hh1->GetXaxis()->SetTitle("x");hh1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh1->GetXaxis()->SetLabelFont(63);hh1->GetYaxis()->SetLabelFont(63);
  hh1->GetXaxis()->SetLabelSize(25);hh1->GetYaxis()->SetLabelSize(25);
  hh1->GetXaxis()->SetLabelOffset(0.01);hh1->GetYaxis()->SetLabelOffset(0.01);
  hh1->GetXaxis()->CenterTitle(1);hh1->GetYaxis()->CenterTitle(1);
  hh1->GetXaxis()->SetTitleSize(40);hh1->GetXaxis()->SetTitleOffset(2.5); 
  hh1->GetYaxis()->SetTitleSize(40);hh1->GetYaxis()->SetTitleOffset(3.0);

  TH2 *hh2 = new TH2D("hh2","",nbins_x,x_bins,nbins_Q2,Q2_bins); //JB Method
  hh2->GetXaxis()->SetTitle("x");hh2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh2->GetXaxis()->SetLabelFont(63);hh2->GetYaxis()->SetLabelFont(63);
  hh2->GetXaxis()->SetLabelSize(25);hh2->GetYaxis()->SetLabelSize(25);
  hh2->GetXaxis()->SetLabelOffset(0.01);hh2->GetYaxis()->SetLabelOffset(0.01);
  hh2->GetXaxis()->CenterTitle(1);hh2->GetYaxis()->CenterTitle(1);
  hh2->GetXaxis()->SetTitleSize(40);hh2->GetXaxis()->SetTitleOffset(2.5); 
  hh2->GetYaxis()->SetTitleSize(40);hh2->GetYaxis()->SetTitleOffset(3.0);
  
  TH2 *hh3 = new TH2D("hh3","",nbins_x,x_bins,nbins_Q2,Q2_bins); //DA Method
  hh3->GetXaxis()->SetTitle("x");hh3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh3->GetXaxis()->SetLabelFont(63);hh3->GetYaxis()->SetLabelFont(63);
  hh3->GetXaxis()->SetLabelSize(25);hh3->GetYaxis()->SetLabelSize(25);
  hh3->GetXaxis()->SetLabelOffset(0.01);hh3->GetYaxis()->SetLabelOffset(0.01);
  hh3->GetXaxis()->CenterTitle(1);hh3->GetYaxis()->CenterTitle(1);
  hh3->GetXaxis()->SetTitleSize(40);hh3->GetXaxis()->SetTitleOffset(2.5); 
  hh3->GetYaxis()->SetTitleSize(40);hh3->GetYaxis()->SetTitleOffset(3.0);

  TH2 *hh4 = new TH2D("hh4","",nbins_x,x_bins,nbins_Q2,Q2_bins); //Sigma Method
  hh4->GetXaxis()->SetTitle("x");hh4->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh4->GetXaxis()->SetLabelFont(63);hh4->GetYaxis()->SetLabelFont(63);
  hh4->GetXaxis()->SetLabelSize(25);hh4->GetYaxis()->SetLabelSize(25);
  hh4->GetXaxis()->SetLabelOffset(0.01);hh4->GetYaxis()->SetLabelOffset(0.01);
  hh4->GetXaxis()->CenterTitle(1);hh4->GetYaxis()->CenterTitle(1);
  hh4->GetXaxis()->SetTitleSize(40);hh4->GetXaxis()->SetTitleOffset(2.5); 
  hh4->GetYaxis()->SetTitleSize(40);hh4->GetYaxis()->SetTitleOffset(3.0);

  TH2 *hh5 = new TH2D("hh5","",nbins_x,x_bins,nbins_Q2,Q2_bins); //eSigma Method
  hh5->GetXaxis()->SetTitle("x");hh5->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh5->GetXaxis()->SetLabelFont(63);hh5->GetYaxis()->SetLabelFont(63);
  hh5->GetXaxis()->SetLabelSize(25);hh5->GetYaxis()->SetLabelSize(25);
  hh5->GetXaxis()->SetLabelOffset(0.01);hh5->GetYaxis()->SetLabelOffset(0.01);
  hh5->GetXaxis()->CenterTitle(1);hh5->GetYaxis()->CenterTitle(1);
  hh5->GetXaxis()->SetTitleSize(40);hh5->GetXaxis()->SetTitleOffset(2.5); 
  hh5->GetYaxis()->SetTitleSize(40);hh5->GetYaxis()->SetTitleOffset(3.0);

  //-------

  //Set constants
  const double Mp = 0.9383;
  const double Me = 0.511e-3;
  const double crossingAngle = -0.025;
  const TLorentzVector ei_const(0,0,-18.,sqrt(18.*18. + Me*Me) );
  const TLorentzVector pi_const(275.*sin(crossingAngle),0,275.*cos(crossingAngle),sqrt(275.*275. + Mp*Mp) );

  //--------------------------------//
  //   Analyse ATHENA Simulation    //
  //--------------------------------//

  //Too many files...subdivide
  for(int ichain=1;ichain<=19;ichain++){
    
    TChain *chain = new TChain("events");
    chain->Reset();

    if(ichain<=5)
      chain->Add(Form("input/v2.1/Q2_1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",ichain));
    else if(ichain<=10)
      chain->Add(Form("input/v2.1/Q2_10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",ichain-5));
    else if(ichain<=15)
      chain->Add(Form("input/v2.1/Q2_100/pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",ichain-10));
    else if(ichain<=19)
      chain->Add(Form("input/v2.1/Q2_1000/pythia8NCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",ichain-15));

    cout<<"-------"<<endl;
    cout<<"Running Chain Number "<<ichain<<endl;
    cout<<"Chain = "<<chain<<endl;
    cout<<"Total Number of events = "<<chain->GetEntries()<<endl;
    cout<<"-------"<<endl;
  
    //Load ROOT Files
    //TFile *f1;
    //TTree *tree;
    //int file_set[] = {1,10,100,1000};
    //for(int iFile=0;iFile<4;iFile++){
    //f1 = new TFile(Form("input/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d.root",iFile));
    //f1 = new TFile(Form("input/pythia8NCDIS_18x275_minQ2_%d_hiDiv_1.root",file_set[iFile]));
  
    TObjArray *fileElements=chain->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while ((chEl=(TChainElement*)next())) {
      TFile* f = new TFile(chEl->GetTitle());
      TTree* tree = (TTree*) f->Get("events");

      //Leaves
      TLeaf *Q2_true,*y_true,*x_true;
      TLeaf *rec_px,*rec_py,*rec_pz,*rec_E;
      TLeaf *rec_charge;
      TLeaf *rec_pid;
      TLeaf *rec_se;

      Q2_true = tree->GetLeaf("InclusiveKinematicsTruth","InclusiveKinematicsTruth.Q2");
      y_true = tree->GetLeaf("InclusiveKinematicsTruth","InclusiveKinematicsTruth.y");
      x_true = tree->GetLeaf("InclusiveKinematicsTruth","InclusiveKinematicsTruth.x");
      rec_px = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.x");
      rec_py = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.y");
      rec_pz = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.p.z");
      rec_E = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.energy");
      rec_charge = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.charge");
      rec_pid = (TLeaf*) tree->GetLeaf("ReconstructedParticles","ReconstructedParticles.pid");
      rec_se = (TLeaf*) tree->GetLeaf("InclusiveKinematicsJB","InclusiveKinematicsJB.scatID.value");
  
      int nevents = tree->GetEntries();
      cout<<"Beginning Analysis of file "<<chEl->GetTitle()<<"!"<<endl;
      cout<<"Total number of events = "<<nevents<<endl;
  
      // Loop over all events
      for(int iEvent=0;iEvent<nevents;iEvent++){
        //if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
        tree->GetEntry(iEvent);

        //Scattered Electron variables
        auto found_se = false;
        TLorentzVector e_lab;
        TLorentzVector e_boosted;

        //Sums in colinear frame
        double pxsum = 0;
        double pysum = 0;
        double pzsum = 0;
        double Esum = 0;
  
        //Reconstructed Particles
        for(int irec = 0; irec < rec_px->GetLen(); irec++){

          if( irec == (int) rec_se->GetValue() ){
            found_se = true;
            e_lab.SetXYZT(rec_px->GetValue(irec),rec_py->GetValue(irec),rec_pz->GetValue(irec),rec_E->GetValue(irec));
            e_boosted = apply_boost(ei_const,pi_const,e_lab);
          }
          //else if(rec_charge->GetValue(irec)!=0){ //only use charged particles
          //else if(rec_pid->GetValue(irec)!=0){ //only use particles with non-zero pid
          else{
            //Lab frame quantities
            double px_lab = rec_px->GetValue(irec);
            double py_lab = rec_py->GetValue(irec);
            double pz_lab = rec_pz->GetValue(irec);
            double E_lab = rec_E->GetValue(irec);

            TLorentzVector hf_lab(px_lab,py_lab,pz_lab,E_lab);

            //Boost to colinear frame using first-order matrix
            TLorentzVector hf_boosted = apply_boost(ei_const,pi_const,hf_lab);
            pxsum += hf_boosted.Px();
            pysum += hf_boosted.Py();
            pzsum += hf_boosted.Pz();
            Esum += hf_boosted.E();
          }
        }//Finished particle Loop

        // DIS kinematics calculations

        //Scattered Electron method
        if(found_se){
          auto q_se = ei_const - e_lab;
          auto y_se = (q_se*pi_const)/(ei_const*pi_const);
          for(int ibin=0;ibin<nybins;ibin++){
            if( y_true->GetValue()>y_low[ibin] && y_true->GetValue()<y_hi[ibin] ){
              h1a[ibin]->Fill(100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );
            }
          }
        
          if(Q2_true->GetValue()>1 && Q2_true->GetValue()<10 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof1->Fill(y_true->GetValue(),100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin1[ibin]->Fill(100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>10 && Q2_true->GetValue()<100 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof1a->Fill(y_true->GetValue(),100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin1a[ibin]->Fill(100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>100 && Q2_true->GetValue()<1000 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof1b->Fill(y_true->GetValue(),100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin1b[ibin]->Fill(100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }

          if( fabs(100.*( y_se - y_true->GetValue() ) / y_true->GetValue()) < res_limit && Q2_true->GetValue()>1 && y_true->GetValue()>y_min && y_true->GetValue()<y_max)
            profh1->Fill(x_true->GetValue(),Q2_true->GetValue(),100.*( y_se - y_true->GetValue() ) / y_true->GetValue() );

        }

        //JB Method
        auto sigma_h = Esum - pzsum;
        auto ptsum = sqrt(pxsum*pxsum + pysum*pysum);

        if(sigma_h>0){
          auto y_jb = sigma_h / (2.*ei_const.E());

          for(int ibin=0;ibin<nybins;ibin++){
            if( y_true->GetValue()>y_low[ibin] && y_true->GetValue()<y_hi[ibin] ){
              h1b[ibin]->Fill(100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );
            }
          }

          if(Q2_true->GetValue()>1 && Q2_true->GetValue()<10 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof2->Fill(y_true->GetValue(),100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin2[ibin]->Fill(100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>10 && Q2_true->GetValue()<100 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof2a->Fill(y_true->GetValue(),100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin2a[ibin]->Fill(100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>100 && Q2_true->GetValue()<1000 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof2b->Fill(y_true->GetValue(),100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin2b[ibin]->Fill(100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }

          if( fabs(100.*( y_jb - y_true->GetValue() ) / y_true->GetValue()) < res_limit && Q2_true->GetValue()>1 && y_true->GetValue()>y_min && y_true->GetValue()<y_max )
            profh2->Fill(x_true->GetValue(),Q2_true->GetValue(),100.*( y_jb - y_true->GetValue() ) / y_true->GetValue() );

        }

        //DA method
        auto theta_e = e_boosted.Theta();
        auto theta_h = 2.*atan(sigma_h/ptsum);

        if(found_se && sigma_h>0){

          auto y_da = tan(theta_h/2.) / ( tan(theta_e/2.) + tan(theta_h/2.) );
          auto y_sigma = sigma_h / ( sigma_h + e_boosted.E()*(1.-cos(theta_e)) );
          auto y_esigma = (2.*ei_const.E()) / ( sigma_h + e_boosted.E()*(1.-cos(theta_e)) ) * y_sigma;

          for(int ibin=0;ibin<nybins;ibin++){
            if( y_true->GetValue()>y_low[ibin] && y_true->GetValue()<y_hi[ibin] ){
              h1c[ibin]->Fill(100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
              h1d[ibin]->Fill(100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
              h1e[ibin]->Fill(100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );
            }
          }

          if(Q2_true->GetValue()>1 && Q2_true->GetValue()<10 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof3->Fill(y_true->GetValue(),100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
            prof4->Fill(y_true->GetValue(),100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
            prof5->Fill(y_true->GetValue(),100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin3[ibin]->Fill(100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
                hbin4[ibin]->Fill(100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
                hbin5[ibin]->Fill(100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>10 && Q2_true->GetValue()<100 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof3a->Fill(y_true->GetValue(),100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
            prof4a->Fill(y_true->GetValue(),100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
            prof5a->Fill(y_true->GetValue(),100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin3a[ibin]->Fill(100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
                hbin4a[ibin]->Fill(100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
                hbin5a[ibin]->Fill(100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }
          if(Q2_true->GetValue()>100 && Q2_true->GetValue()<1000 && y_true->GetValue()>y_min && y_true->GetValue()<y_max){
            prof3b->Fill(y_true->GetValue(),100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
            prof4b->Fill(y_true->GetValue(),100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
            prof5b->Fill(y_true->GetValue(),100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );

            for(int ibin=0;ibin<30;ibin++){
              if(y_true->GetValue()>fine_1d_left[ibin] && y_true->GetValue()<(fine_1d_left[ibin] + fine_1d_width) ){
                hbin3b[ibin]->Fill(100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
                hbin4b[ibin]->Fill(100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
                hbin5b[ibin]->Fill(100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );
                break;
              }
            }
          }

          if( fabs(100.*( y_da - y_true->GetValue() ) / y_true->GetValue()) < res_limit && Q2_true->GetValue()>1 && y_true->GetValue()>y_min && y_true->GetValue()<y_max )
            profh3->Fill(x_true->GetValue(),Q2_true->GetValue(),100.*( y_da - y_true->GetValue() ) / y_true->GetValue() );
          if( fabs(100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue()) < res_limit && Q2_true->GetValue()>1 && y_true->GetValue()>y_min && y_true->GetValue()<y_max )
            profh4->Fill(x_true->GetValue(),Q2_true->GetValue(),100.*( y_sigma - y_true->GetValue() ) / y_true->GetValue() );
          if( fabs(100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue()) < res_limit && Q2_true->GetValue()>1 && y_true->GetValue()>y_min && y_true->GetValue()<y_max )
            profh5->Fill(x_true->GetValue(),Q2_true->GetValue(),100.*( y_esigma - y_true->GetValue() ) / y_true->GetValue() );

        }

      }//Finished Event Loop

      f->Close();
    }//Finished Files Loop within Chain
    chain->Reset();
    delete chain;
  }//Finised Chain Loop


  //Make Latex
  TPaveText* tex_energy = new TPaveText(0.1,0.7,0.9,0.9,"NDCNB");
  tex_energy->AddText("18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
  //tex_energy->AddText("Q^{2} > 1 GeV^{2}");
  tex_energy->AddText("Mix of Q^{2}");
	tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(12);

  TPaveText *tex1_1 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_1->AddText("Electron Method");
  tex1_1->SetFillStyle(4000);tex1_1->SetTextFont(63);tex1_1->SetTextSize(12);

  TPaveText *tex1_2 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_2->AddText("JB Method");
  tex1_2->SetFillStyle(4000);tex1_2->SetTextFont(63);tex1_2->SetTextSize(12);

  TPaveText *tex1_3 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_3->AddText("DA Method");
  tex1_3->SetFillStyle(4000);tex1_3->SetTextFont(63);tex1_3->SetTextSize(12);

  TPaveText *tex1_4 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_4->AddText("#Sigma Method");
  tex1_4->SetFillStyle(4000);tex1_4->SetTextFont(63);tex1_4->SetTextSize(12);

  TPaveText *tex1_5 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_5->AddText("e#Sigma Method");
  tex1_5->SetFillStyle(4000);tex1_5->SetTextFont(63);tex1_5->SetTextSize(12);

  //Make Plots (Constant beam momenta)
  //----------------------------------
  TCanvas *c1a = new TCanvas("c1a");
  c1a->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1a->cd(ibin+1);
    h1a[ibin]->Draw();

    //cout<<h1a[ibin]->GetRMS()<<endl;
  }
  c1a->cd(6);
  tex_energy->Draw();tex1_1->Draw();

  TCanvas *c1b = new TCanvas("c1b");
  c1b->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1b->cd(ibin+1);
    h1b[ibin]->Draw();

    //cout<<h1b[ibin]->GetRMS()<<endl;
  }
  c1b->cd(6);
  tex_energy->Draw();tex1_2->Draw();

  TCanvas *c1c = new TCanvas("c1c");
  c1c->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1c->cd(ibin+1);
    h1c[ibin]->Draw();

    //cout<<h1c[ibin]->GetRMS()<<endl;
  }
  c1c->cd(6);
  tex_energy->Draw();tex1_3->Draw();

  TCanvas *c1d = new TCanvas("c1d");
  c1d->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1d->cd(ibin+1);
    h1d[ibin]->Draw();

    //cout<<h1c[ibin]->GetRMS()<<endl;
  }
  c1d->cd(6);
  tex_energy->Draw();tex1_4->Draw();

  TCanvas *c1e = new TCanvas("c1e");
  c1e->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1e->cd(ibin+1);
    h1e[ibin]->Draw();

    //cout<<h1c[ibin]->GetRMS()<<endl;
  }
  c1e->cd(6);
  tex_energy->Draw();tex1_5->Draw();

  //1D plots vs y
  //------
  //Extract Resolutions from TProfile
  TGraph *gr1 = new TGraph();
  gr1->SetMarkerColor(kRed);gr1->SetLineColor(kRed);gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.5);gr1->SetLineWidth(2);
  TGraph *gr2 = new TGraph();
  gr2->SetMarkerColor(kBlue);gr2->SetLineColor(kBlue);gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.5);gr2->SetLineWidth(2);
  TGraph *gr3 = new TGraph();
  gr3->SetMarkerColor(kGreen);gr3->SetLineColor(kGreen);gr3->SetMarkerStyle(22);
  gr3->SetMarkerSize(1.5);gr3->SetLineWidth(2);
  TGraph *gr4 = new TGraph();
  gr4->SetMarkerColor(kBlack);gr4->SetLineColor(kBlack);gr4->SetMarkerStyle(22);
  gr4->SetMarkerSize(1.5);gr4->SetLineWidth(2);
  TGraph *gr5 = new TGraph();
  gr5->SetMarkerColor(kOrange);gr5->SetLineColor(kOrange);gr5->SetMarkerStyle(22);
  gr5->SetMarkerSize(1.5);gr5->SetLineWidth(2);

  TGraph *gr1a = new TGraph();
  gr1a->SetMarkerColor(kRed);gr1a->SetLineColor(kRed);gr1a->SetMarkerStyle(20);
  gr1a->SetMarkerSize(1.5);gr1a->SetLineWidth(2);
  TGraph *gr2a = new TGraph();
  gr2a->SetMarkerColor(kBlue);gr2a->SetLineColor(kBlue);gr2a->SetMarkerStyle(21);
  gr2a->SetMarkerSize(1.5);gr2a->SetLineWidth(2);
  TGraph *gr3a = new TGraph();
  gr3a->SetMarkerColor(kGreen);gr3a->SetLineColor(kGreen);gr3a->SetMarkerStyle(22);
  gr3a->SetMarkerSize(1.5);gr3a->SetLineWidth(2);
  TGraph *gr4a = new TGraph();
  gr4a->SetMarkerColor(kBlack);gr4a->SetLineColor(kBlack);gr4a->SetMarkerStyle(22);
  gr4a->SetMarkerSize(1.5);gr4a->SetLineWidth(2);
  TGraph *gr5a = new TGraph();
  gr5a->SetMarkerColor(kOrange);gr5a->SetLineColor(kOrange);gr5a->SetMarkerStyle(22);
  gr5a->SetMarkerSize(1.5);gr5a->SetLineWidth(2);

  TGraph *gr1b = new TGraph();
  gr1b->SetMarkerColor(kRed);gr1b->SetLineColor(kRed);gr1b->SetMarkerStyle(20);
  gr1b->SetMarkerSize(1.5);gr1b->SetLineWidth(2);
  TGraph *gr2b = new TGraph();
  gr2b->SetMarkerColor(kBlue);gr2b->SetLineColor(kBlue);gr2b->SetMarkerStyle(21);
  gr2b->SetMarkerSize(1.5);gr2b->SetLineWidth(2);
  TGraph *gr3b = new TGraph();
  gr3b->SetMarkerColor(kGreen);gr3b->SetLineColor(kGreen);gr3b->SetMarkerStyle(22);
  gr3b->SetMarkerSize(1.5);gr3b->SetLineWidth(2);
  TGraph *gr4b = new TGraph();
  gr4b->SetMarkerColor(kBlack);gr4b->SetLineColor(kBlack);gr4b->SetMarkerStyle(22);
  gr4b->SetMarkerSize(1.5);gr4b->SetLineWidth(2);
  TGraph *gr5b = new TGraph();
  gr5b->SetMarkerColor(kOrange);gr5b->SetLineColor(kOrange);gr5b->SetMarkerStyle(22);
  gr5b->SetMarkerSize(1.5);gr5b->SetLineWidth(2);

  int counter1(0),counter1a(0),counter1b(0);
  int counter2(0),counter2a(0),counter2b(0);
  int counter3(0),counter3a(0),counter3b(0);
  int counter4(0),counter4a(0),counter4b(0);
  int counter5(0),counter5a(0),counter5b(0);

  for(int ibin=1;ibin<=prof1->GetNbinsX();ibin++){

    if(prof1->GetBinError(ibin)>0){
      gr1->SetPoint(counter1,prof1->GetBinCenter(ibin),prof1->GetBinError(ibin));
      counter1++;
    }
    
    if(prof2->GetBinError(ibin)>0){
      gr2->SetPoint(counter2,prof2->GetBinCenter(ibin),prof2->GetBinError(ibin));
      counter2++;
    }
    
    if(prof3->GetBinError(ibin)>0){
      gr3->SetPoint(counter3,prof3->GetBinCenter(ibin),prof3->GetBinError(ibin));
      counter3++;
    }

    if(prof4->GetBinError(ibin)>0){
      gr4->SetPoint(counter4,prof4->GetBinCenter(ibin),prof4->GetBinError(ibin));
      counter4++;
    }

    if(prof5->GetBinError(ibin)>0){
      gr5->SetPoint(counter5,prof5->GetBinCenter(ibin),prof5->GetBinError(ibin));
      counter5++;
    }

    if(prof1a->GetBinError(ibin)>0){
      gr1a->SetPoint(counter1a,prof1a->GetBinCenter(ibin),prof1a->GetBinError(ibin));
      counter1a++;
    }
    
    if(prof2a->GetBinError(ibin)>0){
      gr2a->SetPoint(counter2a,prof2a->GetBinCenter(ibin),prof2a->GetBinError(ibin));
      counter2a++;
    }
    
    if(prof3a->GetBinError(ibin)>0){
      gr3a->SetPoint(counter3a,prof3a->GetBinCenter(ibin),prof3a->GetBinError(ibin));
      counter3a++;
    }

    if(prof4a->GetBinError(ibin)>0){
      gr4a->SetPoint(counter4a,prof4a->GetBinCenter(ibin),prof4a->GetBinError(ibin));
      counter4a++;
    }

    if(prof5a->GetBinError(ibin)>0){
      gr5a->SetPoint(counter5a,prof5a->GetBinCenter(ibin),prof5a->GetBinError(ibin));
      counter5a++;
    }

    if(prof1b->GetBinError(ibin)>0){
      gr1b->SetPoint(counter1b,prof1b->GetBinCenter(ibin),prof1b->GetBinError(ibin));
      counter1b++;
    }
    
    if(prof2b->GetBinError(ibin)>0){
      gr2b->SetPoint(counter2b,prof2b->GetBinCenter(ibin),prof2b->GetBinError(ibin));
      counter2b++;
    }
    
    if(prof3b->GetBinError(ibin)>0){
      gr3b->SetPoint(counter3b,prof3b->GetBinCenter(ibin),prof3b->GetBinError(ibin));
      counter3b++;
    }

    if(prof4b->GetBinError(ibin)>0){
      gr4b->SetPoint(counter4b,prof4b->GetBinCenter(ibin),prof4b->GetBinError(ibin));
      counter4b++;
    }

    if(prof5b->GetBinError(ibin)>0){
      gr5b->SetPoint(counter5b,prof5b->GetBinCenter(ibin),prof5b->GetBinError(ibin));
      counter5b++;
    }
 
  }

  //Make plot
  PadMxN *pad2x2_1 = new PadMxN("c2",500,500,150,150,100,125,2,2);
  pad2x2_1->Draw();
  TPad *mypad = {0};
	TH1 *hframe_1 = new TH2F("hframe_1","",100,0.01,0.99,100,0.01,55);

  for(int iCan=0;iCan<4;iCan++){
    mypad = pad2x2_1->GetPad(iCan+1);
    //gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_1->Draw();
    hframe_1->GetXaxis()->SetTitle("y");hframe_1->GetYaxis()->SetTitle("#sigma_{y}/y [%]");
		hframe_1->GetXaxis()->SetLabelFont(63);hframe_1->GetYaxis()->SetLabelFont(63);
    hframe_1->GetXaxis()->SetLabelSize(25);hframe_1->GetYaxis()->SetLabelSize(25);
    hframe_1->GetXaxis()->SetLabelOffset(0.01);hframe_1->GetYaxis()->SetLabelOffset(0.01);
    hframe_1->GetXaxis()->CenterTitle(1);hframe_1->GetYaxis()->CenterTitle(1);
    hframe_1->GetXaxis()->SetTitleSize(40);hframe_1->GetXaxis()->SetTitleOffset(2.5); 
    hframe_1->GetYaxis()->SetTitleSize(40);hframe_1->GetYaxis()->SetTitleOffset(3.0);

    if(iCan==0){ 
      gr1->Draw("PL same");gr2->Draw("PL same");gr3->Draw("PL same");gr4->Draw("PL same");gr5->Draw("PL same");
    }
    else if(iCan==1){
      gr1a->Draw("PL same");gr2a->Draw("PL same");gr3a->Draw("PL same");gr4a->Draw("PL same");gr5a->Draw("PL same");
    }
    else if(iCan==2){
      gr1b->Draw("PL same");gr2b->Draw("PL same");gr3b->Draw("PL same");gr4b->Draw("PL same");gr5b->Draw("PL same");
    }
    
  }

  // get the text pad
  pad2x2_1->GetPad(5);

  TLegend *leg1 = new TLegend(0.525,0.2,0.825,0.45);
  leg1->SetBorderSize(0);leg1->SetTextSize(0.03);
  leg1->AddEntry((TObject*)0, "18 GeV e^{-} on 275 GeV p", "");
  leg1->AddEntry(gr1,"Electron Method","p");
  leg1->AddEntry(gr2,"JB Method","p");
  leg1->AddEntry(gr3,"DA Method","p");
  leg1->AddEntry(gr4,"#Sigma Method","p");
  leg1->AddEntry(gr5,"e#Sigma Method","p");
  leg1->Draw();

  TPaveText* pave1 = new TPaveText(0.15,0.8,0.375,0.9,"NDCNB");
  pave1->AddText("1 < Q^{2} < 10 GeV^{2}");
	pave1->SetFillStyle(4000);
  pave1->SetTextFont(63);pave1->SetTextSize(30);
  pave1->SetTextColor(kBlack);
  pave1->Draw();

  TPaveText* pave2 = new TPaveText(0.5,0.8,0.8,0.9,"NDCNB");
  pave2->AddText("10 < Q^{2} < 100 GeV^{2}");
	pave2->SetFillStyle(4000);
  pave2->SetTextFont(63);pave2->SetTextSize(30);
  pave2->SetTextColor(kBlack);
  pave2->Draw();

  TPaveText* pave3 = new TPaveText(0.15,0.35,0.375,0.55,"NDCNB");
  pave3->AddText("100 < Q^{2} < 1000 GeV^{2}");
	pave3->SetFillStyle(4000);
  pave3->SetTextFont(63);pave3->SetTextSize(30);
  pave3->SetTextColor(kBlack);
  pave3->Draw();

  TPaveText* pave4 = new TPaveText(0.90,0.09,0.94,0.12,"NDCNB");
  pave4->AddText("1");
	pave4->SetFillStyle(4000);
  pave4->SetTextFont(63);pave4->SetTextSize(25);
  pave4->Draw();

  TPaveText* pave5 = new TPaveText(0.1,0.09,0.14,0.12,"NDCNB");
  pave5->AddText("0");
	pave5->SetFillStyle(4000);
  pave5->SetTextFont(63);pave5->SetTextSize(25);
  pave5->Draw();

  TPaveText* pave6 = new TPaveText(0.5,0.09,0.54,0.12,"NDCNB");
  pave6->AddText("0");
	pave6->SetFillStyle(4000);
  pave6->SetTextFont(63);pave6->SetTextSize(25);
  pave6->Draw();

  TPaveText* pave7 = new TPaveText(0.1,0.1,0.12,0.14,"NDCNB");
  pave7->AddText("0");
	pave7->SetFillStyle(4000);
  pave7->SetTextFont(63);pave7->SetTextSize(25);
  pave7->Draw();

  TPaveText* pave8 = new TPaveText(0.1,0.5,0.12,0.53,"NDCNB");
  pave8->AddText("0");
	pave8->SetFillStyle(4000);
  pave8->SetTextFont(63);pave8->SetTextSize(25);
  pave8->Draw();
  //------

  //------
  //Extract Resolutions from TH1s
  int counterf1(0),counterf1a(0),counterf1b(0);
  int counterf2(0),counterf2a(0),counterf2b(0);
  int counterf3(0),counterf3a(0),counterf3b(0);
  int counterf4(0),counterf4a(0),counterf4b(0);
  int counterf5(0),counterf5a(0),counterf5b(0);

  TGraph *grf1 = new TGraph();
  grf1->SetMarkerColor(kRed);grf1->SetLineColor(kRed);grf1->SetMarkerStyle(20);
  grf1->SetMarkerSize(1.5);grf1->SetLineWidth(2);

  TGraph *grf1a = new TGraph();
  grf1a->SetMarkerColor(kRed);grf1a->SetLineColor(kRed);grf1a->SetMarkerStyle(20);
  grf1a->SetMarkerSize(1.5);grf1a->SetLineWidth(2);

  TGraph *grf1b = new TGraph();
  grf1b->SetMarkerColor(kRed);grf1b->SetLineColor(kRed);grf1b->SetMarkerStyle(20);
  grf1b->SetMarkerSize(1.5);grf1b->SetLineWidth(2);

  TGraph *grf2 = new TGraph();
  grf2->SetMarkerColor(kBlue);grf2->SetLineColor(kBlue);grf2->SetMarkerStyle(21);
  grf2->SetMarkerSize(1.5);grf2->SetLineWidth(2);

  TGraph *grf2a = new TGraph();
  grf2a->SetMarkerColor(kBlue);grf2a->SetLineColor(kBlue);grf2a->SetMarkerStyle(21);
  grf2a->SetMarkerSize(1.5);grf2a->SetLineWidth(2);

  TGraph *grf2b = new TGraph();
  grf2b->SetMarkerColor(kBlue);grf2b->SetLineColor(kBlue);grf2b->SetMarkerStyle(21);
  grf2b->SetMarkerSize(1.5);grf2b->SetLineWidth(2);

  TGraph *grf3 = new TGraph();
  grf3->SetMarkerColor(kGreen);grf3->SetLineColor(kGreen);grf3->SetMarkerStyle(22);
  grf3->SetMarkerSize(1.5);grf3->SetLineWidth(2);

  TGraph *grf3a = new TGraph();
  grf3a->SetMarkerColor(kGreen);grf3a->SetLineColor(kGreen);grf3a->SetMarkerStyle(22);
  grf3a->SetMarkerSize(1.5);grf3a->SetLineWidth(2);

  TGraph *grf3b = new TGraph();
  grf3b->SetMarkerColor(kGreen);grf3b->SetLineColor(kGreen);grf3b->SetMarkerStyle(22);
  grf3b->SetMarkerSize(1.5);grf3b->SetLineWidth(2);

  TGraph *grf4 = new TGraph();
  grf4->SetMarkerColor(kBlack);grf4->SetLineColor(kBlack);grf4->SetMarkerStyle(22);
  grf4->SetMarkerSize(1.5);grf4->SetLineWidth(2);

  TGraph *grf4a = new TGraph();
  grf4a->SetMarkerColor(kBlack);grf4a->SetLineColor(kBlack);grf4a->SetMarkerStyle(22);
  grf4a->SetMarkerSize(1.5);grf4a->SetLineWidth(2);

  TGraph *grf4b = new TGraph();
  grf4b->SetMarkerColor(kBlack);grf4b->SetLineColor(kBlack);grf4b->SetMarkerStyle(22);
  grf4b->SetMarkerSize(1.5);grf4b->SetLineWidth(2);

  TGraph *grf5 = new TGraph();
  grf5->SetMarkerColor(kOrange);grf5->SetLineColor(kOrange);grf5->SetMarkerStyle(22);
  grf5->SetMarkerSize(1.5);grf5->SetLineWidth(2);

  TGraph *grf5a = new TGraph();
  grf5a->SetMarkerColor(kOrange);grf5a->SetLineColor(kOrange);grf5a->SetMarkerStyle(22);
  grf5a->SetMarkerSize(1.5);grf5a->SetLineWidth(2);

  TGraph *grf5b = new TGraph();
  grf5b->SetMarkerColor(kOrange);grf5b->SetLineColor(kOrange);grf5b->SetMarkerStyle(22);
  grf5b->SetMarkerSize(1.5);grf5b->SetLineWidth(2);

  TCanvas *ctemp = new TCanvas("ctemp");

  for(int ibin=0;ibin<30;ibin++){

    auto ycenter = fine_1d_left[ibin] + fine_1d_width/2.;

    if(hbin1[ibin]->GetEntries()>50){
      //First fit
      hbin1[ibin]->Fit("gaus","","",hbin1[ibin]->GetBinCenter(hbin1[ibin]->GetMaximumBin()) - 50,
                                    hbin1[ibin]->GetBinCenter(hbin1[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin1[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin1[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin1[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf1->SetPoint(counterf1, ycenter , sigma );
      counterf1++;
    }
    if(hbin1a[ibin]->GetEntries()>50){
      //First fit
      hbin1a[ibin]->Fit("gaus","","",hbin1a[ibin]->GetBinCenter(hbin1a[ibin]->GetMaximumBin()) - 50,
                                     hbin1a[ibin]->GetBinCenter(hbin1a[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin1a[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin1a[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin1a[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf1a->SetPoint(counterf1a, ycenter , sigma );
      counterf1a++;
      
    }
    if(hbin1b[ibin]->GetEntries()>50){
      //First fit
      hbin1b[ibin]->Fit("gaus","","",hbin1b[ibin]->GetBinCenter(hbin1b[ibin]->GetMaximumBin()) - 50,
                                     hbin1b[ibin]->GetBinCenter(hbin1b[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin1b[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin1b[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin1b[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf1b->SetPoint(counterf1b, ycenter , sigma );
      counterf1b++;
    }

    if(hbin2[ibin]->GetEntries()>50){
      //First fit
      hbin2[ibin]->Fit("gaus","","",hbin2[ibin]->GetBinCenter(hbin2[ibin]->GetMaximumBin()) - 50,
                                    hbin2[ibin]->GetBinCenter(hbin2[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin2[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin2[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin2[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf2->SetPoint(counterf2, ycenter , sigma );
      counterf2++;
    }
    if(hbin2a[ibin]->GetEntries()>50){
      //First fit
      hbin2a[ibin]->Fit("gaus","","",hbin2a[ibin]->GetBinCenter(hbin2a[ibin]->GetMaximumBin()) - 50,
                                     hbin2a[ibin]->GetBinCenter(hbin2a[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin2a[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin2a[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin2a[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf2a->SetPoint(counterf2a, ycenter , sigma );
      counterf2a++;
      
    }
    if(hbin2b[ibin]->GetEntries()>50){
      //First fit
      hbin2b[ibin]->Fit("gaus","","",hbin2b[ibin]->GetBinCenter(hbin2b[ibin]->GetMaximumBin()) - 50,
                                     hbin2b[ibin]->GetBinCenter(hbin2b[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin2b[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin2b[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin2b[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf2b->SetPoint(counterf2b, ycenter , sigma );
      counterf2b++;
    }

    if(hbin3[ibin]->GetEntries()>50){
      //First fit
      hbin3[ibin]->Fit("gaus","","",hbin3[ibin]->GetBinCenter(hbin3[ibin]->GetMaximumBin()) - 50,
                                    hbin3[ibin]->GetBinCenter(hbin3[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin3[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin3[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin3[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf3->SetPoint(counterf3, ycenter , sigma );
      counterf3++;
    }
    if(hbin3a[ibin]->GetEntries()>50){
      //First fit
      hbin3a[ibin]->Fit("gaus","","",hbin3a[ibin]->GetBinCenter(hbin3a[ibin]->GetMaximumBin()) - 50,
                                     hbin3a[ibin]->GetBinCenter(hbin3a[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin3a[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin3a[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin3a[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf3a->SetPoint(counterf3a, ycenter , sigma );
      counterf3a++;
      
    }
    if(hbin3b[ibin]->GetEntries()>50){
      //First fit
      hbin3b[ibin]->Fit("gaus","","",hbin3b[ibin]->GetBinCenter(hbin3b[ibin]->GetMaximumBin()) - 50,
                                     hbin3b[ibin]->GetBinCenter(hbin3b[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin3b[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin3b[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin3b[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf3b->SetPoint(counterf3b, ycenter , sigma );
      counterf3b++;
    }

    if(hbin4[ibin]->GetEntries()>50){
      //First fit
      hbin4[ibin]->Fit("gaus","","",hbin4[ibin]->GetBinCenter(hbin4[ibin]->GetMaximumBin()) - 50,
                                    hbin4[ibin]->GetBinCenter(hbin4[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin4[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin4[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin4[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf4->SetPoint(counterf4, ycenter , sigma );
      counterf4++;
    }
    if(hbin4a[ibin]->GetEntries()>50){
      //First fit
      hbin4a[ibin]->Fit("gaus","","",hbin4a[ibin]->GetBinCenter(hbin4a[ibin]->GetMaximumBin()) - 50,
                                     hbin4a[ibin]->GetBinCenter(hbin4a[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin4a[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin4a[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin4a[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf4a->SetPoint(counterf4a, ycenter , sigma );
      counterf4a++;
      
    }
    if(hbin4b[ibin]->GetEntries()>50){
      //First fit
      hbin4b[ibin]->Fit("gaus","","",hbin4b[ibin]->GetBinCenter(hbin4b[ibin]->GetMaximumBin()) - 50,
                                     hbin4b[ibin]->GetBinCenter(hbin4b[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin4b[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin4b[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin4b[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf4b->SetPoint(counterf4b, ycenter , sigma );
      counterf4b++;
    }

    if(hbin5[ibin]->GetEntries()>50){
      //First fit
      hbin5[ibin]->Fit("gaus","","",hbin5[ibin]->GetBinCenter(hbin5[ibin]->GetMaximumBin()) - 50,
                                    hbin5[ibin]->GetBinCenter(hbin5[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin5[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin5[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin5[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf5->SetPoint(counterf5, ycenter , sigma );
      counterf5++;
    }
    if(hbin5a[ibin]->GetEntries()>50){
      //First fit
      hbin5a[ibin]->Fit("gaus","","",hbin5a[ibin]->GetBinCenter(hbin5a[ibin]->GetMaximumBin()) - 50,
                                     hbin5a[ibin]->GetBinCenter(hbin5a[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin5a[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin5a[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin5a[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf5a->SetPoint(counterf5a, ycenter , sigma );
      counterf5a++;
      
    }
    if(hbin5b[ibin]->GetEntries()>50){
      //First fit
      hbin5b[ibin]->Fit("gaus","","",hbin5b[ibin]->GetBinCenter(hbin5b[ibin]->GetMaximumBin()) - 50,
                                     hbin5b[ibin]->GetBinCenter(hbin5b[ibin]->GetMaximumBin()) + 50);

      //Second fit
      TF1 *func = (TF1*)hbin5b[ibin]->GetListOfFunctions()->FindObject("gaus");
      hbin5b[ibin]->Fit("gaus","","",func->GetParameter(1)+2*func->GetParameter(2),func->GetParameter(1)-2*func->GetParameter(2));
      func = (TF1*)hbin5b[ibin]->GetListOfFunctions()->FindObject("gaus");
      auto mean = func->GetParameter(1);
      auto sigma = func->GetParameter(2);

      grf5b->SetPoint(counterf5b, ycenter , sigma );
      counterf5b++;
    }

  }

  //Make plot
  PadMxN *pad2x2f = new PadMxN("c2f",500,500,150,150,100,125,2,2);
  pad2x2f->Draw();

  for(int iCan=0;iCan<4;iCan++){
    mypad = pad2x2f->GetPad(iCan+1);
    //gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_1->Draw();

    if(iCan==0){ 
      grf1->Draw("PL same");grf2->Draw("PL same");grf3->Draw("PL same");grf4->Draw("PL same");grf5->Draw("PL same");
    }
    else if(iCan==1){
      grf1a->Draw("PL same");grf2a->Draw("PL same");grf3a->Draw("PL same");grf4a->Draw("PL same");grf5a->Draw("PL same");
    }
    else if(iCan==2){
      grf1b->Draw("PL same");grf2b->Draw("PL same");grf3b->Draw("PL same");grf4b->Draw("PL same");grf5b->Draw("PL same");
    }
    
  }

  // get the text pad
  pad2x2f->GetPad(5);
  leg1->Draw();
  pave1->Draw();
  pave2->Draw();
  pave3->Draw();
  pave4->Draw();
  pave5->Draw();
  pave6->Draw();
  pave7->Draw();
  pave8->Draw();
  //------

  //2D histogram plots
  for(int ibinx=1;ibinx<=profh1->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=profh1->GetNbinsY();ibiny++){
      
      hh1->SetBinContent(ibinx,ibiny,profh1->GetBinError(ibinx,ibiny) );
      hh2->SetBinContent(ibinx,ibiny,profh2->GetBinError(ibinx,ibiny) );
      hh3->SetBinContent(ibinx,ibiny,profh3->GetBinError(ibinx,ibiny) );
      hh4->SetBinContent(ibinx,ibiny,profh4->GetBinError(ibinx,ibiny) );
      hh5->SetBinContent(ibinx,ibiny,profh5->GetBinError(ibinx,ibiny) );

      cout<<hh1->GetXaxis()->GetBinCenter(ibinx)<<"   "<<hh1->GetYaxis()->GetBinCenter(ibiny)<<"   "
          <<hh1->GetBinContent(ibinx,ibiny)<<"   "<<hh2->GetBinContent(ibinx,ibiny)<<"   "<<hh3->GetBinContent(ibinx,ibiny)
          <<"   "<<hh4->GetBinContent(ibinx,ibiny)<<"   "<<hh5->GetBinContent(ibinx,ibiny)
          <<endl;

    }
  }

  PadMxN *pad3x1 = new PadMxN("c3",500,500,175,50,175,50,3,1);
  pad3x1->Draw();

  const Double_t min = 0.;
  const Double_t max = 55.;
  
  const Int_t nLevels = 999;
  Double_t levels[nLevels];
  gStyle->SetNumberContours(nLevels);

  for(int i = 1; i < nLevels; i++) {
    levels[i] = min + (max - min) / (nLevels - 1) * (i);
  }
  //levels[0] = 0.01;
  levels[0] = -1; //This also works as we want

  for(int iCan=0;iCan<3;iCan++){
    mypad = pad3x1->GetPad(iCan+1);
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    if(iCan==0){
      hh1->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh1->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh1->Draw("col");
    }
    else if(iCan==1){
      hh2->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh2->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh2->Draw("col");
    }
    else if(iCan==2){
      hh3->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh3->GetZaxis()->SetRangeUser(min, max); // ... set the range ..
      hh3->Draw("colz");
    }
  }

  // get the text pad
  pad3x1->GetPad(4);

  TPaveText* pave9 = new TPaveText(0.45,0.0,0.55,0.125,"NDCNB");
  pave9->AddText("x");
	pave9->SetFillStyle(4000);
  pave9->SetTextFont(63);pave9->SetTextSize(45);
  pave9->Draw();

  TPaveText* pave10 = new TPaveText(0.0,0.45,0.1,0.55,"NDCNB");
  auto text10 = pave10->AddText("Q^{2} [GeV^{2}]");
	pave10->SetFillStyle(4000);
  pave10->SetTextFont(63);pave10->SetTextSize(45);
  text10->SetTextAngle(90.);
  pave10->Draw();

  TPaveText* pave11 = new TPaveText(0.92,0.70,1.0,0.85,"NDCNB");
  auto text11 = pave11->AddText("#sigma_{y}/y [%]");
	pave11->SetFillStyle(4000);
  pave11->SetTextFont(63);pave11->SetTextSize(35);
  text11->SetTextAngle(90.);
  pave11->Draw();

  TPaveText* pave12 = new TPaveText(0.075,0.65,0.25,0.85,"NDCNB");
  pave12->AddText("Electron Method");
	pave12->SetFillStyle(4000);
  pave12->SetTextFont(42);pave12->SetTextSize(0.055);
  pave12->Draw();

  TPaveText* pave13 = new TPaveText(0.335,0.65,0.485,0.85,"NDCNB");
  pave13->AddText("JB Method");
	pave13->SetFillStyle(4000);
  pave13->SetTextFont(42);pave13->SetTextSize(0.055);
  pave13->Draw();
  
  TPaveText* pave14 = new TPaveText(0.615,0.65,0.765,0.85,"NDCNB");
  pave14->AddText("DA Method");
	pave14->SetFillStyle(4000);
  pave14->SetTextFont(42);pave14->SetTextSize(0.055);
  pave14->Draw();

  TPaveText* pave15 = new TPaveText(0.11,0.80,0.26,0.9,"NDCNB");
  pave15->AddText("18 GeV e^{-} on 275 GeV p");
	pave15->SetFillStyle(4000);
  pave15->SetTextFont(42);pave15->SetTextSize(0.055);
  pave15->Draw();

  TPaveText* pave16 = new TPaveText(0.3605,0.1,0.3705,0.155,"NDCNB");
  pave16->AddText("10^{-5}");
  pave16->SetFillStyle(1001);pave16->SetFillColor(0);
	//pave16->SetFillStyle(4000);
  pave16->SetTextFont(63);pave16->SetTextSize(25);
  pave16->Draw();

  TPaveText* pave17 = new TPaveText(0.631,0.1,0.641,0.155,"NDCNB");
  pave17->AddText("10^{-5}");
  pave17->SetFillStyle(1001);pave17->SetFillColor(kWhite);
	//pave17->SetFillStyle(4000);
  pave17->SetTextFont(63);pave17->SetTextSize(25);
  pave17->Draw();

  //gROOT->ProcessLine("c3->Modified();");
  //gROOT->ProcessLine("c3->Update();");
  //gROOT->ProcessLine("c3->SaveAs(\"mytest.C\")");

  PadMxN *pad3x1_1 = new PadMxN("c4",500,500,175,50,175,50,3,1);
  pad3x1_1->Draw();

  for(int iCan=0;iCan<3;iCan++){
    mypad = pad3x1_1->GetPad(iCan+1);
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    if(iCan==0){
      hh1->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh1->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh1->Draw("col");
    }
    else if(iCan==1){
      hh4->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh4->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh4->Draw("col");
    }
    else if(iCan==2){
      hh5->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh5->GetZaxis()->SetRangeUser(min, max); // ... set the range ..
      hh5->Draw("colz");
    }
  }

  // get the text pad
  pad3x1_1->GetPad(4);

  pave9->Draw();pave10->Draw();pave11->Draw();pave12->Draw();
  pave15->Draw();pave16->Draw();pave17->Draw();

  TPaveText* pave18 = new TPaveText(0.335,0.65,0.485,0.85,"NDCNB");
  pave18->AddText("#Sigma Method");
	pave18->SetFillStyle(4000);
  pave18->SetTextFont(42);pave18->SetTextSize(0.055);
  pave18->Draw();
  
  TPaveText* pave19 = new TPaveText(0.615,0.65,0.765,0.85,"NDCNB");
  pave19->AddText("e#Sigma Method");
	pave19->SetFillStyle(4000);
  pave19->SetTextFont(42);pave19->SetTextSize(0.055);
  pave19->Draw();

  //Print to file
  c1a->Print("./plots/athena_Rec_y.pdf[");
  c1a->Print("./plots/athena_Rec_y.pdf");
  c1b->Print("./plots/athena_Rec_y.pdf");
  c1c->Print("./plots/athena_Rec_y.pdf");
  c1d->Print("./plots/athena_Rec_y.pdf");
  c1e->Print("./plots/athena_Rec_y.pdf");
  gROOT->ProcessLine("c2->Print(\"./plots/athena_Rec_y.pdf\");");
  //gROOT->ProcessLine("c2->Print(\"./plots/athena_Rec_y.pdf]\");");
  gROOT->ProcessLine("c2f->Print(\"./plots/athena_Rec_y.pdf\");");
  gROOT->ProcessLine("c3->Print(\"./plots/athena_Rec_y.pdf\");"); //Crashes on my home computer...works on SDCC machine
  gROOT->ProcessLine("c4->Print(\"./plots/athena_Rec_y.pdf\");");
  gROOT->ProcessLine("c4->Print(\"./plots/athena_Rec_y.pdf]\");");

}
