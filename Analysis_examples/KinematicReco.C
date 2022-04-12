#include <TTreeReader.h>
#include <TTreeReaderArray.h>

TLorentzVector TransformLabToHeadOnFrame(TLorentzVector eBeam, TLorentzVector pBeam, TLorentzVector Lvec);

void KinematicReco(){
  // Random number generator for smearing
  TRandom3 *rndm = new TRandom3(7);
  
  TChain *chain = new TChain("events");
  chain->Add("PythiaFiles_22_11_21/*.root");

  
  // Truth
  TTreeReader tr(chain);
  TTreeReaderArray<Int_t>    mcparticles_pdgID(tr,     "mcparticles.pdgID");
  TTreeReaderArray<Double_t> mcparticles_psx(tr,       "mcparticles.ps.x");
  TTreeReaderArray<Double_t> mcparticles_psy(tr,       "mcparticles.ps.y");
  TTreeReaderArray<Double_t> mcparticles_psz(tr,       "mcparticles.ps.z");
  TTreeReaderArray<Int_t>    mcparticles_status(tr,    "mcparticles.status");
  TTreeReaderArray<Int_t>    mcparticles_genStatus(tr, "mcparticles.genStatus");
  TTreeReaderArray<Double_t> mcparticles_mass(tr,      "mcparticles.mass");
  TTreeReaderArray<Int_t> mcparticles_charge(tr,      "mcparticles.charge");

  // Generated
  TTreeReaderArray<int> GeneratedParticles_pid(tr, "GeneratedParticles.pid");
  TTreeReaderArray<float> GeneratedParticles_energy(tr, "GeneratedParticles.energy");
  TTreeReaderArray<float> GeneratedParticles_p_x(tr, "GeneratedParticles.p.x");
  TTreeReaderArray<float> GeneratedParticles_p_y(tr, "GeneratedParticles.p.y");
  TTreeReaderArray<float> GeneratedParticles_p_z(tr, "GeneratedParticles.p.z");
  TTreeReaderArray<float> GeneratedParticles_mass(tr, "GeneratedParticles.mass");

  // Reco
  TTreeReaderArray<int> ReconstructedParticles_pid(tr, "ReconstructedParticles.pid");
  TTreeReaderArray<float> ReconstructedParticles_energy(tr, "ReconstructedParticles.energy");
  TTreeReaderArray<float> ReconstructedParticles_p_x(tr, "ReconstructedParticles.p.x");
  TTreeReaderArray<float> ReconstructedParticles_p_y(tr, "ReconstructedParticles.p.y");
  TTreeReaderArray<float> ReconstructedParticles_p_z(tr, "ReconstructedParticles.p.z");
  TTreeReaderArray<float> ReconstructedParticles_mass(tr, "ReconstructedParticles.mass");
  TTreeReaderArray<Int_t> ReconstructedParticles_mcID_value(tr, "ReconstructedParticles.mcID.value");
  TTreeReaderArray<Short_t> ReconstructedParticles_charge(tr, "ReconstructedParticles.charge");

  // Reco
  TTreeReaderArray<int> ReconstructedChargedParticles_pid(tr, "ReconstructedChargedParticles.pid");
  TTreeReaderArray<float> ReconstructedChargedParticles_energy(tr, "ReconstructedChargedParticles.energy");
  TTreeReaderArray<float> ReconstructedChargedParticles_p_x(tr, "ReconstructedChargedParticles.p.x");
  TTreeReaderArray<float> ReconstructedChargedParticles_p_y(tr, "ReconstructedChargedParticles.p.y");
  TTreeReaderArray<float> ReconstructedChargedParticles_p_z(tr, "ReconstructedChargedParticles.p.z");
  TTreeReaderArray<float> ReconstructedChargedParticles_mass(tr, "ReconstructedChargedParticles.mass");
  TTreeReaderArray<Int_t> ReconstructedChargedParticles_mcID_value(tr, "ReconstructedChargedParticles.mcID.value");

  // Inclusive Kinematics Truth
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_x(tr, "InclusiveKinematicsTruth.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_y(tr, "InclusiveKinematicsTruth.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_Q2(tr, "InclusiveKinematicsTruth.Q2");
  TTreeReaderArray<Int_t> InclusiveKinematicsTruth_scatID_value(tr, "InclusiveKinematicsTruth.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsTruth_scatID_source(tr, "InclusiveKinematicsTruth.scatID.source");
  
  // Inclusive Kinematics Electron 
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_x(tr, "InclusiveKinematicsElectron.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_y(tr, "InclusiveKinematicsElectron.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_Q2(tr, "InclusiveKinematicsElectron.Q2");
  TTreeReaderArray<Int_t> InclusiveKinematicsElectron_scatID_value(tr, "InclusiveKinematicsElectron.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsElectron_scatID_source(tr, "InclusiveKinematicsElectron.scatID.source");

  // Inclusive Kinematics JB
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_x(tr, "InclusiveKinematicsJB.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_y(tr, "InclusiveKinematicsJB.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_Q2(tr, "InclusiveKinematicsJB.Q2");
  TTreeReaderArray<Int_t> InclusiveKinematicsJB_scatID_value(tr, "InclusiveKinematicsJB.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsJB_scatID_source(tr, "InclusiveKinematicsJB.scatID.source");

  // Inclusive Kinematics DA
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_x(tr, "InclusiveKinematicsDA.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_y(tr, "InclusiveKinematicsDA.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_Q2(tr, "InclusiveKinematicsDA.Q2");
  TTreeReaderArray<Int_t> InclusiveKinematicsDA_scatID_value(tr, "InclusiveKinematicsDA.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsDA_scatID_source(tr, "InclusiveKinematicsDA.scatID.source");

  // //HcalEndcap
  // TTreeReaderArray<float>  HcalEndcapPClusters_energy(tr, "HcalEndcapPClusters.energy");
  // TTreeReaderArray<float> HcalEndcapPClusters_x(tr,      "HcalEndcapPClusters.position.x");
  // TTreeReaderArray<float> HcalEndcapPClusters_y(tr,      "HcalEndcapPClusters.position.y");
  // TTreeReaderArray<float> HcalEndcapPClusters_z(tr,      "HcalEndcapPClusters.position.z");
  // TTreeReaderArray<float>  HcalEndcapPClusters_theta(tr,  "HcalEndcapPClusters.polar.theta");
  // TTreeReaderArray<float>  HcalEndcapPClusters_phi(tr,    "HcalEndcapPClusters.polar.phi");
  // 
  // //HcalEndcap
  // TTreeReaderArray<float> HcalEndcapNClusters_energy(tr, "HcalEndcapNClusters.energy");
  // TTreeReaderArray<float> HcalEndcapNClusters_x(tr,      "HcalEndcapNClusters.position.x");
  // TTreeReaderArray<float> HcalEndcapNClusters_y(tr,      "HcalEndcapNClusters.position.y");
  // TTreeReaderArray<float> HcalEndcapNClusters_z(tr,      "HcalEndcapNClusters.position.z");
  // TTreeReaderArray<float> HcalEndcapNClusters_theta(tr,  "HcalEndcapNClusters.polar.theta");
  // TTreeReaderArray<float> HcalEndcapNClusters_phi(tr,    "HcalEndcapNClusters.polar.phi");
  // 
  // //HcalBarrel
  // TTreeReaderArray<float> HcalBarrelClusters_energy(tr, "HcalBarrelClusters.energy");
  // TTreeReaderArray<float> HcalBarrelClusters_x(tr,      "HcalBarrelClusters.position.x");
  // TTreeReaderArray<float> HcalBarrelClusters_y(tr,      "HcalBarrelClusters.position.y");
  // TTreeReaderArray<float> HcalBarrelClusters_z(tr,      "HcalBarrelClusters.position.z");
  // TTreeReaderArray<float> HcalBarrelClusters_theta(tr,  "HcalBarrelClusters.polar.theta");
  // TTreeReaderArray<float> HcalBarrelClusters_phi(tr,    "HcalBarrelClusters.polar.phi");
  // 
  // //Ecal
  // TTreeReaderArray<float> EcalEndcapPClusters_energy(tr, "EcalEndcapPClusters.energy");
  // TTreeReaderArray<float> EcalEndcapPClusters_x(tr,      "EcalEndcapPClusters.position.x");
  // TTreeReaderArray<float> EcalEndcapPClusters_y(tr,      "EcalEndcapPClusters.position.y");
  // TTreeReaderArray<float> EcalEndcapPClusters_z(tr,      "EcalEndcapPClusters.position.z");
  // TTreeReaderArray<float> EcalEndcapPClusters_theta(tr,  "EcalEndcapPClusters.polar.theta");
  // TTreeReaderArray<float> EcalEndcapPClusters_phi(tr,    "EcalEndcapPClusters.polar.phi");
  // 
  // TTreeReaderArray<float> EcalEndcapNClusters_energy(tr, "EcalEndcapNClusters.energy");
  // TTreeReaderArray<float> EcalEndcapNClusters_x(tr,      "EcalEndcapNClusters.position.x");
  // TTreeReaderArray<float> EcalEndcapNClusters_y(tr,      "EcalEndcapNClusters.position.y");
  // TTreeReaderArray<float> EcalEndcapNClusters_z(tr,      "EcalEndcapNClusters.position.z");
  // //watch out the branch name
  // TTreeReaderArray<float>  EcalEndcapNClusters_theta(tr,  "EcalEndcapNClusters.polar.theta");
  // TTreeReaderArray<float>  EcalEndcapNClusters_phi(tr,    "EcalEndcapNClusters.polar.phi");
  // 
  // TTreeReaderArray<float> EcalBarrelClusters_energy(tr, "EcalBarrelImagingClusters.energy");
  // TTreeReaderArray<float> EcalBarrelClusters_x(tr,      "EcalBarrelImagingClusters.position.x");
  // TTreeReaderArray<float> EcalBarrelClusters_y(tr,      "EcalBarrelImagingClusters.position.y");
  // TTreeReaderArray<float> EcalBarrelClusters_z(tr,      "EcalBarrelImagingClusters.position.z");
  // TTreeReaderArray<float> EcalBarrelClusters_theta(tr,  "EcalBarrelImagingClusters.polar.theta");
  // TTreeReaderArray<float>  EcalBarrelClusters_phi(tr,    "EcalBarrelImagingClusters.polar.phi");


  // Setup tree to store kinematic variables
  Double_t x_true, x_e_ecal, x_e_track, x_jb, x_da, x_sig, x_esig;
  Double_t y_true, y_e_ecal, y_e_track, y_jb, y_da, y_sig, y_esig;
  Double_t Q2_true, Q2_e_ecal, Q2_e_track, Q2_jb, Q2_da, Q2_sig, Q2_esig;


  TFile *varTreeFile = new TFile("TreeOutput.root", "RECREATE");
  TTree *varTree = new TTree("event_var", "event_var");
  varTree->Branch("x_true", &x_true, "x_true/D");
  //varTree->Branch("x_e_ecal", &x_e_nm_s, "x_e_nm_s/D");
  varTree->Branch("x_e_track", &x_e_track, "x_e_track/D");
  varTree->Branch("x_jb", &x_jb, "x_jb/D");
  varTree->Branch("x_da", &x_da, "x_da/D");
  varTree->Branch("x_sig", &x_sig, "x_sig/D");
  varTree->Branch("x_esig", &x_esig, "x_esig/D");

  varTree->Branch("y_true", &y_true, "y_true/D");
  //varTree->Branch("y_e_ecal", &y_e_nm_s, "y_e_nm_s/D");
  varTree->Branch("y_e_track", &y_e_track, "y_e_track/D");
  varTree->Branch("y_jb", &y_jb, "y_jb/D");
  varTree->Branch("y_da", &y_da, "y_da/D");
  varTree->Branch("y_sig", &y_sig, "y_sig/D");
  varTree->Branch("y_esig", &y_esig, "y_esig/D");

  varTree->Branch("Q2_true", &Q2_true, "Q2_true/D");
  //varTree->Branch("Q2_e_ecal", &Q2_e_nm_s, "Q2_e_nm_s/D");
  varTree->Branch("Q2_e_track", &Q2_e_track, "Q2_e_track/D");
  varTree->Branch("Q2_jb", &Q2_jb, "Q2_jb/D");
  varTree->Branch("Q2_da", &Q2_da, "Q2_da/D");
  varTree->Branch("Q2_sig", &Q2_sig, "Q2_sig/D");
  varTree->Branch("Q2_esig", &Q2_esig, "Q2_esig/D");

  int true_index (0), reco_index (0);
  double E_scattered (0.0);


  //TFile *testTreeFile = new TFile("TreeOutput.root", "RECREATE");
  //TTree *varTree = new TTree("event", "event");
  //varTree->Branch("true_index", &true_index, "true_index/I");
  //varTree->Branch("reco_index", &reco_index, "reco_index/I");
  //varTree->Branch("E_scattered", &E_scattered, "E_scattered/D");
  
  TTreeReader::EEntryStatus entrystats = tr.SetEntry(0);
  
  TLorentzVector ei;
  TLorentzVector ef;
  TLorentzVector pni;
  TLorentzVector q_e;
  
  TLorentzVector vecElectron;
  TLorentzVector vecHadron;
  
  
  Double_t m_e(0.000511);
  Double_t Q2_e;
  Double_t x_e;
  Double_t y_e;
  Double_t y_h;
  
  Double_t sigmah;
  
  TH1D *Q2Hist = new TH1D("Q2Hist", "Q2Hist", 100, 0, 200);
  TH1D *xHist = new TH1D("xHist", "xHist", 100, 0, 1.2);
  TH1D *yHist = new TH1D("yHist", "yHist", 100, 0, 1.2);
  TH1D *yh_noBoost = new TH1D("yh_noBoost", "yh_noBoost", 1000, -2, 2);
  TH1D *yh_boosted = new TH1D("yh_boosted", "yh_boosted", 1000, -2, 2);
  TH2D *eta_vs_Q2 = new TH2D("eta_vs_Q2", ";Q^{2};#eta_{e}", 1000, 0, 300, 1000, -4, 4);
  
  int counter {0};
  while (tr.Next()){

    if (!InclusiveKinematicsTruth_x.GetSize()) continue;

    counter++;
    //if (counter > 10) break;
    //cout << typeid((float)InclusiveKinematicsTruth_x.At(0)).name() << endl;
    //cout << InclusiveKinematicsTruth_x.At(0) << endl;

    x_true = InclusiveKinematicsTruth_x[0];
    y_true = InclusiveKinematicsTruth_y[0];
    Q2_true = InclusiveKinematicsTruth_Q2[0];

    if (!InclusiveKinematicsElectron_x.GetSize()) continue;
    x_e_track = InclusiveKinematicsElectron_x[0];
    y_e_track = InclusiveKinematicsElectron_y[0];
    Q2_e_track = InclusiveKinematicsElectron_Q2[0];

    if (!InclusiveKinematicsJB_x.GetSize()) continue;
    x_jb= InclusiveKinematicsJB_x[0];
    y_jb = InclusiveKinematicsJB_y[0];
    Q2_jb = InclusiveKinematicsJB_Q2[0];

    if (!InclusiveKinematicsDA_x.GetSize()) continue;
    x_da= InclusiveKinematicsDA_x[0];
    y_da = InclusiveKinematicsDA_y[0];
    Q2_da = InclusiveKinematicsDA_Q2[0];

    //cout << InclusiveKinematicsTruth_x[0] << " " << InclusiveKinematicsTruth_y[0]  << " " << InclusiveKinematicsTruth_Q2[0] << endl;

    Double_t mc_neutral_E(0),mc_neutral_px(0), mc_neutral_py(0), mc_neutral_pz(0);
    Double_t max_electron_energy = 0;
    for (int i=0; i<mcparticles_pdgID.GetSize(); i++){
      int pid = mcparticles_pdgID[i];
      double psx = mcparticles_psx[i];
      double psy = mcparticles_psy[i];
      double psz = mcparticles_psz[i];
      double pT = TMath::Sqrt(psx*psx + psy*psy);
      double mass = mcparticles_mass[i];
      int charge = mcparticles_charge[i];
      double energy = TMath::Sqrt(psx*psx + psy*psy + psz*psz + mass*mass);
      double smeared_energy(0);
      double eta = TMath::ATanH(psz/TMath::Sqrt(energy*energy - mass*mass));

      double stochastic_a(0), constant_c(0), sigma_energy(0);


      if (mcparticles_genStatus[i] == 4 && pid == 11){
	ei.SetPxPyPzE(psx, psy, psz, energy);
      }
      // Find beam p
      if (mcparticles_genStatus[i] == 4 && pid == 2212){
	pni.SetPxPyPzE(psx, psy, psz, energy);
      }
      if (i == InclusiveKinematicsTruth_scatID_value[0]){
	ef.SetPxPyPzE(psx, psy, psz, energy);
	true_index = i;
      }
      if (mcparticles_genStatus[i] == 1 && pid == 11){
	if (energy > max_electron_energy){
	  max_electron_energy = energy;
	  reco_index = i;
	}
      }
      if (mcparticles_genStatus[i] == 1 && charge == 0 && TMath::Abs(eta) < 4){
	if (eta > -1 && eta < 1){
	  constant_c = 0.37;
	}
	else if (eta > 1 && eta < 4){
	  stochastic_a = 0.4;
	  constant_c = 0.06;
	}
	sigma_energy = TMath::Sqrt(stochastic_a*stochastic_a*energy + constant_c*constant_c*energy*energy);
	smeared_energy = rndm->Gaus(energy, sigma_energy);
	psx = psx*smeared_energy/energy;
	psy = psy*smeared_energy/energy;
	psz = psz*smeared_energy/energy;
	
	mc_neutral_E+=smeared_energy;
	mc_neutral_px+=psx;
	mc_neutral_py+=psy;
	mc_neutral_pz+=psz;
      }
    }
    if (max_electron_energy == 0) continue;

    
    E_scattered = ef.E();
    

    eta_vs_Q2->Fill(InclusiveKinematicsTruth_Q2[0], ef.Eta());
    //Q2Hist->Fill(InclusiveKinematicsTruth_Q2[0]);
    //xHist->Fill(InclusiveKinematicsTruth_x[0]);
    //yHist->Fill(InclusiveKinematicsTruth_y[0]);

    // track loop
    // int electron_index(-1);
    // double hpx=0;
    // double hpy=0;
    // double hpz=0;
    // double hE=0;
    // double maxP=0;
    // for(int itrk=0; itrk<(int)ReconstructedChargedParticles_pid.GetSize(); itrk++)
    // {
    // // FIXME: pid is using the true information
    // // Add PID smearing                       
    // int pid_ = ReconstructedChargedParticles_pid[itrk];
    // 
    // // pid==0: reconstructed tracks with no matching truth pid
    // if(pid_ == 0) continue;
    // 
    // double reco_E = ReconstructedChargedParticles_energy[itrk];
    // double reco_px = ReconstructedChargedParticles_p_x[itrk];
    // double reco_py = ReconstructedChargedParticles_p_y[itrk];
    // double reco_pz = ReconstructedChargedParticles_p_z[itrk];
    // double reco_mass = ReconstructedChargedParticles_mass[itrk];
    // double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
    // hpx += reco_px;
    // hpy += reco_py;
    // hpz += reco_pz;
    // hE += reco_E;
    // 
    // if (pid_ ==  11 && ReconstructedChargedParticles_mcID_value[itrk] == true_index){
    // vecElectron.SetPxPyPzE(reco_px, reco_py, reco_pz, reco_E);
    // electron_index = itrk;
    // }
    // }
    // if (electron_index == -1) continue;

    int electron_index(-1);
    double hpx=0;
    double hpy=0;
    double hpz=0;
    double hE=0;
    double maxP=0;
    for(int itrk=0; itrk<(int)ReconstructedParticles_pid.GetSize(); itrk++)
      {
	// FIXME: pid is using the true information
	// Add PID smearing                       
	int pid_ = ReconstructedParticles_pid[itrk];
	
	// pid==0: reconstructed tracks with no matching truth pid
	// if(pid_ == 0) continue;
	// if(ReconstructedParticles_charge[itrk] == 0) continue;
	
	double reco_E = ReconstructedParticles_energy[itrk];
	double reco_px = ReconstructedParticles_p_x[itrk];
	double reco_py = ReconstructedParticles_p_y[itrk];
	double reco_pz = ReconstructedParticles_p_z[itrk];
	double reco_mass = ReconstructedParticles_mass[itrk];
	double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
	hpx += reco_px;
	hpy += reco_py;
	hpz += reco_pz;
	hE += reco_E;
	
	if (pid_ ==  11 && ReconstructedParticles_mcID_value[itrk] == true_index){
	  vecElectron.SetPxPyPzE(reco_px, reco_py, reco_pz, reco_E);
	  electron_index = itrk;
	}
      }
    if (electron_index == -1) continue;
    if (ReconstructedParticles_pid.GetSize() < 2) continue;
      
    //vecHadron.SetPxPyPzE(hpx+mc_neutral_px, hpy+mc_neutral_py, hpz+mc_neutral_pz, hE+mc_neutral_E);
    vecHadron.SetPxPyPzE(hpx, hpy, hpz, hE);
    vecHadron -= vecElectron;

    auto BoostedVecElectron = TransformLabToHeadOnFrame(ei, pni, vecElectron);
    auto BoostedVecHadron = TransformLabToHeadOnFrame(ei, pni, vecHadron);

    sigmah = BoostedVecHadron.E() - BoostedVecHadron.Pz();

    //Q2_e_track = 4.*ei.E()*BoostedVecElectron.E()*TMath::Cos(BoostedVecElectron.Theta()/2)*TMath::Cos(BoostedVecElectron.Theta()/2);
    //y_e_track = 1. - ((0.5*BoostedVecElectron.E()/ei.E())*(1 - TMath::Cos(BoostedVecElectron.Theta())));
    //x_e_track = Q2_e_track/(4*ei.E()*pni.E()*y_e_track);
    //q_e = ei - vecElectron;
    //Q2_e_track = -1.0*(q_e*q_e);
    //y_e_track = (pni*q_e)/(pni*ei);
    //x_e_track = Q2_e_track/(2*pni*q_e); 


    // y_jb = sigmah/(2*ei.E());
    // Q2_jb = (BoostedVecHadron.Pt()*BoostedVecHadron.Pt())/(1. - y_jb);
    // x_jb = Q2_jb/(4*ei.E()*pni.E()*y_jb);

    y_sig = sigmah/(sigmah + BoostedVecElectron.E()*(1 - TMath::Cos(BoostedVecElectron.Theta())));
    Q2_sig = BoostedVecElectron.Pt()*BoostedVecElectron.Pt()/(1 - y_sig);
    x_sig = Q2_sig/(4*ei.E()*pni.E()*y_sig);

    Q2_esig = Q2_e_track;
    x_esig = x_sig;
    y_esig = Q2_esig/(4*ei.E()*pni.E()*x_esig);

    //Q2_da = 4.*ei.E()*ei.E()*( 1./TMath::Tan(BoostedVecElectron.Theta()/2.) )*
    //( 1./(TMath::Tan(BoostedVecElectron.Theta()/2.)+TMath::Tan(BoostedVecHadron.Theta()/2.)) );
    //y_da = TMath::Tan(BoostedVecHadron.Theta()/2)/(TMath::Tan(BoostedVecElectron.Theta()/2)+TMath::Tan(BoostedVecHadron.Theta()/2));
    //x_da = Q2_da/(4*ei.E()*pni.E()*y_da);


    //cout << x_true << " " << y_true << " " << Q2_true << endl;

    //  
    //  
    //  
    // q_e = ei - ef;
    // Q2_e = -1.0*(q_e*q_e);
    // Q2Hist->Fill(Q2_e);
    // y_e = (pni*q_e)/(pni*ei);
    // x_e = Q2_e/(2*pni*q_e); 
    // yHist->Fill(y_e);
    // xHist->Fill(x_e);
    //  
    //  
    //  
    //vecHadron.SetPxPyPzE(hpx,hpy,hpz,hE);
    //vecHadron -= vecElectron;
    // 
    //sigmah = vecHadron.E() - vecHadron.Pz();
    //yh_noBoost->Fill((sigmah/36));
    // 
    //vecHadron = TransformLabToHeadOnFrame(ei, pni, vecHadron);
    //sigmah = vecHadron.E() - vecHadron.Pz();
    //yh_boosted->Fill((sigmah/36));
    
    
    
    
    
    
    varTree->Fill();
  }//end event loop 

  varTreeFile->Write();
  eta_vs_Q2->Draw("colz");
}

TLorentzVector TransformLabToHeadOnFrame(TLorentzVector eBeam, TLorentzVector pBeam, TLorentzVector Lvec) {
  TLorentzVector cmBoost = eBeam + pBeam;
  // Define boost for back to back beams
  TLorentzVector boost(-cmBoost[0],-cmBoost[1],-cmBoost[2],cmBoost[3]);
  TVector3 b = boost.BoostVector();
  
  // Define boost to restore original energies
  TLorentzVector boostBack(0.0,0.0,cmBoost[2],cmBoost[3]);
  TVector3 bb = boostBack.BoostVector();
  
  // Define angles for rotations
  TLorentzVector p = pBeam;
  p.Boost(b);
  double rotAboutY = -1.0*TMath::ATan2(p.Px(),p.Pz());
  double rotAboutX = -1.0*TMath::ATan2(p.Py(),p.Pz());
  
  Lvec.Boost(b);
  Lvec.RotateY(rotAboutY);
  Lvec.RotateX(rotAboutX);
  Lvec.Boost(bb);
  return Lvec;
}
