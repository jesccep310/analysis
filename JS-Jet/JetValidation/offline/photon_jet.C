void photon_jet() 
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  //get our files
  TChain * ct = new TChain("T");  
  for(int i = 0; i < 15000; i++){
    ct->Add(Form("/sphenix/user/vbailey/MDC2pythiajets_photonjet/output%i_r04.root",i));
  }


  //read in the trees
  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *pt = 0;
  vector<float> *e = 0;
  vector<float> *subtracted_et = 0;
  vector<float> *truthEta = 0;
  vector<float> *truthPhi = 0;
  vector<float> *truthPt = 0;
  vector<float> *truthE = 0;
  vector<float> *clusterEta = 0;
  vector<float> *clusterPhi = 0;
  vector<float> *clusterPt = 0;
  vector<float> *clusterE = 0;
  vector<float> *photonProb = 0;
  vector<float> *photonChi = 0;
  vector<float> *photonIso = 0;
  int cent;

  ct->SetBranchAddress("eta",&eta);
  ct->SetBranchAddress("phi",&phi);
  ct->SetBranchAddress("pt",&pt);
  ct->SetBranchAddress("e",&e);
  ct->SetBranchAddress("subtracted_et",&subtracted_et);
  ct->SetBranchAddress("truthEta",&truthEta);
  ct->SetBranchAddress("truthPhi",&truthPhi);
  ct->SetBranchAddress("truthPt",&truthPt);
  ct->SetBranchAddress("truthE",&truthE);
  ct->SetBranchAddress("cent", &cent);

  ct->SetBranchAddress("clusterEta",&clusterEta);
  ct->SetBranchAddress("clusterPhi",&clusterPhi);
  ct->SetBranchAddress("clusterPt",&clusterPt);
  ct->SetBranchAddress("clusterE",&clusterE);
  ct->SetBranchAddress("photonProb",&photonProb);
  ct->SetBranchAddress("photonChi",&photonChi);
  ct->SetBranchAddress("photonIso",&photonIso);

  //set up binning for the histograms
  const Float_t pt_range[] = {10,15,20,25,30,35,40,45,50,60,80};
  const int pt_N = sizeof(pt_range)/sizeof(float) - 1;
  const int resp_N = 500;
  Float_t resp_bins[resp_N];
  for(int i = 0; i < resp_N+1; i++)
    {
      resp_bins[i] = 5.0/resp_N * i;
    }
  const int eta_N = 40;
  Float_t eta_bins[eta_N];
  for(int i = 0; i < eta_N+1; i++)
    {
      eta_bins[i] = -1.1 + 2.2/eta_N * i;
    }
  const int phi_N = 64;
  Float_t phi_bins[phi_N];
  for(int i = 0; i < phi_N+1; i++)
    {
      phi_bins[i] = TMath::Pi()/phi_N * i;
    }
    
  //create some histograms
  TH1F *h_jet_pt = new TH1F("h_jet_pt","",pt_N, pt_range);
  TH1F *h_cluster_pt = new TH1F("h_cluster_pt","",pt_N, pt_range);
  TH2F *h_response = new TH2F("h_response","",pt_N, pt_range, resp_N, resp_bins);
  TH2F *h_dPhi = new TH2F("h_dPhi","",pt_N, pt_range, phi_N, phi_bins);

  //loop through each event and fill the histograms
  int nentries = ct->GetEntries();
  for(int i = 0; i < nentries; i++){
    ct->GetEntry(i);
    
    int nrecojets = pt->size();  
    int nclusters = clusterPt->size();
    
    //get the leading cluster
    int maxCluster = -1;
    float maxClusterPt= 0;
    for(int cl = 0; cl < nclusters; cl++)
      {
	if(photonChi->at(cl) > 3) continue;
	if(clusterPt->at(cl) > maxClusterPt)
          {
            maxCluster = cl;
            maxClusterPt = clusterPt->at(cl);
          }
      }
    if(maxCluster < 0) continue;
    h_cluster_pt->Fill(clusterPt->at(maxCluster)); //fill the pt of the leading cluster
  
    int maxJet = -1;
    float maxJetPt = 0;
  
    int njets = truthPt->size();
    for(int tj = 0; tj < njets; tj++){

      int nrecojets = pt->size();
      float dR;

    
      //we already found the max cluster!
      /*float maxClusterPt = 900; 
      int nclusters = clusterPt->size();
      for(int cl = 0; cl < nclusters; cl++)
	{
	  if(photonChi->at(cl) > 3) continue;
	  if(clusterPt->at(cl) > maxClusterPt)
	    {
	      maxCluster = cl;
	      maxClusterPt = clusterPt->at(cl);
	    }
	  if(clusterPt->at(cl) < 0.5*truthPt->at(tj)) continue;
	  }*/

      //instead do something like this
      float closestClusterdR = 999;
      int closestCluster = -1;
      //find the closest cluster to the truth jet
      for(int cl = 0; cl < nclusters; cl++)
	{
	  if(photonChi->at(cl) > 3) continue;
	  float dEta = truthEta->at(tj) - clusterEta->at(cl);
	  float dPhi = truthPhi->at(tj) - clusterPhi->at(cl);
	  while(dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
	  while(dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();
	  dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	  if(dR < closestClusterdR)
	    {
	      closestClusterdR = dR;
	      closestCluster = cl;
	    }
	}
      //skip this truth jet if its got a cluster with >0.5*its pt within dR of 0.2
      if(closestClusterdR < 0.2 && clusterPt->at(closestCluster) > 0.5*truthPt->at(tj)) continue;

 
      //reco to truth jet matching
      float matchEta, matchPhi, matchPt, matchE, matchsubtracted_et;
      float dRMax = 100;
      int matchJet;
      
      for(int rj = 0; rj < nrecojets; rj++){
	float dEta = truthEta->at(tj) - eta->at(rj);
	float dPhi = truthPhi->at(tj) - phi->at(rj);
	while(dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
	while(dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();
	dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	if(dR < dRMax){
	  matchEta = eta->at(rj);
	  matchPhi = phi->at(rj);
	  matchPt = pt->at(rj);
	  matchE = e->at(rj);
	  matchsubtracted_et = subtracted_et->at(rj);
	  dRMax = dR;
	  matchJet = rj;
	}
      }
      //now we have a reco jet which we think is not a photon
      //check if its the highest pt jet in the event
      if(matchPt > maxJetPt)
	{
	  maxJetPt = matchPt;
	  maxJet = matchJet;
	}
    }
  
    if(maxCluster < 0 || maxJet < 0) continue;
   
    float dEta = clusterEta->at(maxCluster) - eta->at(maxJet);
    float dPhi = clusterPhi->at(maxCluster) - phi->at(maxJet);
    
    //get dphi between -pi and pi
    while (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
    while (dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();

    float dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

    h_dPhi->Fill(clusterPt->at(maxCluster), abs(dPhi));

    if(abs(dPhi) > 2.9)
      {
	float response = pt->at(maxJet)/clusterPt->at(maxCluster);
	h_response->Fill(clusterPt->at(maxCluster), response);
      }

  }






  TCanvas *c = new TCanvas("c","c");
  c->Print("h1D1plotrecotruthjet.pdf(");

  for (int i = 0; i < pt_N; i++) {
    h_response->GetXaxis()->SetRange(i+1,i+1);
    TH1F *h1D = (TH1F*) h_response->ProjectionY();
    h1D->Draw();
    c->Print("h1D1plotrecotruthjet.pdf");
  }
  c->Print("h1D1plotrecotruthjet.pdf)");

  //write out histograms to the output file
  TFile *f_out = new TFile("photon_jet_histsrecotruthjetdPhi29.root","RECREATE");
  h_jet_pt->Write();
  h_cluster_pt->Write();
  h_response->Write();
  h_dPhi->Write();
}
