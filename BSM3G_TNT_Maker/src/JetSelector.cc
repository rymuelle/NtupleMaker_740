// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"

JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the JetSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
//  jetToken_       = iConfig.getParameter<edm::InputTag>("jets");
  jetToken_       = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
//  puppi_jetToken_ = iConfig.getParameter<edm::InputTag>("jetsPUPPI");
  puppi_jetToken_ = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsPUPPI"));
  AK8_puppi_jetToken_ = iCC.consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAk8PUPPI"));
//  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  _vertexInputTag = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

JetSelector::~JetSelector(){
  delete tree_;
}

void JetSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  // grab handle to the jet collections
  edm::Handle<pat::JetCollection> jets;                                       
//  iEvent.getByLabel(jetToken_, jets);                                         
  iEvent.getByToken(jetToken_, jets);                                         
  edm::Handle<pat::JetCollection> puppijets;                                       
//  iEvent.getByLabel(puppi_jetToken_, puppijets);                                         
  iEvent.getByToken(puppi_jetToken_, puppijets);    

  edm::Handle<pat::JetCollection> ak8puppijets;                                       
//  iEvent.getByLabel(puppi_jetToken_, puppijets);                                         
  iEvent.getByToken(AK8_puppi_jetToken_, ak8puppijets);                                        

  if(debug_) std::cout << "     JetSelector: Cleared the vectors, grabbed the jet collection handle, and looping over jets." << std::endl;
  
  // loop over "standard" ak4 CHS PF jets
  for (const pat::Jet &j : *jets) { 
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_pt.push_back(j.pt());         
    Jet_eta.push_back(j.eta());       
    Jet_phi.push_back(j.phi());       
    Jet_energy.push_back(j.energy());
    Jet_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_bDiscriminator_CISVV2.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet_vtxMass.push_back(j.userFloat("vtxMass"));
    //Jet_decayLength.push_back(j.userFloat("vtx3DVal"));
    //Jet_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    //Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_partonFlavour.push_back(j.partonFlavour());
    Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_mass.push_back(j.mass()); 
      Jet_electronEnergy.push_back(j.electronEnergy());                               
      Jet_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }
  } 
  
  // loop over ak4 PF jets reconstructed with the PUPPI pileup mitigation algorithms
  for (const pat::Jet &j : *puppijets) { 
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_puppi_pt.push_back(j.pt());         
    Jet_puppi_eta.push_back(j.eta());       
    Jet_puppi_phi.push_back(j.phi());       
    Jet_puppi_energy.push_back(j.energy());
    Jet_puppi_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_puppi_bDiscriminator_CISVV2.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_puppi_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet_puppi_vtxMass.push_back(j.userFloat("vtxMass"));
    //Jet_puppi_decayLength.push_back(j.userFloat("vtx3DVal"));
    //Jet_puppi_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    //Jet_puppi_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_puppi_partonFlavour.push_back(j.partonFlavour());
    Jet_puppi_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_puppi_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_puppi_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_puppi_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_puppi_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_puppi_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_puppi_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_puppi_mass.push_back(j.mass()); 
      Jet_puppi_electronEnergy.push_back(j.electronEnergy());                               
      Jet_puppi_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_puppi_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }
  } 

  int jet_counter = 0;
  for (const pat::Jet &j : *ak8puppijets) { 
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_AK8_pt.push_back(j.pt());         
    Jet_AK8_eta.push_back(j.eta());       
    Jet_AK8_phi.push_back(j.phi());       
    Jet_AK8_energy.push_back(j.energy());
    Jet_AK8_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_AK8_bDiscriminator_CISVV2.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_AK8_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    //Jet_AK8_vtxMass.push_back(j.userFloat("vtxMass"));
    //Jet_AK8_decayLength.push_back(j.userFloat("vtx3DVal"));
    //Jet_AK8_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    //Jet_AK8_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_AK8_partonFlavour.push_back(j.partonFlavour());
    Jet_AK8_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_AK8_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_AK8_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_AK8_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_AK8_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_AK8_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_AK8_chargedMultiplicity.push_back(j.chargedMultiplicity());
    /*if(j.subjets(0).size() > 0 ){
      Jet_AK8_vSubjetIndex0.push_back(j.subjets(0).at(0).key());
    } else {
      Jet_AK8_vSubjetIndex0.push_back(-1);
    }*/

    Jet_AK8_puppi_pt.push_back(j.userFloat("ak8PFJetsPuppiValueMap:pt"));
    Jet_AK8_puppi_eta.push_back(j.userFloat("ak8PFJetsPuppiValueMap:eta"));
    Jet_AK8_puppi_mass.push_back(j.userFloat("ak8PFJetsPuppiValueMap:mass"));
    Jet_AK8_puppi_phi.push_back(j.userFloat("ak8PFJetsPuppiValueMap:phi"));
    Jet_AK8_puppi_tau1.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));
    Jet_AK8_puppi_tau2.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
    Jet_AK8_puppi_tau3.push_back(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3"));

    //vector <double> Jet_AK8_puppi_pt, Jet_AK8_puppi_energy,Jet_AK8_puppi_phi, Jet_AK8_puppi_eta, Jet_AK8_puppi_mass, Jet_AK8_puppi_tau2, Jet_AK8_puppi_tau3, Jet_AK8_puppi_;
    //vector <double>  Jet_AK8_subjet0_pt, Jet_AK8_subjet1_pt, Jet_AK8_subjet0_eta, Jet_AK8_subjet1_eta, Jet_AK8_subjet0_phi, Jet_AK8_subjet1_phi, Jet_AK8_subjet0_energy, Jet_AK8_subjet1_energy,  Jet_AK8_subjet0_CSVv2, Jet_AK8_subjet1_CSVv2 ;
   

    auto const & sdSubjetsPuppi = j.subjets("SoftDropPuppi");


    if(sdSubjetsPuppi.size() > 0){
      Jet_AK8_subjet0_pt.push_back(sdSubjetsPuppi[0]->correctedP4(0).pt());
      Jet_AK8_subjet0_phi.push_back(sdSubjetsPuppi[0]->correctedP4(0).phi());
      Jet_AK8_subjet0_eta.push_back(sdSubjetsPuppi[0]->correctedP4(0).eta());
      Jet_AK8_subjet0_mass.push_back(sdSubjetsPuppi[0]->correctedP4(0).mass());
      Jet_AK8_subjet0_energy.push_back(sdSubjetsPuppi[0]->correctedP4(0).energy());
      Jet_AK8_subjet0_CSVv2.push_back(sdSubjetsPuppi[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    } else {
      Jet_AK8_subjet0_pt.push_back(-100);
      Jet_AK8_subjet0_phi.push_back(-100);
      Jet_AK8_subjet0_eta.push_back(-100);
      Jet_AK8_subjet0_mass.push_back(-100);
      Jet_AK8_subjet0_energy.push_back(-100);
      Jet_AK8_subjet0_CSVv2.push_back(-100);
    }

    if(sdSubjetsPuppi.size() > 1){
      Jet_AK8_subjet1_pt.push_back(sdSubjetsPuppi[1]->correctedP4(0).pt());
      Jet_AK8_subjet1_phi.push_back(sdSubjetsPuppi[1]->correctedP4(0).phi());
      Jet_AK8_subjet1_eta.push_back(sdSubjetsPuppi[1]->correctedP4(0).eta());
      Jet_AK8_subjet1_mass.push_back(sdSubjetsPuppi[1]->correctedP4(0).mass());
      Jet_AK8_subjet1_energy.push_back(sdSubjetsPuppi[1]->correctedP4(0).energy());
      Jet_AK8_subjet1_CSVv2.push_back(sdSubjetsPuppi[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    } else {
      Jet_AK8_subjet1_pt.push_back(-100);
      Jet_AK8_subjet1_phi.push_back(-100);
      Jet_AK8_subjet1_eta.push_back(-100);
      Jet_AK8_subjet1_mass.push_back(-100);
      Jet_AK8_subjet1_energy.push_back(-100);
      Jet_AK8_subjet1_CSVv2.push_back(-100);
    }



    if(j.genJet()){
      Jet_AK8_GEN_pt.push_back(j.genJet()->pt());
      Jet_AK8_GEN_phi.push_back(j.genJet()->phi());
      Jet_AK8_GEN_eta.push_back(j.genJet()->eta());
      Jet_AK8_GEN_mass.push_back(j.genJet()->mass());
      Jet_AK8_GEN_energy.push_back(j.genJet()->energy());
      //Jet_AK8_GEN_parton.push_back(genJet().pdgId());
     // Jet_AK8_GEN_mother.push_back(genJet().mother().pdgId());
      Jet_AK8_GEN_parton.push_back(-100);
      Jet_AK8_GEN_mother.push_back(-100);

    } else {
      Jet_AK8_GEN_pt.push_back(-100);
      Jet_AK8_GEN_phi.push_back(-100);
      Jet_AK8_GEN_eta.push_back(-100);
      Jet_AK8_GEN_mass.push_back(-100);
      Jet_AK8_GEN_parton.push_back(-100);
      Jet_AK8_GEN_mother.push_back(-100);
    }


    /*vector<double> subjet_pt;
    Jet_AK8_subjet_pt.push_back(subjet_pt);

    auto const & sdSubjetsPuppi = j.subjets("SoftDropPuppi");
    for ( auto const & it : sdSubjetsPuppi ) {
      Jet_AK8_subjet_pt[jet_counter].push_back(it->correctedP4(0).pt());
      std::cout << "subjet stuff " << it->correctedP4(0).pt() << " " << j.pt() << std::endl;
    }*/

    //for ( auto const & it : sbSubjetsPuppi ) {
    //for ( unsigned int it =0; it  < sdSubjetsPuppi.size(); it++ ) {
      //Jet_AK8_subjet_pt.push_back(it->correctedP4(0).pt());
    //  Jet_AK8_subjet_pt.push_back(sdSubjetsPuppi.at(it).correctedP4(0).pt());
    //}

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_AK8_mass.push_back(j.mass()); 
      Jet_AK8_electronEnergy.push_back(j.electronEnergy());                               
      Jet_AK8_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_AK8_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }
  }
  jet_counter++;

}

void JetSelector::SetBranches(){
  if(debug_) std::cout << "     JetSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Jet_pt,                  		"Jet_pt");
  AddBranch(&Jet_eta,                 		"Jet_eta");
  AddBranch(&Jet_phi,                 		"Jet_phi");
  AddBranch(&Jet_energy,             		"Jet_energy");
  AddBranch(&Jet_bDiscriminator,      		"Jet_bDiscriminator");
  AddBranch(&Jet_bDiscriminator_CISVV2,      	"Jet_bDiscriminator_CISVV2");
  AddBranch(&Jet_bDiscriminator_pfCISVV2,      	"Jet_bDiscriminator_pfCISVV2");
  //AddBranch(&Jet_vtxMass,            		"Jet_vtxMass");
  //AddBranch(&Jet_decayLength,            	"Jet_decayLength");
  //AddBranch(&Jet_decayLengthSignificance,	"Jet_decayLengthSignificance");
  //AddBranch(&Jet_pileupId,            		"Jet_pileupId");
  AddBranch(&Jet_partonFlavour,       		"Jet_partonFlavour");
  AddBranch(&Jet_neutralHadEnergyFraction,    	"Jet_neutralHadEnergyFraction");
  AddBranch(&Jet_neutralEmEmEnergyFraction,   	"Jet_neutralEmEmEnergyFraction");
  AddBranch(&Jet_chargedHadronEnergyFraction, 	"Jet_chargedHadronEnergyFraction");
  AddBranch(&Jet_chargedEmEnergyFraction,     	"Jet_chargedEmEnergyFraction");
  AddBranch(&Jet_muonEnergyFraction,          	"Jet_muonEnergyFraction");
  AddBranch(&Jet_numberOfConstituents,		"Jet_numberOfConstituents");
  AddBranch(&Jet_chargedMultiplicity, 		"Jet_chargedMultiplicity");
  AddBranch(&Jet_puppi_pt,             		"Jet_puppi_pt");
  AddBranch(&Jet_puppi_eta,              	"Jet_puppi_eta");
  AddBranch(&Jet_puppi_phi,                 	"Jet_puppi_phi");
  AddBranch(&Jet_puppi_energy,             	"Jet_puppi_energy");
  AddBranch(&Jet_puppi_bDiscriminator,      	"Jet_puppi_bDiscriminator");
  AddBranch(&Jet_puppi_bDiscriminator_CISVV2,   "Jet_puppi_bDiscriminator_CISVV2");
  AddBranch(&Jet_puppi_bDiscriminator_pfCISVV2, "Jet_puppi_bDiscriminator_pfCISVV2");
  //AddBranch(&Jet_puppi_vtxMass,       		"Jet_puppi_vtxMass");
  //AddBranch(&Jet_puppi_decayLength,            	"Jet_puppi_decayLength");
  //AddBranch(&Jet_puppi_decayLengthSignificance,	"Jet_puppi_decayLengthSignificance");
  //AddBranch(&Jet_puppi_pileupId,            	"Jet_puppi_pileupId");
  AddBranch(&Jet_puppi_partonFlavour,       	"Jet_puppi_partonFlavour");
  AddBranch(&Jet_puppi_neutralHadEnergyFraction,"Jet_puppi_neutralHadEnergyFraction");
  AddBranch(&Jet_puppi_neutralEmEmEnergyFraction,"Jet_puppi_neutralEmEmEnergyFraction");
  AddBranch(&Jet_puppi_chargedHadronEnergyFraction,"Jet_puppi_chargedHadronEnergyFraction");
  AddBranch(&Jet_puppi_chargedEmEnergyFraction, "Jet_puppi_chargedEmEnergyFraction");
  AddBranch(&Jet_puppi_muonEnergyFraction,      "Jet_puppi_muonEnergyFraction");
  AddBranch(&Jet_puppi_numberOfConstituents,	"Jet_puppi_numberOfConstituents");
  AddBranch(&Jet_puppi_chargedMultiplicity, 	"Jet_puppi_chargedMultiplicity");


  AddBranch(&Jet_AK8_pt,                "Jet_AK8_pt");
  AddBranch(&Jet_AK8_eta,               "Jet_AK8_eta");
  AddBranch(&Jet_AK8_phi,                   "Jet_AK8_phi");
  AddBranch(&Jet_AK8_energy,              "Jet_AK8_energy");
  AddBranch(&Jet_AK8_bDiscriminator,        "Jet_AK8_bDiscriminator");
  AddBranch(&Jet_AK8_bDiscriminator_CISVV2,   "Jet_AK8_bDiscriminator_CISVV2");
  AddBranch(&Jet_AK8_bDiscriminator_pfCISVV2, "Jet_AK8_bDiscriminator_pfCISVV2");
  //AddBranch(&Jet_AK8_vtxMass,           "Jet_AK8_vtxMass");
  //AddBranch(&Jet_AK8_decayLength,             "Jet_AK8_decayLength");
  //AddBranch(&Jet_AK8_decayLengthSignificance, "Jet_AK8_decayLengthSignificance");
  //AddBranch(&Jet_AK8_pileupId,              "Jet_AK8_pileupId");
  AddBranch(&Jet_AK8_partonFlavour,         "Jet_AK8_partonFlavour");
  AddBranch(&Jet_AK8_neutralHadEnergyFraction,"Jet_AK8_neutralHadEnergyFraction");
  AddBranch(&Jet_AK8_neutralEmEmEnergyFraction,"Jet_AK8_neutralEmEmEnergyFraction");
  AddBranch(&Jet_AK8_chargedHadronEnergyFraction,"Jet_AK8_chargedHadronEnergyFraction");
  AddBranch(&Jet_AK8_chargedEmEnergyFraction, "Jet_AK8_chargedEmEnergyFraction");
  AddBranch(&Jet_AK8_muonEnergyFraction,      "Jet_AK8_muonEnergyFraction");
  AddBranch(&Jet_AK8_numberOfConstituents,  "Jet_AK8_numberOfConstituents");
  AddBranch(&Jet_AK8_chargedMultiplicity,   "Jet_AK8_chargedMultiplicity");

  AddBranch(&Jet_AK8_puppi_pt, "Jet_AK8_puppi_pt");

 
  AddBranch(&Jet_AK8_puppi_eta,"Jet_AK8_puppi_eta");
  AddBranch(&Jet_AK8_puppi_mass,"Jet_AK8_puppi_mass");
  AddBranch(&Jet_AK8_puppi_phi,"Jet_AK8_puppi_phi");
  AddBranch(&Jet_AK8_puppi_tau1,"Jet_AK8_puppi_tau1");
  AddBranch(&Jet_AK8_puppi_tau2,"Jet_AK8_puppi_tau2");
  AddBranch(&Jet_AK8_puppi_tau3,"Jet_AK8_puppi_tau3");
  AddBranch(&Jet_AK8_GEN_pt, "Jet_AK8_GEN_pt");
  AddBranch(&Jet_AK8_GEN_phi,"Jet_AK8_GEN_phi");
  AddBranch(&Jet_AK8_GEN_eta,"Jet_AK8_GEN_eta");
  AddBranch(&Jet_AK8_GEN_mass,"Jet_AK8_GEN_mass");
  AddBranch(&Jet_AK8_GEN_energy,"Jet_AK8_GEN_energy");
  AddBranch(&Jet_AK8_GEN_parton,"Jet_AK8_GEN_parton");
  AddBranch(&Jet_AK8_GEN_mother,"Jet_AK8_GEN_mother");

  AddBranch(&Jet_AK8_subjet0_pt, "Jet_AK8_subjet0_pt");
  AddBranch(&Jet_AK8_subjet0_phi, "Jet_AK8_subjet0_phi");
  AddBranch(&Jet_AK8_subjet0_eta, "Jet_AK8_subjet0_eta");
  AddBranch(&Jet_AK8_subjet0_mass, "Jet_AK8_subjet0_mass");
  AddBranch(&Jet_AK8_subjet0_energy, "Jet_AK8_subjet0_energy");
  AddBranch(&Jet_AK8_subjet0_CSVv2, "Jet_AK8_subjet0_CSVv2");

  AddBranch(&Jet_AK8_subjet1_pt, "Jet_AK8_subjet1_pt");
  AddBranch(&Jet_AK8_subjet1_phi, "Jet_AK8_subjet1_phi");
  AddBranch(&Jet_AK8_subjet1_eta, "Jet_AK8_subjet1_eta");
  AddBranch(&Jet_AK8_subjet1_mass, "Jet_AK8_subjet1_mass");
  AddBranch(&Jet_AK8_subjet1_energy, "Jet_AK8_subjet1_energy");
  AddBranch(&Jet_AK8_subjet1_CSVv2, "Jet_AK8_subjet1_CSVv2");

  if(!_super_TNT){
    AddBranch(&Jet_mass,               		"Jet_mass");
    AddBranch(&Jet_electronEnergy,      	"Jet_electronEnergy");
    AddBranch(&Jet_photonEnergy,        	"Jet_photonEnergy");
    AddBranch(&UncorrJet_pt,            	"UncorrJet_pt");
    AddBranch(&Jet_puppi_mass,         		"Jet_puppi_mass");
    AddBranch(&Jet_puppi_electronEnergy,      	"Jet_puppi_electronEnergy");
    AddBranch(&Jet_puppi_photonEnergy,        	"Jet_puppi_photonEnergy");
    AddBranch(&UncorrJet_puppi_pt,            	"UncorrJet_puppi_pt");

    AddBranch(&Jet_AK8_mass,            "Jet_AK8_mass");
    AddBranch(&Jet_AK8_electronEnergy,        "Jet_AK8_electronEnergy");
    AddBranch(&Jet_AK8_photonEnergy,          "Jet_AK8_photonEnergy");
    AddBranch(&UncorrJet_AK8_pt,              "UncorrJet_AK8_pt");
  }

  if(debug_) std::cout << "     JetSelector: Finished setting branches." << std::endl;
}

void JetSelector::Clear(){
  
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  Jet_bDiscriminator.clear();
  Jet_bDiscriminator_CISVV2.clear();
  Jet_bDiscriminator_pfCISVV2.clear();
  //Jet_vtxMass.clear();
  //Jet_decayLength.clear();
  //Jet_decayLengthSignificance.clear();
  //Jet_pileupId.clear();
  Jet_partonFlavour.clear();
  Jet_mass.clear();
  Jet_neutralHadEnergyFraction.clear();
  Jet_neutralEmEmEnergyFraction.clear();
  Jet_chargedHadronEnergyFraction.clear();
  Jet_chargedEmEnergyFraction.clear();
  Jet_muonEnergyFraction.clear();
  Jet_numberOfConstituents.clear();
  Jet_chargedMultiplicity.clear();
  Jet_electronEnergy.clear();
  Jet_photonEnergy.clear();
  UncorrJet_pt.clear();
  Jet_puppi_pt.clear();
  Jet_puppi_eta.clear();
  Jet_puppi_phi.clear();
  Jet_puppi_energy.clear();
  Jet_puppi_bDiscriminator.clear();
  Jet_puppi_bDiscriminator_CISVV2.clear();
  Jet_puppi_bDiscriminator_pfCISVV2.clear();
  //Jet_puppi_vtxMass.clear();
  //Jet_puppi_decayLength.clear();
  //Jet_puppi_decayLengthSignificance.clear();
  //Jet_puppi_pileupId.clear();
  Jet_puppi_partonFlavour.clear();
  Jet_puppi_mass.clear();
  Jet_puppi_neutralHadEnergyFraction.clear();
  Jet_puppi_neutralEmEmEnergyFraction.clear();
  Jet_puppi_chargedHadronEnergyFraction.clear();
  Jet_puppi_chargedEmEnergyFraction.clear();
  Jet_puppi_muonEnergyFraction.clear();
  Jet_puppi_numberOfConstituents.clear();
  Jet_puppi_chargedMultiplicity.clear();
  Jet_puppi_electronEnergy.clear();
  Jet_puppi_photonEnergy.clear();
  UncorrJet_puppi_pt.clear();



  Jet_AK8_pt.clear();
  Jet_AK8_eta.clear();
  Jet_AK8_phi.clear();
  Jet_AK8_energy.clear();
  Jet_AK8_bDiscriminator.clear();
  Jet_AK8_bDiscriminator_CISVV2.clear();
  Jet_AK8_bDiscriminator_pfCISVV2.clear();
  //Jet_AK8_vtxMass.clear();
  //Jet_AK8_decayLength.clear();
  //Jet_AK8_decayLengthSignificance.clear();
  //Jet_AK8_pileupId.clear();
  Jet_AK8_partonFlavour.clear();
  Jet_AK8_mass.clear();
  Jet_AK8_neutralHadEnergyFraction.clear();
  Jet_AK8_neutralEmEmEnergyFraction.clear();
  Jet_AK8_chargedHadronEnergyFraction.clear();
  Jet_AK8_chargedEmEnergyFraction.clear();
  Jet_AK8_muonEnergyFraction.clear();
  Jet_AK8_numberOfConstituents.clear();
  Jet_AK8_chargedMultiplicity.clear();
  Jet_AK8_electronEnergy.clear();
  Jet_AK8_photonEnergy.clear();
  UncorrJet_AK8_pt.clear();
  Jet_AK8_puppi_pt.clear();
  Jet_AK8_subjet0_pt.clear(); 
  Jet_AK8_puppi_eta.clear();
  Jet_AK8_puppi_mass.clear();
  Jet_AK8_puppi_phi.clear();
  Jet_AK8_puppi_tau1.clear();
  Jet_AK8_puppi_tau2.clear();
  Jet_AK8_puppi_tau3.clear();
  Jet_AK8_subjet0_pt.clear();
  Jet_AK8_subjet0_phi.clear();
  Jet_AK8_subjet0_eta.clear();
  Jet_AK8_subjet0_mass.clear();
  Jet_AK8_subjet0_energy.clear();
  Jet_AK8_subjet0_CSVv2.clear();
  Jet_AK8_subjet1_pt.clear();
  Jet_AK8_subjet1_phi.clear();
  Jet_AK8_subjet1_eta.clear();
  Jet_AK8_subjet1_mass.clear();
  Jet_AK8_subjet1_energy.clear();
  Jet_AK8_subjet1_CSVv2.clear();
  Jet_AK8_GEN_pt.clear();  
  Jet_AK8_GEN_phi.clear();
  Jet_AK8_GEN_eta.clear();
  Jet_AK8_GEN_mass.clear();
  Jet_AK8_GEN_energy.clear();
  Jet_AK8_GEN_parton.clear();
  Jet_AK8_GEN_mother.clear();

}
