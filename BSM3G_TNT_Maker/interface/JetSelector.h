// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __JET_MU_H_
#define __JET_MU_H_

#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>                                                                    
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "baseTree.h"

using namespace std;
using namespace pat;
using namespace edm;

//
// class declaration
//

class JetSelector : public  baseTree{

public:

  JetSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~JetSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
 private:
  JetSelector(){};
  
  vector <double> Jet_pt, Jet_eta, Jet_phi,Jet_energy,  Jet_bDiscriminator,  Jet_mass;
  vector <double> Jet_bDiscriminator_pfCISVV2, Jet_bDiscriminator_CISVV2, Jet_vtxMass, Jet_decayLength, Jet_decayLengthSignificance;
  vector <double> Jet_pileupId, JetIDPU, Jetpass_pileupJetId, Jet_neutralHadEnergyFraction, Jet_neutralEmEmEnergyFraction; 
  vector <double> Jet_chargedHadronEnergyFraction, Jet_chargedEmEnergyFraction, Jet_muonEnergyFraction; 
  vector <double> Jet_electronEnergy,Jet_photonEnergy, UncorrJet_pt; 
  vector <int> Jet_numberOfConstituents;
  vector <int> Jet_chargedMultiplicity, Jet_partonFlavour;

  vector <double> Jet_puppi_pt, Jet_puppi_eta, Jet_puppi_phi, Jet_puppi_energy, Jet_puppi_bDiscriminator, Jet_puppi_mass;
  vector <double> Jet_puppi_bDiscriminator_pfCISVV2, Jet_puppi_bDiscriminator_CISVV2, Jet_puppi_vtxMass, Jet_puppi_decayLength, Jet_puppi_decayLengthSignificance;
  vector <double> Jet_puppi_pileupId, Jet_puppi_neutralHadEnergyFraction, Jet_puppi_neutralEmEmEnergyFraction; 
  vector <double> Jet_puppi_chargedHadronEnergyFraction, Jet_puppi_chargedEmEnergyFraction, Jet_puppi_muonEnergyFraction; 
  vector <double> Jet_puppi_electronEnergy, Jet_puppi_photonEnergy, UncorrJet_puppi_pt; 
  vector <int> Jet_puppi_numberOfConstituents;
  vector <int> Jet_puppi_chargedMultiplicity, Jet_puppi_partonFlavour;


  vector <double> Jet_AK8_pt, Jet_AK8_eta, Jet_AK8_phi, Jet_AK8_energy, Jet_AK8_bDiscriminator, Jet_AK8_mass;
  vector <double> Jet_AK8_bDiscriminator_pfCISVV2, Jet_AK8_bDiscriminator_CISVV2, Jet_AK8_vtxMass, Jet_AK8_decayLength, Jet_AK8_decayLengthSignificance;
  vector <double> Jet_AK8_pileupId, Jet_AK8_neutralHadEnergyFraction, Jet_AK8_neutralEmEmEnergyFraction; 
  vector <double> Jet_AK8_chargedHadronEnergyFraction, Jet_AK8_chargedEmEnergyFraction, Jet_AK8_muonEnergyFraction; 
  vector <double> Jet_AK8_electronEnergy, Jet_AK8_photonEnergy, UncorrJet_AK8_pt; 
  vector <int> Jet_AK8_numberOfConstituents;
  vector <int> Jet_AK8_chargedMultiplicity, Jet_AK8_partonFlavour;

  
  vector <double> Jet_AK8_tau2, Jet_AK8_tau3;
  vector <double> Jet_AK8_puppi_pt, Jet_AK8_puppi_energy,Jet_AK8_puppi_phi, Jet_AK8_puppi_eta, Jet_AK8_puppi_mass, Jet_AK8_puppi_tau1, Jet_AK8_puppi_tau2, Jet_AK8_puppi_tau3;
  vector <double>  Jet_AK8_subjet0_pt, Jet_AK8_subjet1_pt, Jet_AK8_subjet0_eta, Jet_AK8_subjet1_eta, Jet_AK8_subjet0_phi, Jet_AK8_subjet1_phi, Jet_AK8_subjet0_energy, Jet_AK8_subjet1_energy, Jet_AK8_subjet0_mass, Jet_AK8_subjet1_mass, Jet_AK8_subjet0_CSVv2, Jet_AK8_subjet1_CSVv2 ;

  vector <double> Jet_AK8_GEN_pt, Jet_AK8_GEN_energy,Jet_AK8_GEN_phi, Jet_AK8_GEN_eta, Jet_AK8_GEN_mass, Jet_AK8_GEN_parton,Jet_AK8_GEN_mother;


  /*


   fChain->SetBranchStatus("jetAK8Puppi_size"             , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_Pt"             , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_Eta"            , 1 ); 
   fChain->SetBranchStatus("jetAK8Puppi_Phi"            , 1 ); 
   fChain->SetBranchStatus("jetAK8Puppi_E"            , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_CSVv2"            , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_GenJetE"            , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_vSubjetIndex0"            , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_vSubjetIndex1"            , 1 ); 
   fChain->SetBranchStatus("jetAK8Puppi_tau2"            , 1 );  
   fChain->SetBranchStatus("jetAK8Puppi_tau3"            , 1 ); 
   fChain->SetBranchStatus("jetAK8Puppi_filteredMass"            , 1 ); 

   fChain->SetBranchStatus("jetAK4CHS_size"             , 1 ); 
   fChain->SetBranchStatus("jetAK4CHS_Pt"             , 1 );  
   fChain->SetBranchStatus("jetAK4CHS_Eta"            , 1 ); 
   fChain->SetBranchStatus("jetAK4CHS_Phi"            , 1 ); 
   fChain->SetBranchStatus("jetAK4CHS_E"            , 1 );  
   fChain->SetBranchStatus("jetAK4CHS_CSVv2"            , 1 );



   fChain->SetBranchStatus("subjetAK8Puppi_size"            , 1 );
   fChain->SetBranchStatus("subjetAK8Puppi_Pt"            , 1 );
   fChain->SetBranchStatus("subjetAK8Puppi_Eta"            , 1 );
   fChain->SetBranchStatus("subjetAK8Puppi_Phi"            , 1 );
   fChain->SetBranchStatus("subjetAK8Puppi_CSVv2"            , 1 );
   fChain->SetBranchStatus("subjetAK8Puppi_GenJetE"            , 1 );

   fChain->SetBranchStatus("mu_size"            , 1 );
   fChain->SetBranchStatus("mu_Pt"            , 1 );
   fChain->SetBranchStatus("mu_Eta"            , 1 );
   fChain->SetBranchStatus("mu_Phi"            , 1 );





  */


  bool _super_TNT;
  // Jet cuts
//  edm::InputTag jetToken_;
//  edm::InputTag puppi_jetToken_;
//  edm::InputTag _vertexInputTag;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::JetCollection> puppi_jetToken_;
  edm::EDGetTokenT<pat::JetCollection> AK8_puppi_jetToken_;
  edm::EDGetTokenT<reco::VertexCollection> _vertexInputTag;

  double _Jet_pt_min;
};

#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
