import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")

#needs to be commented out or else it complains about: 'two EventSetup Producers want to deliver type="CaloSubdetectorGeometry" label="ZDC"'
#process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50ns'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#Spring16_25nsV6_MC

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #'/store/mc/RunIISpring16MiniAODv2/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/00A945CC-6236-E611-A06B-0002C94CD11E.root'
    #'/store/mc/RunIISpring16MiniAODv2/ZprimeToTauTau_M-1750_TuneCUETP8M1_13TeV-pythia8-tauola/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/02AA02E0-1843-E611-B889-00266CFEFE08.root'
	'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/Scalar_MonoTop_LO_Mphi-2700_Mchi-100_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/62E4658F-B3CD-E611-B825-00259073E412.root'
  )
)

# START ELECTRON ID SECTION
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)

process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
    # boolean variables
    debug_ = cms.bool(False),
    filltriggerinfo     = cms.bool(True),
    fillmuoninfo        = cms.bool(True),
    fillelectronpatinfo = cms.bool(True),
    filltauinfo         = cms.bool(True),
    fillgeninfo         = cms.bool(True),
    fillgenweightinfo   = cms.bool(True),
    fillruninfo         = cms.bool(True),
    fillPVinfo          = cms.bool(True),
    filljetinfo         = cms.bool(True),
    fillMETinfo         = cms.bool(True),
    fillphotoninfo      = cms.bool(True),   

    # make a super tiny ntuple, only with a few branches?
    super_TNT  = cms.bool(False),
    # is data or MC?
    is_data     = cms.bool(False),
 
    # input tags 
    triggerResults      = cms.InputTag( 'TriggerResults', '', 'HLT' ),
    triggerPrescales    = cms.InputTag("patTrigger"),
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot            = cms.InputTag("offlineBeamSpot"),
    pileupInfo          = cms.InputTag("slimmedAddPileupInfo"),
    genparts		= cms.InputTag("prunedGenParticles"),
    genweights		= cms.InputTag("generator"),
    muons               = cms.InputTag("slimmedMuons"),
    patElectrons        = cms.InputTag("slimmedElectrons"),
    taus                = cms.InputTag("slimmedTaus"),
    photons             = cms.InputTag("slimmedPhotons"),
    jets                = cms.InputTag("slimmedJets"),
    jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
    mets                = cms.InputTag("slimmedMETs"),
    metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
    metsNoHF            = cms.InputTag("slimmedMETsNoHF"),
    metsCov             = cms.InputTag("METSignificance","METCovariance"),
    packedPFCandidates  = cms.InputTag("packedPFCandidates"),
    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    electronMVA_wp1     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
    electronMVA_wp2     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),

    # muon cuts
    Muon_pt_min              = cms.double(8.0),
    Muon_eta_max             = cms.double(2.4),
    Muon_vtx_ndof_min        = cms.int32(4),
    Muon_vtx_rho_max         = cms.int32(2),
    Muon_vtx_position_z_max  = cms.double(24.),

    # electron cuts
    patElectron_pt_min       = cms.double(10.0),
    patElectron_eta_max      = cms.double(2.4),
    patElectron_vtx_ndof_min        = cms.int32(4),
    patElectron_vtx_rho_max         = cms.int32(2),
    patElectron_vtx_position_z_max  = cms.double(24.),

    # tau cuts
    Tau_pt_min  = cms.double(20.),
    Tau_eta_max = cms.double(2.3),
    Tau_vtx_ndof_min        = cms.int32(4),
    Tau_vtx_rho_max         = cms.int32(2),
    Tau_vtx_position_z_max  = cms.double(24.),

    # jet cuts
    Jet_pt_min   = cms.double(30.),

    # photon cuts 
    Photon_pt_min   = cms.double(10.0),
    Photon_eta_max  = cms.double(5.0),    

    # primary vertex cuts
    Pvtx_ndof_min   = cms.double(4.),
    Pvtx_vtx_max  = cms.double(24.),
    Pvtx_vtxdxy_max = cms.double(24.),

)

# include bad muon filter
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

# include bad charged hadron filter
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

# met filters using triggerResults
process.load("NtupleMaker.BSM3G_TNT_Maker.triggerReq_cfi")

# re-calculate MET significance and the covariance matrix
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

#process.p = cms.Path(process.goodVerticesFilterRECO * 
#                     process.CSCTightHaloFilterRECO *
#                     process.eeBadScFilterRECO *
#                     process.EcalDeadCellTriggerPrimitiveFilterRECO *
#                     process.HBHENoiseFilterRECO * 
#                     process.HBHENoiseIsoFilterRECO * 
#		     process.BadPFMuonFilter *
#		     process.BadChargedCandidateFilter *
#                     process.egmGsfElectronIDSequence * 
#                     process.METSignificance * 
#                     process.TNT
#)

process.p = cms.Path(process.goodVerticesFilterPAT * 
                     process.CSCTightHaloFilterPAT *
                     process.eeBadScFilterPAT *
                     process.EcalDeadCellTriggerPrimitiveFilterPAT *
                     process.HBHENoiseFilterPAT * 
                     process.HBHENoiseIsoFilterPAT * 
		     process.BadPFMuonFilter *
		     process.BadChargedCandidateFilter *
                     process.egmGsfElectronIDSequence * 
                     process.METSignificance * 
                     process.TNT
)

#process.p = cms.Path(process.TNT)
