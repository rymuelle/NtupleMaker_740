from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = '80X_Scalar_MonoTop_LO_Mphi-1500_Mchi-100_13TeV_BSM3G'

config.section_('JobType')
config.JobType.psetName = 'miniAODv2_8010_25ns_MC.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['OutTree.root']

config.section_('Data')
config.Data.inputDataset = '/Scalar_MonoTop_LO_Mphi-1500_Mchi-100_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 1
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/rymuelle/'
config.Data.publication = False

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


