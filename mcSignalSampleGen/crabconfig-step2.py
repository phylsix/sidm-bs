from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DIGI_RAW_HLT-SIDM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step2.py'
config.Data.inputDataset = '/ZpEE_100MeV_e-13GeV/wsi-ZpEE_100MeV_e-13GeV_GEN_SIM-0cd17f6dccc217863c4c100da1ad457f/USER'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1000
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outputDatasetTag = 'ZpEE_100MeV_e-13GeV_DIGI_RAW_HLT'
config.Data.ignoreLocality = True

config.Site.storageSite ='T3_US_FNALLPC'
