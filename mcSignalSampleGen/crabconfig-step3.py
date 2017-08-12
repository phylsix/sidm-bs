from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'RECO_MINIAOD-SIDM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3.py'
config.Data.inputDataset = '/ZpEE-1GeV-1cm/wsi-ZpEE-1GeV-1cm_DIGI_RAW_HLT-dd9e1b8e9f1afcca61a04ea517251c8c/USER'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1000
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outputDatasetTag = 'ZpEE-1GeV-1cm_RECO_MINIAOD'
config.Data.ignoreLocality = True

config.Site.storageSite ='T3_US_FNALLPC' 
