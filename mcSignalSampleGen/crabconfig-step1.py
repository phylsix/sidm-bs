from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'GEN_SIM-SIDM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'step1.py'
#config.JobType.inputFiles = ['sidm_events.lhe']

config.Data.outputPrimaryDataset = 'ZpEE-1GeV-1cm'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1000
NJOBS = 300  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'ZpEE-1GeV-1cm_GEN_SIM'

config.Site.storageSite = 'T3_US_FNALLPC'
