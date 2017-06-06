import FWCore.ParameterSet.Config as cms

process = cms.Process("treeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/DE395EB0-3CC8-E611-881D-20CF3027A5C5.root'
    )
)

process.treeMaker = cms.EDAnalyzer('TreeMaker',
                                   GenParticleTag = cms.untracked.InputTag("prunedGenParticles")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root"))

process.p = cms.Path(process.treeMaker)
