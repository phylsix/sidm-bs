import FWCore.ParameterSet.Config as cms

process = cms.Process("electronFinder")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_1.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_10.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_2.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_3.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_4.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_5.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_6.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_7.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_8.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-13GeV/ZpEE-100MeV_e-13GeV_RECO_MINIAOD/170822_015241/0000/step3_inMINIAODSIM_9.root'
    ])
)

process.electronFinder = cms.EDAnalyzer("sidm::electronFinder",
            GenParticleTag_ = cms.untracked.InputTag('prunedGenParticles'),
            PatElectronTag_ = cms.untracked.InputTag("slimmedElectrons"),
                 PatJetTag_ = cms.untracked.InputTag("slimmedJets"),
                  PatPfTag_ = cms.untracked.InputTag("packedPFCandidates"),
              PatPkdGenTag_ = cms.untracked.InputTag("packedGenParticles")
)
process.TFileService = cms.Service("TFileService",
                                   fileName
                                   = cms.string("$CMSSW_BASE/src/sidm-bs/jetComponent/mydata/electronFinder_100MeV_10E-13GeV.root"))

process.p = cms.Path(process.electronFinder)
