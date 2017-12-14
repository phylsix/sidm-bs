import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_1.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_10.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_2.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_3.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_4.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_5.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_6.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_7.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_8.root',
        'root://cmseos.fnal.gov//store/user/wsi/ZpMuMu_1GeV_2e-11GeV/ZpMuMu_1GeV_2e-11GeV_RECO_MINIAOD/171124_154609/0000/step3_inMINIAODSIM_9.root'
    ])
)

process.MuonAnalysis = cms.EDAnalyzer("sidm::MuonAnalysis",
            GenParticleTag_ = cms.untracked.InputTag("prunedGenParticles"),
              PatPkdGenTag_ = cms.untracked.InputTag("packedGenParticles"),
                  PatPfTag_ = cms.untracked.InputTag("packedPFCandidates"),
            PatElectronTag_ = cms.untracked.InputTag("slimmedElectrons"),
                PatMuonTag_ = cms.untracked.InputTag("slimmedMuons"),
                 PatJetTag_ = cms.untracked.InputTag("slimmedJets"),
)
process.TFileService = cms.Service("TFileService",
                                   fileName
                                   = cms.string("$CMSSW_BASE/src/sidm-bs/jetComponent/mydata/MuonAnalysis.root"))

process.p = cms.Path(process.MuonAnalysis)
