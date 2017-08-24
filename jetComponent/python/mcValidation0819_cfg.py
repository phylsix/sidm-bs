import FWCore.ParameterSet.Config as cms

process = cms.Process("mcValidation")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_1.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_10.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_2.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_3.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_4.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_5.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_6.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_7.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_8.root',
        'root://cms-xrd-global.cern.ch//store/user/wsi/ZpEE_100MeV_e-12GeV/ZpEE-100MeV_e-12GeV_RECO_MINIAOD/170819_203908/0000/step3_inMINIAODSIM_9.root'
    ])
)

mcValidationParams = cms.PSet(
    ### parameter set

    # The mass window would be [ (1-ZpMassSideBand)*M, (1+ZpMassSideBand)*M ]
    ZpMassSideBand = cms.untracked.double(0.1) 
)

process.mcValidation = cms.EDAnalyzer('sidm::mcValidation',
                                      mcValidationParams,

                                      GenParticleTag_ = cms.untracked.InputTag('prunedGenParticles'),
                                      SsVerticeTag_ = cms.untracked.InputTag('slimmedSecondaryVertices'),
                                      PatElectronTag_ = cms.untracked.InputTag('slimmedElectrons'),
                                      PatJetTag_ = cms.untracked.InputTag('slimmedJets')
)
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("$CMSSW_BASE/src/sidm-bs/jetComponent/mydata/mcvalidation0819.root"))

process.p = cms.Path(process.mcValidation)
