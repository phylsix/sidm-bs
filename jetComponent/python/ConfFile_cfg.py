import FWCore.ParameterSet.Config as cms

process = cms.Process("jetComponent")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
	    'root://cms-xrd-global.cern.ch//store/user/wsi/Zp-1GeV_ct-1cm/TEST-SIDM-MC/170810_045144/0000/step3_inMINIAODSIM_1.root'
    ])
)

jetComponentParams = cms.PSet(
    # parameter set
)

process.jetComponent = cms.EDAnalyzer('sidm::jetComponent',
                                      jetComponentParams,
                                      GenJetTag_1 = cms.untracked.InputTag('slimmedGenJets'),
                                      GenJetTag_2 = cms.untracked.InputTag('slimmedGenJetsAK8'),

                                      PatJetTag_1 = cms.untracked.InputTag('slimmedJets'),
                                      PatJetTag_2 = cms.untracked.InputTag('slimmedJetsAK8'),
                                      PatJetTag_3 = cms.untracked.InputTag('slimmedJetsPuppi'),
                                      PatJetTag_4 = cms.untracked.InputTag('slimmedJetsAK8PFCHSSoftDropPacked'),
                                      PatJetTag_5 = cms.untracked.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked')
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("$CMSSW_BASE/src/sidm-bs/jetComponent/data/output.root"))

process.p = cms.Path(process.jetComponent)
