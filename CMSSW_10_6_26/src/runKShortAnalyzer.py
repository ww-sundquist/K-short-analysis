import FWCore.ParameterSet.Config as cms
process = cms.Process("KShortProcess")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) ) #-1 in Jaebak example
process.source = cms.Source("PoolSource",
  #fileNames = cms.untracked.vstring('file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/06ABA08A-E2D8-3D4B-9538-76A1098320A5.root')
  fileNames = cms.untracked.vstring('file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/06ABA08A-E2D8-3D4B-9538-76A1098320A5.root', 'file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/0AA3317C-FC80-4D46-A791-0EBB2969DDB4.root', 'file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/2517F2DB-3563-A340-91AD-5C6FBB457B63.root', 'file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/32FF633C-0FA3-1B4E-AB3C-A8189254F0483.root')
  #fileNames = cms.untracked.vstring( 'file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/0AA3317C-FC80-4D46-A791-0EBB2969DDB4.root','file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/06ABA08A-E2D8-3D4B-9538-76A1098320A5.root','file:///net/cms11/cms11r0/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2560000/2517F2DB-3563-A340-91AD-5C6FBB457B63.root') 
)
process.kshort = cms.EDAnalyzer('KShortLifetimeAnalyzer',
  tracks = cms.untracked.InputTag("generalTracks"),
  packedPFCandidates = cms.untracked.InputTag("packedPFCandidates"),
  slimmedKshortVertices = cms.untracked.InputTag("slimmedKshortVertices"),
  offlineSlimmedPrimaryVertices = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  slimmedLambdaVertices = cms.untracked.InputTag("slimmedLambdaVertices"),
#  offlineSlimmedSecondaryVertices = cms.untracked.InputTag("offlineSlimmedSecondaryVertices")
  prunedGenParticles = cms.untracked.InputTag("prunedGenParticles"),
  packedGenParticles = cms.untracked.InputTag("packedGenParticles"),
)
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('kshortAnalyze.root')
)
process.p = cms.Path(process.kshort)
