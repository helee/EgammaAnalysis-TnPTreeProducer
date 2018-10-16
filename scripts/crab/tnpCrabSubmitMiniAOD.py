from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "HN_EleTnP_80X_ID_Trig"
doEleTree = 'doEleID=True'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=True'
calibEn   = 'useCalibEn=False'

mainOutputDir = '/store/user/%s/%s' % (getUsernameFromSiteDB(),submitVersion)

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/user/h/helee/CMSSW_8_0_27/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
config.Data.allowNonValidInputDataset = True
config.JobType.sendExternalFolder     = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 

config.Site.storageSite = 'T2_KR_KNU'



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    ##### submit MC
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'mc')
    config.Data.splitting     = 'FileBased'
#    config.Data.unitsPerJob   = 8
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['isMC=True',doEleTree,doPhoTree,doHLTTree,calibEn,'GT=80X_mcRun2_asymptotic_2016_TrancheIV_v6']

#    config.General.requestName  = 'ttbar_madgraph'
#    config.Data.inputDataset    = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

 #   sys.exit(0)

    
    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    submit(config)

#    config.General.requestName  = 'WJets_madgraph'
#    config.Data.inputDataset    = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM'
#    submit(config)


#    config.General.requestName  = 'DYToEE_powheg_m50_120'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)
    
#    config.General.requestName  = 'DYToEE_powheg_m120_200'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)
    
#    config.General.requestName  = 'DYToEE_powheg_m200_400'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

#    config.General.requestName  = 'DYToEE_powheg_m400_800'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

#    config.General.requestName  = 'DYToEE_powheg_m800_1400'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)
    
#    config.General.requestName  = 'DYToEE_powheg_m1400_2300'
#    config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

#    sys.exit(0)

    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    config.Data.lumiMask      = '/afs/cern.ch/user/h/helee/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['isMC=False',doEleTree,doPhoTree,doHLTTree,calibEn,'GT=80X_dataRun2_2016LegacyRepro_v4']
 
    config.General.requestName  = '07Aug17_RunB'
    config.Data.inputDataset    = '/SingleElectron/Run2016B-07Aug17_ver2-v2/MINIAOD'
#    submit(config)
    config.General.requestName  = '07Aug17_RunC'
    config.Data.inputDataset    = '/SingleElectron/Run2016C-07Aug17-v1/MINIAOD'
#    submit(config)
    config.General.requestName  = '07Aug17_RunD'
    config.Data.inputDataset    = '/SingleElectron/Run2016D-07Aug17-v1/MINIAOD'
#    submit(config)
    config.General.requestName  = '07Aug17_RunE'
    config.Data.inputDataset    = '/SingleElectron/Run2016E-07Aug17-v1/MINIAOD'
#    submit(config)
    config.General.requestName  = '07Aug17_RunF'
    config.Data.inputDataset    = '/SingleElectron/Run2016F-07Aug17-v1/MINIAOD'
#    submit(config)

#    config.JobType.pyCfgParams  = ['isMC=False',doEleTree,doPhoTree,doHLTTree,calibEn,'GT=80X_dataRun2_Prompt_v16']
    config.General.requestName  = '07Aug17_RunG'
    config.Data.inputDataset    = '/SingleElectron/Run2016G-07Aug17-v1/MINIAOD'
#    submit(config)
    config.General.requestName  = '07Aug17_RunH'
    config.Data.inputDataset    = '/SingleElectron/Run2016H-07Aug17-v1/MINIAOD'
#    submit(config)


