from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "HN_EleTnP_80X_Trig_tagHEfilter"
doEleTree = 'doEleID=False'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=True'
#calibEn   = 'CalibEn=False'

#mainOutputDir = '/store/user/%s/%s' % (getUsernameFromSiteDB(),submitVersion)
mainOutputDir = '/store/user/helee/%s' % submitVersion
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
    config.Data.unitsPerJob   = 8
#    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['isMC=True',doEleTree,doPhoTree,doHLTTree,'GT=80X_mcRun2_asymptotic_2016_TrancheIV_v6']
#    config.JobType.pyCfgParams  = ['isMC=True',doEleTree,doPhoTree,doHLTTree,'GT=80X_mcRun2_asymptotic_end2016_forEGM_v0']   #only for DYToEE

    config.General.requestName  = 'DYToEE_NNPDF30_13TeV-powheg-pythia8'
    config.Data.inputDataset    = '/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-EGM0_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

 #   sys.exit(0)

    
    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    submit(config)

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    submit(config)
    
    config.General.requestName  = 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8'
    config.Data.inputDataset    = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'
    config.Data.inputDataset    = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'
    config.Data.inputDataset    = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'WW_TuneCUETP8M1_13TeV-pythia8'
    config.Data.inputDataset    = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'WZ_TuneCUETP8M1_13TeV-pythia8'
    config.Data.inputDataset    = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'ZZ_TuneCUETP8M1_13TeV-pythia8'
    config.Data.inputDataset    = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    submit(config)    

    config.General.requestName  = 'GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    submit(config)

#    sys.exit(0)

    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    config.Data.lumiMask      = '/afs/cern.ch/user/h/helee/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['isMC=False',doEleTree,doPhoTree,doHLTTree,'GT=80X_dataRun2_2016SeptRepro_v7']
 
    config.General.requestName  = '03Feb2017_RunB'
    config.Data.inputDataset    = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunC'
    config.Data.inputDataset    = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunD'
    config.Data.inputDataset    = '/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunE'
    config.Data.inputDataset    = '/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunF'
    config.Data.inputDataset    = '/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunG'
    config.Data.inputDataset    = '/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD'
    submit(config)

    config.JobType.pyCfgParams  = ['isMC=False',doEleTree,doPhoTree,doHLTTree,'GT=80X_dataRun2_Prompt_v16']
    config.General.requestName  = '03Feb2017_RunH_v2'
    config.Data.inputDataset    = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    submit(config)
    config.General.requestName  = '03Feb2017_RunH_v3'
    config.Data.inputDataset    = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    submit(config)


#/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-EGM0_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM (GT : 80X_mcRun2_asymptotic_end2016_forEGM_v0)
#/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM 
#/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM 
#/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
