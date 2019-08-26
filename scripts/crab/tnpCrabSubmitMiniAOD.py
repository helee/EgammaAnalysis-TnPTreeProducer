from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys

'''
# this will use CRAB client API
from CRABAPI.RawCommand import crabCommand

# talk to DBS to get list of files in this dataset
from dbs.apis.dbsClient import DbsApi
dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')

dataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISpring18MiniAOD-100X_upgrade2018_realistic_v10-v2/MINIAODSIM'
fileDictList=dbs.listFiles(dataset=dataset)

print ("dataset %s has %d files" % (dataset, len(fileDictList)))

# DBS client returns a list of dictionaries, but we want a list of Logical File Names
lfnList = [ dic['logical_file_name'] for dic in fileDictList ]

# this now standard CRAB configuration

from WMCore.Configuration import Configuration
'''

config = config()

submitVersion ="2016MiniAODv3_Trig2"
doEleTree = 'doEleID=False'
doPhoTree = 'doPhoID=False'
doHLTTree = 'doTrigger=True'
#calibEn   = 'useCalibEn=False'

mainOutputDir = '/store/user/helee/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/user/h/helee/tnp_electron/CMSSW_10_2_5/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
#config.Data.allowNonValidInputDataset = False
config.JobType.sendExternalFolder     = True

config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
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

#    config.Data.splitting     = 'FileBased'
#    config.Data.unitsPerJob   = 10
#    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISpring18MiniAOD-100X_upgrade2018_realistic_v10-v2/MINIAODSIM'

    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'mc')
    config.JobType.pyCfgParams  = ['isMC=True',doEleTree,doHLTTree,'GT=94X_mcRun2_asymptotic_v3']
#    config.Data.userInputFiles = lfnList
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 10

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
#    submit(config)

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'
#    submit(config)


    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    config.Data.lumiMask      = '/afs/cern.ch/user/h/helee/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['isMC=False',doEleTree,doHLTTree,'GT=94X_dataRun2_v10']
 
    config.General.requestName  = '17Jul2018_Run2016B'
    config.Data.inputDataset    = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#    submit(config)    

    config.General.requestName  = '17Jul2018_Run2016C'
    config.Data.inputDataset    = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
#    submit(config)    

    config.General.requestName  = '17Jul2018_Run2016D'
    config.Data.inputDataset    = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
#    submit(config)

    config.General.requestName  = '17Jul2018_Run2016E'
    config.Data.inputDataset    = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
#    submit(config) 

    config.General.requestName  = '17Jul2018_Run2016F-vol2'
    config.Data.inputDataset    = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
    submit(config) 

    config.General.requestName  = '17Jul2018_Run2016G-vol2'
    config.Data.inputDataset    = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
#    submit(config) 

    config.General.requestName  = '17Jul2018_Run2016H-vol2'
    config.Data.inputDataset    = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
#    submit(config) 
