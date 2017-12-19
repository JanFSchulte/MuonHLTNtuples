from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple.root']#, 'DQMIO.root']

config.Data.unitsPerJob     = 5000
config.Data.totalUnits      = 400000
config.Data.splitting       = 'EventAwareLumiBased'

config.Data.useParent       = True #!!!!
#config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T2_CH_CERN'
config.JobType.numCores     = 4
config.JobType.maxMemoryMB  = 4000
if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    tag = 'IterL3_V3'

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea   = tag
    config.Data.outLFNDirBase = '/store/group/phys_muon/jschulte/13TeV/MC/' + tag

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    datasets = {}
    
#    datasets['DY']        = ('/DYToLL_M_1_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/AODSIM','muonHLTMCRelaxedPixel.py')
    datasets['DY']        = ('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','iterL3MC.py')
   # datasets['data']        = ('/SingleMuon/Run2017F-ZMu-PromptReco-v1/RAW-RECO','muonHLTData.py')
    #datasets['dataD']        = ('/SingleMuon/Run2017D-ZMu-PromptReco-v1/RAW-RECO','muonHLTData.py')
    datasets['DY_Zpt']    = ('/DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','iterL3MC.py')    
    datasets['Displaced'] = ('/DisplacedSUSY_StopToBL_M-400_CTau-10_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM','iterL3MC.py')
    datasets['JPsi']      = ('/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/AODSIM', 'iterL3MC.py')

    for k, v in datasets.iteritems():
        print v[0]
        config.JobType.psetName    = v[1]
        config.General.requestName = k
        config.Data.inputDataset   = v[0]
        config.Data.outputDatasetTag   = 'IterL3_V3_'+k
        print 'submitting config:'
        print config
        submit(config)


