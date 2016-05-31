from optparse import OptionParser
import sys
import ConfigParser
import ROOT
import os
import multiprocessing
import array
 
training_vars_float = [
#   "SubJet_csv",
  "z_ratio",#*(SV_vtx_pt_0+SV_vtx_pt_1)/(SV_mass_0+SV_mass_1)",
#  "tau_dot",
#  "SV_mass_0",
#  "SV_vtx_EnergyRatio_0",
#  "SV_vtx_EnergyRatio_1",
#  "SV_vtx_deltaR_0", 
  "trackSipdSig_3",
  "trackSipdSig_2",
  "trackSipdSig_1",
  "trackSipdSig_0",
  "trackSipdSig_1_0",
  "trackSipdSig_0_0",
  "trackSipdSig_1_1",
  "trackSipdSig_0_1", 
  #"Track_distance_TwoHighest3DSig",
  "trackSip2dSigAboveCharm_0",
  #"trackSip2dSigAboveCharm_1", 	
  "trackSip2dSigAboveBottom_0",
  "trackSip2dSigAboveBottom_1", 
  "tau0_trackEtaRel_0",
  "tau0_trackEtaRel_1",
  "tau0_trackEtaRel_2",
  "tau1_trackEtaRel_0",
  "tau1_trackEtaRel_1",
  "tau1_trackEtaRel_2",
  "tau_vertexMass_0",
  #"tau_vertexMass_corrected_0",	
  "tau_vertexEnergyRatio_0",
  "tau_vertexDeltaR_0",
  "tau_flightDistance2dSig_0",
  "tau_vertexMass_1",
  #"tau_vertexMass_corrected_1",
  #"TagVarCSV1_vertexEnergyRatio",
  "tau_vertexEnergyRatio_1",
  #"tau_vertexDeltaR_1",
  "tau_flightDistance2dSig_1",
  ]
 
training_vars_int = [
  "jetNTracks",
  #"tau_vertexNsv_1", 
  #"tau_vertexNsv_0", 
  #"jetNTracks_ip",	
 # "TagVarCSV1_jetNTracksEtaRel",
  "nSV",
  #"tau_vertexNtrk_0",
  #"tau_vertexNtrk_1",
  ]


training_vars_float2 = [
   "SubJet_csv",
  "PFLepton_ratio",
  "PFLepton_ptrel",
  "trackSip3dSig_3",
  "trackSip3dSig_2",
  "trackSip3dSig_1",
  "trackSip3dSig_0",
  "TagVarCSV1_trackSip2dSigAboveCharm",
  "trackEtaRel_0",
  "trackEtaRel_1",
  "trackEtaRel_2",
  "TagVarCSV1_vertexMass",
  "TagVarCSV1_vertexEnergyRatio",
  "TagVarCSV1_vertexJetDeltaR",
  "TagVarCSV1_flightDistance2dSig",
  ]

training_vars_int2 = [
  "nSL_3",              
  "TagVarCSV1_jetNTracks",
  "TagVarCSV1_jetNSecondaryVertices",
  "TagVarCSV1_vertexNTracks"
  ]



argv = sys.argv
parser = OptionParser()
parser.add_option("-c", "--categories", dest="categories", default=False, action="store_true",
                              help="train in pt-eta categories")
parser.add_option("-w", "--weight", dest="weight", default=False, action="store_true",
                              help="pt-eta reweight")
parser.add_option("-g", "--gluonsplitting", dest="gluonsplitting", default=False, action="store_true",
                              help="train bb vs. gsp")
parser.add_option("-C", "--charm", dest="charm", default=False, action="store_true",
                              help="train bb vs. charm") 
parser.add_option("-p", "--usePT", dest="usePT", default=False, action="store_true",
                              help="use pT in training") 
parser.add_option("-a", "--useALL", dest="useALL", default=False, action="store_true",
                              help="use all signal samples in training")                       
parser.add_option("-f", "--file", dest="filename",
                  help="write to FILE", metavar="FILE")                                                                                                                                                                               			      			      			      			      
(opts, args) = parser.parse_args(argv)  

def train(bdtoptions):

  outFile = ROOT.TFile('TMVA_%s.root'%opts.filename, 'RECREATE')
  print "Printing output to %s" %outFile.GetName()

  factory = ROOT.TMVA.Factory(
                               "TMVAClassification", 
                               outFile, 
                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification"
                             )
  
  TMVA_tools = ROOT.TMVA.Tools.Instance()
  treeS = ROOT.TChain('Fjets')
  treeB = ROOT.TChain('Fjets')
  
  
############ OLD LUMI-REWEIGHTING SCHEME!! DO NOT USE! Reweight before passing ############
# files = [
#     "../weighted_rootfiles/QCD170-300_forTraining.root"  , "../weighted_rootfiles/QCD300-470_forTraining.root"  ,"../weighted_rootfiles/QCD470-600_forTraining.root","../weighted_rootfiles/QCD600-800_forTraining.root","../weighted_rootfiles/QCD800-1000_forTraining.root"  ,"../weighted_rootfiles/QCD1000-1400_forTraining.root","../weighted_rootfiles/QCD1400-1800_forTraining.root"
#   ]
#
#   xSec = [12030.,
#             7475.,
#             587.1,
#             167.0,
#             28.25,
#             8.195,
#             0.7346
#           ]
#   genEv = [2001169.,
#            1986177.,
#            2001071.,
#            1997744.,
#            1000065.,
#            500550.0,
#            199627.0,
#          ]
#
#   treelist = []
#   for f in files:
#     print 'Opening file %s' %f
#     file = ROOT.TFile.Open(f)
#     tree = file.Get("Fjets")
#     treelist.append(tree)
#     treeB.Add(f)
#
#   for i in range(0,len(treelist)):
#     weight = xSec[i]/genEv[i]
#     tree = treelist[i]
#     print i
#     factory.AddBackgroundTree(tree, weight)
#     print 'Setting weight to %f (xSec = %2f, #genevents = %i)' %(weight, xSec[i], genEv[i])
#     # tree.SetWeight(weight)

 #  treeB.Add('QCD1000-1400_forTraining.root')
 #  treeB.Add('QCD170-300_forTraining.root')
 #  treeB.Add('QCD300-470_forTraining.root')
 #  treeB.Add('QCD470-600_forTraining.root')
 #  treeB.Add('QCD600-800_forTraining.root')
 #  treeB.Add('QCD800-1000_forTraining.root')
 #  treeB.Add('QCD1000-1400_forTraining.root')

  
  
  
  #if(opts.useALL):
  treeS.Add('bulkGrav800-4000_forTraining.root')
  #else:
  #  treeS.Add('../weighted_rootfiles/r800_forTraining.root')
 
#     treeS.Add('../weighted_rootfiles/r1000_forTraining.root')
#     treeS.Add('../weighted_rootfiles/r1600_forTraining.root')
#     treeS.Add('../weighted_rootfiles/r2000_forTraining.root')

  treeB.Add('qcd_forTraining.root')#weighted_rootfiles/qcd_forTraining.root')
  
  signal_selection = '(flavour==5||flavour==-5) && nbHadrons>1'#abs(flavour==5) && nbHadrons>1' # bb massPruned>80 && massPruned<150
  print "Signal selection = %s" %signal_selection
  
  if(opts.gluonsplitting):
    background_selection = 'abs(flavour==5) && nbHadrons>1' #gsp massPruned>80 && massPruned<150 &&
  elif(opts.charm):
    background_selection = 'massPruned>70 && massPruned<150 && abs(flavour==4)' # charm
  else:
    background_selection = 'massPruned>70 && massPruned<200' # no b
  
  print "Bkg selection = %s" %background_selection
  num_pass = treeS.GetEntries(signal_selection)
  num_fail = treeB.GetEntries(background_selection)

  print 'N events signal', num_pass
  print 'N events background', num_fail
  

  for var in training_vars_float:
    print "Adding variable: %s" %var
    #if var =="SV_vtx_EnergyRatio_0" or var=="SV_vtx_EnergyRatio_1":
#	factory.AddVariable(var, 'F', -1, 10)
#    else:	 	
    factory.AddVariable(var, 'F') # add float variable
  for var in training_vars_int:
    factory.AddVariable(var, 'I') # add integer variable
  
  if(opts.usePT):  
    factory.AddVariable("ptPruned", 'F')
    factory.AddVariable("etaPruned")	
    
  factory.AddSpectator("massPruned")
  factory.AddSpectator("flavour")
  factory.AddSpectator("nbHadrons")
  #factory.AddSpectator("SubJet_csv")
  #factory.AddSpectator("FatJet_csv") 
  
  if not opts.usePT: 
     factory.AddSpectator("ptPruned")
     factory.AddSpectator("etaPruned")	
    
  if (opts.weight):
    # factory.AddSpectator("weight_etaPt")
    factory.SetWeightExpression('weight')
  else:
    factory.SetWeightExpression('1.')
    
  factory.AddSignalTree(treeS, 1.)
  factory.AddBackgroundTree(treeB, 1.)

  
  factory.PrepareTrainingAndTestTree( ROOT.TCut(signal_selection), ROOT.TCut(background_selection), 
     "nTrain_Signal=0::nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:!V" )
       #"nTrain_Signal=30000:nTest_Signal=12000:nTrain_Background=30000:nTest_Background=50000:SplitMode=Random:!V" )
      
  # factory.BookMethod( ROOT.TMVA.Types.kFisher, "Fisher", "!H:!V:Fisher" )
  factory.BookMethod( ROOT.TMVA.Types.kBDT,
                      "BDTG",
	#		"!H:!V:NTrees=750:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=35:MaxDepth=4:PruneMethod=CostComplexity:PruneStrength=3"
			"!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.2:SeparationType=MisClassificationError:nCuts=20:PruneMethod=CostComplexity:PruneStrength=2"
                      # "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedGrad:GradBaggingFraction=0.9:SeparationType=GiniIndex:nCuts=500:NNodesMax=5"
                      #":".join(bdtoptions)
                    )
  
  if (opts.categories):
    
    
    if(opts.usePT):  
      theCat1Vars = "PFLepton_ptrel:z_ratio:tau_dot:SV_mass_0:SV_vtx_EnergyRatio_0:SV_vtx_EnergyRatio_1:PFLepton_IP2D:tau2/tau1:nSV:nSL:ptPruned:etaPruned"
    else:
      theCat1Vars = "PFLepton_ptrel:z_ratio:tau_dot:SV_mass_0:SV_vtx_EnergyRatio_0:SV_vtx_EnergyRatio_1:PFLepton_IP2D:tau2/tau1:nSV:nSL"
    # theCat1Vars =  "massPruned:tau2/tau1:SV_flight2D_0:SV_flight2D_1:SV_flight2DErr_0:SV_flight2DErr_1:SV_totCharge_0:SV_totCharge_1:SV_mass_0:SV_mass_1:SV_vtx_pt_0:SV_vtx_pt_1:SV_vtx_EnergyRatio_0:SV_vtx_EnergyRatio_1:SV_vtx_deltaR_0:SV_vtx_deltaR_1:trackSip3dSig_3:trackPtRel_3:trackEtaRel_0:trackEtaRel_1:trackEtaRel_2:PFLepton_deltaR:PFLepton_ptrel:PFLepton_ratioRel:PFLepton_IP2D:nSV:SV_nTrk_0:SV_nTrk_1:jetNTracksEtaRel:nSL"
    # mcat = factory.BookMethod( ROOT.TMVA.Types.kCategory, "BDTCat4","" )
    # cuts = [
    #   'abs(etaPruned)<=1.4 && ptPruned < 400', 'abs(etaPruned)<=1.4 && ptPruned >= 400', 'abs(etaPruned)>1.4 &&  ptPruned <400', 'abs(etaPruned)>1.4 && ptPruned >= 400'
    # ]
    #
    # for cut in cuts:
    #   print "Training in category %s" %cut
    #   mcat.AddMethod( ROOT.TCut(cut), theCat1Vars, ROOT.TMVA.Types.kBDT, "Category_BDT_4","!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4" )
    
    mcat2 = factory.BookMethod( ROOT.TMVA.Types.kCategory, "BDTCat8","" )  
    cuts2 = [ #'ptPruned <500' ,  'ptPruned >=500'
      'abs(etaPruned)<=2.1 && ptPruned <450', 'abs(etaPruned)<=1.4 && ptPruned >= 450 && ptPruned < 800', 'abs(etaPruned)<=1.4 && ptPruned >= 800',# && ptPruned < 800', 'abs(etaPruned)<=1.4 && ptPruned >= 800',
      'abs(etaPruned)>2.1 && ptPruned <450', 'abs(etaPruned)>1.4 && ptPruned >= 450 && ptPruned < 800', 'abs(etaPruned)>1.4 && ptPruned >= 800',# && ptPruned < 800', 'abs(etaPruned)>1.4 && ptPruned >= 800'
    
      # 'abs(etaPruned)<=1.2 && ptPruned <400', 'abs(etaPruned)<=1.2 && ptPruned >= 400 && ptPruned < 600', 'abs(etaPruned)<=1.2 && ptPruned >= 600 && ptPruned < 800', 'abs(etaPruned)<=1.2 && ptPruned >= 800',
      # 'abs(etaPruned)>1.2 && abs(etaPruned)<=2.1 && ptPruned < 400', 'abs(etaPruned)>1.2 && abs(etaPruned)<=2.1 && ptPruned >= 400 && ptPruned < 600', 'abs(etaPruned)>1.2 && abs(etaPruned)<=2.1 && ptPruned >= 600 && ptPruned < 800', 'abs(etaPruned)>1.2 && abs(etaPruned)<=2.1 && ptPruned >= 800',
      # 'abs(etaPruned)>2.1 && ptPruned <400', 'abs(etaPruned)>2.1 && ptPruned >= 400 && ptPruned < 600', 'abs(etaPruned)>2.1 && ptPruned >= 600 && ptPruned < 800','abs(etaPruned)>2.1 && ptPruned > 800'
    ]
 
    for cut in cuts2:
      print "Training in category %s" %cut
      mcat2.AddMethod( ROOT.TCut(cut), theCat1Vars, ROOT.TMVA.Types.kBDT, "BDTCat8", "!H:!V:NTrees=350:MinNodeSize=2.5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.2:SeparationType=MisClassificationError:nCuts=25:PruneMethod=CostComplexity:PruneStrength=3")
#MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=25:PruneMethod=CostComplexity:PruneBeforeBoost=False:PruneStrength=3:MaxDepth=4" )
  
  
  # (ROOT.TMVA.gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 2
  factory.TrainAllMethods()

  factory.TestAllMethods()

  factory.EvaluateAllMethods()

  outFile.Close()

  # ROOT.gROOT.LoadMacro('$ROOTSYS/tmva/test/TMVAGui.C')
  # ROOT.TMVAGui('TMVA_classification.root')
  # raw_input("Press Enter to continue...")

# # def read(inDirName, inFileName):
def read(inDirName, inFileName):
  print "Reading", inFileName
  print "################################"

  TMVA_tools = ROOT.TMVA.Tools.Instance()

  tree = ROOT.TChain('Fjets')
  tree.Add('%s%s' %(inDirName,inFileName))
  print '%s%s' %(inDirName,inFileName)
  print "################################"
  print "################################"
  print "################################"
  reader = ROOT.TMVA.Reader('TMVAClassification_BDTG')
  reader2 = ROOT.TMVA.Reader('TMVAClassification_BDTG')
  etaPruned = array.array('f',[0])
  ptPruned = array.array('f',[0])
  flavour = array.array('f',[0])
  nbHadrons = array.array('f',[0])
  SubJet_csv=array.array('f',[0])
  FatJet_csv=array.array('f',[0])
  massPruned = array.array('f',[0])
  
  
  reader2.AddSpectator("massPruned", massPruned)
  reader2.AddSpectator("flavour", flavour)
  reader2.AddSpectator("nbHadrons", nbHadrons)
  reader.AddSpectator("massPruned", massPruned)
  reader.AddSpectator("flavour", flavour)
  reader.AddSpectator("nbHadrons", nbHadrons)
  #reader.AddSpectator("SubJet_csv",SubJet_csv)
  #reader.AddSpectator("FatJet_csv",FatJet_csv) 
  if not opts.usePT:
    reader.AddSpectator("ptPruned", ptPruned)
    reader.AddSpectator("etaPruned", etaPruned)	
    reader2.AddSpectator("ptPruned", ptPruned)
    reader2.AddSpectator("etaPruned", etaPruned)
  varDict = {}
  varDict2 = {} 
  for var in training_vars_float:
    varDict[var] = array.array('f',[0])
    reader.AddVariable(var, varDict[var])
  for var in training_vars_int:
    varDict[var] = array.array('f',[0])
    reader.AddVariable(var, varDict[var])
  for var in training_vars_float2:
    varDict2[var] = array.array('f',[0])
    reader2.AddVariable(var, varDict2[var])
  for var in training_vars_int2:
    varDict2[var] = array.array('f',[0])
    reader2.AddVariable(var, varDict2[var])
    
  if(opts.usePT):
    reader.AddVariable("ptPruned",ptPruned)
    reader.AddVariable("etaPruned", etaPruned)	
    #reader2.AddVariable("ptPruned",ptPruned)
    #reader2.AddVariable("etaPruned", etaPruned)	

  
  reader.BookMVA("BDTG","weights/TMVAClassification_BDTG.weights.xml")
  #reader2.BookMVA("BDTboost","weights/boost.weights.xml")
  if(opts.categories):
  	reader.BookMVA("BDTCat8","weights/TMVAClassification_BDTCat8.weights.xml")

  bdtOuts = []
  bdtOuts2 = []
  # bdtOutsCat4 = []
  bdtOutsCat8 = []
  flavours = []
  ptPruneds = []
  etaPruneds = []
  massPruneds = []
  nbHadronss=[] 
  SubJet_csvs =[]
  FatJet_csvs =[]
  hBDTGDisc = ROOT.TH1F("hBDTGDisc","",1000,-5,5)
  hBDTCat8Disc = ROOT.TH1F("hBDTCat8Disc","",1000,-5,5)
  

  for jentry in xrange(tree.GetEntries()):

    ientry = tree.LoadTree(jentry)
    nb = tree.GetEntry(jentry)

    for var in varDict:
      if var.find("tau2/tau1") != -1:
        varDict[var][0] = getattr(tree, "tau2")
        varDict[var][0] /= getattr(tree, "tau1")
      else:
        varDict[var][0] = getattr(tree, var)

    bdtOutput = reader.EvaluateMVA("BDTG")

    for var in varDict2:
      if var.find("tau2/tau1") != -1:
        varDict2[var][0] = getattr(tree, "tau2")
        varDict2[var][0] /= getattr(tree, "tau1")
      else:
        varDict2[var][0] = getattr(tree, var)

    #bdtOutput2 = reader2.EvaluateMVA("BDTboost")
    # bdtOutputCat4 = reader.EvaluateMVA("BDTCat4")
    if(opts.categories):
    	bdtOutputCat8 = reader.EvaluateMVA("BDTCat8")
        bdtOutsCat8.append(bdtOutputCat8)

    flavour = tree.flavour
    bdtOuts.append(bdtOutput)
    #bdtOuts2.append(bdtOutput2)
    # bdtOutsCat4.append(bdtOutputCat4)
    flavours.append(flavour)
    ptPruneds.append(tree.ptPruned)
    etaPruneds.append(tree.etaPruned)
    massPruneds.append(tree.massPruned)
    nbHadronss.append(tree.nbHadrons)	
    SubJet_csvs.append(tree.SubJet_csv)
    FatJet_csvs.append(tree.FatJet_csv)
    hBDTGDisc.Fill(bdtOutput)
    if(opts.categories):
       hBDTCat8Disc.Fill(bdtOutputCat8)

    if jentry%10000 == 0:
      print jentry, bdtOutput, flavour

  writeSmallTree = True

  if writeSmallTree:
    print "Writing small tree"

    BDTG = array.array('f',[0])
   # BDTG2 = array.array('f',[0])
    # BDTCat4 = array.array('f',[0])
    if(opts.categories): 
      BDTCat8 = array.array('f',[0])
    flav = array.array('f',[0])
    etaPruned = array.array('f',[0])
    ptPruned = array.array('f',[0])  
    massPruned = array.array('f',[0])
    nbHadrons = array.array('i',[0])
    SubJet_csv= array.array('f',[0])	
    FatJet_csv= array.array('f',[0]) 
    fout = ROOT.TFile('validation_%s_%s.root'%(opts.filename,inFileName.replace(".root","")), 'RECREATE')
    outTree = ROOT.TTree( 'tree', 'b-tagging training tree' )
    outTree.Branch('BDTG', BDTG, 'BDTG/F')
    #outTree.Branch('BDTG2', BDTG2, 'BDTG2/F')
    # outTree.Branch('BDTCat4', BDTCat4, 'BDTCat4/F')
    if(opts.categories):
    	outTree.Branch('BDTCat8', BDTCat8, 'BDTCat8/F')
    outTree.Branch('flavour', flav, 'flavour/F')
    outTree.Branch('etaPruned', etaPruned, 'etaPruned/F')
    outTree.Branch('ptPruned', ptPruned, 'ptPruned/F')
    outTree.Branch('nbHadrons',nbHadrons,'nbHadrons/I')
    outTree.Branch('SubJet_csv',SubJet_csv,'SubJet_csv/F')
    outTree.Branch('FatJet_csv',FatJet_csv,'FatJet_csv/F')
    outTree.Branch('massPruned', massPruned, 'massPruned/F')


    for i in range(len((bdtOuts))):
      BDTG[0] = bdtOuts[i]
      #BDTG2[0] = bdtOuts2[i]
      # BDTCat4[0] = bdtOutsCat4[i]
      if(opts.categories):
      	BDTCat8[0] = bdtOutsCat8[i]
      flav[0] = flavours[i]
      etaPruned[0] = etaPruneds[i]
      ptPruned[0] = ptPruneds[i]
      massPruned[0] = massPruneds[i]
      nbHadrons[0]=nbHadronss[i]
      SubJet_csv[0]=SubJet_csvs[i]
      FatJet_csv[0]=FatJet_csvs[i]
      if i%10000==0:
        print i, bdtOuts[i], flavours[i]
        # print i, bdtOutsCat4[i], flavours[i]
        if(opts.categories):
		print i, bdtOutsCat8[i], flavours[i]
      outTree.Fill()
      
      # treeout.Write()
    fout.Write()
    hBDTGDisc.Write()
    if(opts.categories):
      hBDTCat8Disc.Write()
      del hBDTCat8Disc
    del hBDTGDisc
    fout.Close()
  print "done", inFileName

def readParallel():

  print "start readParallel()"
  ROOT.gROOT.SetBatch(True)
  parallelProcesses = multiprocessing.cpu_count()

  inDirName="/afs/cern.ch/work/c/cvernier/testRegression/CMSSW_7_2_2_patch2/src/hlt_pileup/PuIdTreeProducer/test_HiggsTagging/codeTraining_june/HiggsTagger/plot-PAS/"
  files = [
    'qcd_forTraining.root',
    #'r1000_forTraining.root',
    #'r1600_forTraining.root',
    #'r2000_forTraining.root',
    #'r800_forTraining.root',
    'grav_ALL.root',	
# #     'QCD1000-1400_forTraining.root',
# #     'QCD120-170_forTraining.root',
# #     'QCD1400-1800_forTraining.root',
# #     'QCD170-300_forTraining.root',
# #     'QCD300-470_forTraining.root',
# #     'QCD470-600_forTraining.root',
# #     'QCD600-800_forTraining.root',
# #     'QCD800-1000_forTraining.root',
    ]
  
  # for inFileName in os.listdir(inDirName):
  #  if inFileName.endswith(".root"):
  #    files.append(inFileName)

  # create Pool
  p = multiprocessing.Pool(parallelProcesses)
  print "Using %i parallel processes" %parallelProcesses

  for f in files:
    # debug
    read(inDirName, f)
     # read(f)
    # break
    # run jobs
    #p.apply_async(read, args = (inDirName, f,))

  p.close()
  p.join()


#"!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4"
if __name__ == '__main__':
    bdtoptions = [ "!H",
                                 "!V",
                                 "NTrees=100",
                                 "MinNodeSize=2.5%",
                                 "BoostType=Grad",
                                 "Shrinkage=0.10",
                                 "UseBaggedBoost",
                                 "GradBaggingFraction=0.5",
                                 "nCuts=20",
                                 "MaxDepth=4",
                               ]
    train(bdtoptions)
    # read("/shome/thaarres/HiggsTagger/rootfiles/", "r800_forTraining.root")
    # trainMultiClass()
    inDirName="/afs/cern.ch/work/c/cvernier/testRegression/CMSSW_7_2_2_patch2/src/hlt_pileup/PuIdTreeProducer/test_HiggsTagging/codeTraining_june/HiggsTagger/plot-PAS/"
    files2 = ['zz1000_forTraining.root']
    #files = ['ttbar_forTraining.root',
    files = ['bulkGrav800-4000_forTraining.root',	
        'qcd_forTraining.root',
        #'r800_forTraining.root',
	#'r1200_forTraining.root',
	#'r1400_forTraining.root',
	#'r1600_forTraining.root',
	#'r600_forTraining.root',
	#'r1800_forTraining.root',
        #'r1000_forTraining.root',
    	#'r1600_forTraining.root',
        #'r2000_forTraining.root', 
#	'zz1000_forTraining.root'
        ]
    for f in files:
      read(inDirName, f)
# readParallel()

