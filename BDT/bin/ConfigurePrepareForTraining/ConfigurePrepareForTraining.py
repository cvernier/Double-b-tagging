#!/usr/bin/env python
# -*- coding: utf-8 -*

filelist = open ("runPrepareForTraining/filelist", "r")
run_all = open ("runPrepareForTraining.sh", "wa")

for aline in filelist:
    sample = aline.split('.')

    run_file = open ("runPrepareForTraining/runPrepareForTraining"+sample[0]+".c", "wa")

    print >> run_file, "{"
    print >> run_file, "gROOT->ProcessLine(\".L PrepareForTraining.cc+g\");"
    print >> run_file, "PrepareForTraining(\"file:/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch3/src/"+sample[0]+".root\",\"fortraining/fortraining_"+sample[0]+".root\");"
    print >> run_file, "gROOT->ProcessLine(\".q\");"
    print >> run_file, "}"

    run_file.close()

    print >> run_all, "root -l -b  runPrepareForTraining/runPrepareForTraining"+sample[0]+".c &"

run_all.close()
filelist.close()
    
