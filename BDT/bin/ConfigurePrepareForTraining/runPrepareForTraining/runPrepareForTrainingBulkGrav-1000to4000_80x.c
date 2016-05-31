{
gROOT->ProcessLine(".L PrepareForTraining.cc+g");
PrepareForTraining("/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch3/src/BulkGrav-1000to4000_80x.root","fortraining/fortraining_BulkGrav-1000to4000_80x.root");
gROOT->ProcessLine(".q");
}
