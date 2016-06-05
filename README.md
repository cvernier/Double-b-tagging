# Double-b-tagging

bin/ contains the code to process the BTagAnalyzer output and produce flat tree containing the input variables per jet
It also applies basic selection cut:
- AK08 jets
  pt>300 GeV
  70 < mj_pruned < 200 GeV
- CA15 jets
  pt>170 GeV
  70 < mj_pruned < 200 GeV

Training/ contains all the code needed to perform the training and evaluate the performance.
Plots macro are also stored here. 
