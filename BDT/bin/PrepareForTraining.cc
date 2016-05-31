//================================================================================================
// Input arguments:
//   argv[1] => input file name
//   argv[2] => output file name
//________________________________________________________________________________________________

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TMath.h>

// Other C++ headers
#include <string>
#include <vector>
#include "iostream"
#include "cmath"
#include "cassert"

using namespace std;

double etaRel(TVector3 &dir, TVector3 &track)
     {
             dir = dir.Unit();
             double momPar = dir.Dot(track);
             double energy = std::sqrt(track.Mag2() + 0.13957*0.13957);

            return 0.5 * std::log((energy + momPar) / (energy - momPar));
}



//=== MAIN =======================================================================================================

//int main(int argc, char **argv)
void PrepareForTraining(const std::string input, const std::string output)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================  

  // handle input arguments
//  const std::string infilename   = argv[1];
//  const std::string outfilename  = argv[2];
   const std::string infilename   = input;
  const std::string outfilename  = output;
  

  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================

  TFile *outf  = new TFile(outfilename.c_str(),"RECREATE");
  TTree *trOUT = new TTree("Fjets","Fjets");
  TH1F *h_multiplicity = new TH1F("h_multiplicity", "h_multiplicity", 6, 0, 6);

  int nSV, SV_nTrk_0, SV_nTrk_1;	
  float ptPruned, etaPruned, phiPruned, tau2, tau1, massPruned, tau1IVF, tau2IVF,  deltaR_tau2SV, deltaR_tau1SV,SubJet_csv, tau_dot;
  float  SV_flight_0 ,  SV_flightErr_0 , SV_vtx_pt_0 , SV_flight2D_0 , SV_flight2DErr_0 ,  SV_totCharge_0 ,  SV_mass_0, SV_dirX_1, SV_dirX_0, SV_dirY_1, SV_dirY_0 , SV_dirZ_1, SV_dirZ_0, SV_deltaR_axis_0, SV_deltaR_axis_1;
  float  SV_flight_1 ,  SV_flightErr_1 , SV_vtx_pt_1 , SV_flight2D_1 , SV_flight2DErr_1 ,  SV_totCharge_1 ,  SV_mass_1 ;
  int jetNTracks, jetNTracks_ip          ,jetNTracksEtaRel, jetNSecondaryVertices, vertexNTracks, vertexCategory ;
  float trackSumJetEtRatio  ,trackSumJetDeltaR   ,trackSip2dValAboveCharm , trackSip2dSigAboveCharm , trackSip3dValAboveCharm , trackSip3dSigAboveCharm ,vertexMass      ,vertexEnergyRatio   ,vertexJetDeltaR     ,flightDistance2dVal ,flightDistance2dSig ,flightDistance3dVal ,flightDistance3dSig ;
  float  trackSip2dSig_0,  trackSip2dSig_1,     trackSip2dSig_2,    trackSip2dSig_3,     trackSip3dSig_0,    trackSip3dSig_1,       trackSip3dSig_3,     trackPtRel_0 ,           trackPtRel_1 ,           trackPtRel_2 ,           trackPtRel_3 ,          trackEtaRel_0 , trackEtaRel_1 ,   trackEtaRel_2, trackEtaRel_3 , trackSip3dSig_2, z_ratio=0, dR_SV, z_ratio1=0, z_ratio2=0,tau1_trackEtaRel_0 , tau1_trackEtaRel_1 ,   tau1_trackEtaRel_2 ,  tau0_trackEtaRel_0 , tau0_trackEtaRel_1 ,   tau0_trackEtaRel_2, trackSip2dSigAboveCharm_0, trackSip2dSigAboveCharm_1, trackSip2dSigAboveBottom_0, trackSip2dSigAboveBottom_1;
  float trackSipdSig_0,  trackSipdSig_1,     trackSipdSig_2,    trackSipdSig_3 ;
  float trackSipdSig_0_0,  trackSipdSig_0_1,     trackSipdSig_0_2,    trackSipdSig_0_3 ;
  float trackSipdSig_1_0,  trackSipdSig_1_1,     trackSipdSig_1_2,    trackSipdSig_1_3 ;
  float trackEtaRel_0_0 , trackEtaRel_0_1 ,   trackEtaRel_0_2, trackEtaRel_1_0 , trackEtaRel_1_1 ,   trackEtaRel_1_2;
  float tau01_trackEtaRel_0 , tau01_trackEtaRel_1 ,   tau01_trackEtaRel_2 ,  tau00_trackEtaRel_0 , tau00_trackEtaRel_1 ,   tau00_trackEtaRel_2;
  float DoubleB =-1., BDTG=-1.;	
  int SubJet1_flavour  ,SubJet1_nbHadrons, TagVarCSV1_jetNTracks  ,    TagVarCSV1_jetNTracksEtaRel   ,TagVarCSV1_vertexCategory   ,TagVarCSV1_jetNSecondaryVertices  ,TagVarCSV1_vertexNTracks  ;

  float SubJet1_pt ,SubJet1_eta ,SubJet1_phi , SubJet1_mass  ,SubJet1_JP , SubJet1_JBP ,SubJet1_CSV ,SubJet1_CSVIVF   ,TagVarCSV1_trackSumJetEtRatio  , TagVarCSV1_trackSumJetDeltaR , TagVarCSV1_trackSip2dValAboveCharm, TagVarCSV1_trackSip2dSigAboveCharm, TagVarCSV1_trackSip3dValAboveCharm, TagVarCSV1_trackSip3dSigAboveCharm, TagVarCSV1_vertexMass , TagVarCSV1_vertexEnergyRatio ,TagVarCSV1_vertexJetDeltaR  ,TagVarCSV1_flightDistance2dVal ,TagVarCSV1_flightDistance2dSig ,TagVarCSV1_flightDistance3dVal ,TagVarCSV1_flightDistance3dSig ,TagVarCSV1_trackEtaRel_0 ,TagVarCSV1_trackEtaRel_1 ,TagVarCSV1_trackEtaRel_2 ;

  int SubJet2_flavour  ,SubJet2_nbHadrons, TagVarCSV2_jetNTracks  ,    TagVarCSV2_jetNTracksEtaRel   ,TagVarCSV2_vertexCategory   ,TagVarCSV2_jetNSecondaryVertices  ,TagVarCSV2_vertexNTracks , nbHadrons, flavour ;
  float PFLepton_deltaR=-1, PFLepton_ptrel=-1, PFLepton_ratioRel=-1,PFLepton_IP2D=-1,SV_vtx_EnergyRatio_1=-1,SV_vtx_EnergyRatio_0=-1,SV_vtx_deltaR_1,SV_vtx_deltaR_0, PFLepton_pt=-1, PFLepton_ratio=-1, PFLepton_IP=-1;
  float PFLepton_deltaR_0 , PFLepton_ptrel_0, PFLepton_ratioRel_0, PFLepton_IP2D_0;
  float PFLepton_deltaR_1 , PFLepton_ptrel_1, PFLepton_ratioRel_1, PFLepton_IP2D_1;
  float SV_vtx_EnergyRatio_pseudo_1=-1,SV_vtx_EnergyRatio_tighter_1=-1,SV_vtx_EnergyRatio_pseudo_0=-1,SV_vtx_EnergyRatio_tighter_0=-1;
  float SV_mass_correctedT_0=-1,SV_mass_corrected_0=-1,SV_mass_correctedT_1=-1,SV_mass_corrected_1=-1;

  float SubJet2_pt ,SubJet2_eta ,SubJet2_phi , SubJet2_mass  ,SubJet2_JP , SubJet2_JBP ,SubJet2_CSV ,SubJet2_CSVIVF   ,TagVarCSV2_trackSumJetEtRatio  , TagVarCSV2_trackSumJetDeltaR , TagVarCSV2_trackSip2dValAboveCharm, TagVarCSV2_trackSip2dSigAboveCharm, TagVarCSV2_trackSip3dValAboveCharm, TagVarCSV2_trackSip3dSigAboveCharm, TagVarCSV2_vertexMass , TagVarCSV2_vertexEnergyRatio ,TagVarCSV2_vertexJetDeltaR  ,TagVarCSV2_flightDistance2dVal ,TagVarCSV2_flightDistance2dSig ,TagVarCSV2_flightDistance3dVal ,TagVarCSV2_flightDistance3dSig ,TagVarCSV2_trackEtaRel_0 ,TagVarCSV2_trackEtaRel_1 ,TagVarCSV2_trackEtaRel_2 ;
  int cc1=0, cc0 =0 , cc2=0, cc3=0, cc4=0, nsel=0 ;
  int nSL=-1, nSE=-1, nSM=-1; int nSL_1, nSL_2, nSL_3;
  float PFMuon_deltaR=-1, PFMuon_ptrel, PFMuon_ratioRel, PFMuon_IP2D;
  float PFElectron_deltaR=-1, PFElectron_ptrel, PFElectron_ratioRel, PFElectron_IP2D;	
  float FatJet_csv=-1;
  float tau_vertexMass_0=-3.,tau_vertexMass_1=-3.,tau_vertexMass_corrected_0=-3.,tau_vertexMass_corrected_1=-3.,tau_vertexEnergyRatio_1=-3.,tau_vertexEnergyRatio_0=-3.,tau_vertexDeltaR_1=-3.,	tau_vertexDeltaR_0=-3.,Track_distance_TwoHighest3DSig=-3., tau_flightDistance2dSig_1=-3. ,tau_flightDistance2dSig_0=-3.,tau_vertexFlight2DErr_1=-3.,tau_vertexFlight2DErr_0=-3.;
  int tau_vertexNtrk_1=-3, tau_vertexNtrk_0=-3, tau_vertexNsv_1=-3, tau_vertexNsv_0=-3;

  float SV_axis_deltaR, SV_axis_deltaPhi,SV_axis_dot, SV_axis_cross; 

  trOUT->Branch("tau_dot",&tau_dot,"tau_dot/F");

  trOUT->Branch("tau_vertexMass_0", &tau_vertexMass_0, "tau_vertexMass_0/F");
  trOUT->Branch("tau_vertexMass_corrected_0", &tau_vertexMass_corrected_0, "tau_vertexMass_corrected_0/F");
  trOUT->Branch("tau_vertexMass_corrected_1", &tau_vertexMass_corrected_1, "tau_vertexMass_corrected_1/F");
  trOUT->Branch("tau_vertexMass_1", &tau_vertexMass_1, "tau_vertexMass_1/F");

  trOUT->Branch("tau_vertexEnergyRatio_1", &tau_vertexEnergyRatio_1, "tau_vertexEnergyRatio_1/F");
  trOUT->Branch("tau_vertexEnergyRatio_0", &tau_vertexEnergyRatio_0, "tau_vertexEnergyRatio_0/F");

  trOUT->Branch("tau_vertexDeltaR_1", &tau_vertexDeltaR_1, "tau_vertexDeltaR_1/F");
  trOUT->Branch("tau_vertexDeltaR_0", &tau_vertexDeltaR_0, "tau_vertexDeltaR_0/F");
  trOUT->Branch("Track_distance_TwoHighest3DSig", &Track_distance_TwoHighest3DSig, "Track_distance_TwoHighest3DSig/F");


  trOUT->Branch("tau_vertexNtrk_1", &tau_vertexNtrk_1, "tau_vertexNtrk_1/I");
  trOUT->Branch("tau_vertexNtrk_0", &tau_vertexNtrk_0, "tau_vertexNtrk_0/I");
  trOUT->Branch("tau_vertexNsv_1", &tau_vertexNsv_1, "tau_vertexNsv_1/I");
  trOUT->Branch("tau_vertexNsv_0", &tau_vertexNsv_0, "tau_vertexNsv_0/I");
  trOUT->Branch("tau_flightDistance2dSig_1", &tau_flightDistance2dSig_1, "tau_flightDistance2dSig_1/F");
  trOUT->Branch("tau_flightDistance2dSig_0", &tau_flightDistance2dSig_0, "tau_flightDistance2dSig_0/F");
  trOUT->Branch("tau_vertexFlight2DErr_1", &tau_vertexFlight2DErr_1, "tau_vertexFlight2DErr_1/F");
  trOUT->Branch("tau_vertexFlight2DErr_0", &tau_vertexFlight2DErr_0, "tau_vertexFlight2DErr_0/F");

  trOUT->Branch("SV_deltaR_axis_0",&SV_deltaR_axis_0,"SV_deltaR_axis_0/F");
  trOUT->Branch("SV_deltaR_axis_1",&SV_deltaR_axis_1,"SV_deltaR_axis_1/F");
  trOUT->Branch("SubJet_csv",&SubJet_csv,"SubJet_csv/F");
  trOUT->Branch("FatJet_csv",&FatJet_csv,"FatJet_csv/F");
  trOUT->Branch("SubJet2_pt ", &SubJet2_pt, "SubJet2_pt/F");
  trOUT->Branch("SubJet2_eta",   &SubJet2_eta, "SubJet2_eta/F");
  trOUT->Branch("SubJet2_phi      ",   &SubJet2_phi, "SubJet2_phi/F");
  trOUT->Branch("SubJet2_mass     ",   &SubJet2_mass, "SubJet2_mass/F");
  trOUT->Branch("SubJet2_flavour  ",   &SubJet2_flavour, "SubJet2_flavour/I");
  trOUT->Branch("SubJet2_nbHadrons",   &SubJet2_nbHadrons, "SubJet2_nbHadrons/I");
  trOUT->Branch("SubJet2_JP       ",   &SubJet2_JP, "SubJet2_JP/F");
  trOUT->Branch("SubJet2_JBP      ",   &SubJet2_JBP, "SubJet2_JBP/F");
  trOUT->Branch("SubJet2_CSV      ",   &SubJet2_CSV, "SubJet2_CSV/F");
  trOUT->Branch("SubJet2_CSVIVF   ",   &SubJet2_CSVIVF, "SubJet2_CSVIVF/F");
  trOUT->Branch("TagVarCSV2_jetNTracks             ",   &TagVarCSV2_jetNTracks, "TagVarCSV2_jetNTracks/I");
  trOUT->Branch("TagVarCSV2_jetNTracksEtaRel       ",   &TagVarCSV2_jetNTracksEtaRel, "TagVarCSV2_jetNTracksEtaRel/I");
  trOUT->Branch("TagVarCSV2_trackSumJetEtRatio     ",   &TagVarCSV2_trackSumJetEtRatio, "TagVarCSV2_trackSumJetEtRatio/F");
  trOUT->Branch("TagVarCSV2_trackSumJetDeltaR      ",   &TagVarCSV2_trackSumJetDeltaR, "TagVarCSV2_trackSumJetDeltaR/F");
  trOUT->Branch("TagVarCSV2_trackSip2dValAboveCharm",   &TagVarCSV2_trackSip2dValAboveCharm, "TagVarCSV2_trackSip2dValAboveCharm/F");
  trOUT->Branch("TagVarCSV2_trackSip2dSigAboveCharm",   &TagVarCSV2_trackSip2dSigAboveCharm, "TagVarCSV2_trackSip2dSigAboveCharm/F");
  trOUT->Branch("TagVarCSV2_trackSip3dValAboveCharm",   &TagVarCSV2_trackSip3dValAboveCharm, "TagVarCSV2_trackSip3dValAboveCharm/F");
  trOUT->Branch("TagVarCSV2_trackSip3dSigAboveCharm",   &TagVarCSV2_trackSip3dSigAboveCharm, "TagVarCSV2_trackSip3dSigAboveCharm/F");
  trOUT->Branch("TagVarCSV2_vertexCategory         ",   &TagVarCSV2_vertexCategory, "TagVarCSV2_vertexCategory/i");
  trOUT->Branch("TagVarCSV2_jetNSecondaryVertices  ",   &TagVarCSV2_jetNSecondaryVertices, "TagVarCSV2_jetNSecondaryVertices/I");
  trOUT->Branch("TagVarCSV2_vertexMass             ",   &TagVarCSV2_vertexMass, "TagVarCSV2_vertexMass/F");
  trOUT->Branch("TagVarCSV2_vertexNTracks          ",   &TagVarCSV2_vertexNTracks, "TagVarCSV2_vertexNTracks/I");
  trOUT->Branch("TagVarCSV2_vertexEnergyRatio      ",   &TagVarCSV2_vertexEnergyRatio, "TagVarCSV2_vertexEnergyRatio/F");
  trOUT->Branch("TagVarCSV2_vertexJetDeltaR        ",   &TagVarCSV2_vertexJetDeltaR, "TagVarCSV2_vertexJetDeltaR/F");
  trOUT->Branch("TagVarCSV2_flightDistance2dVal    ",   &TagVarCSV2_flightDistance2dVal, "TagVarCSV2_flightDistance2dVal/F");
  trOUT->Branch("TagVarCSV2_flightDistance2dSig    ",   &TagVarCSV2_flightDistance2dSig, "TagVarCSV2_flightDistance2dSig/F");
  trOUT->Branch("TagVarCSV2_flightDistance3dVal    ",   &TagVarCSV2_flightDistance3dVal, "TagVarCSV2_flightDistance3dVal/F");
  trOUT->Branch("TagVarCSV2_flightDistance3dSig    ",   &TagVarCSV2_flightDistance3dSig, "TagVarCSV2_flightDistance3dSig/F");
  trOUT->Branch("TagVarCSV2_trackEtaRel_0 ",   &TagVarCSV2_trackEtaRel_0, "TagVarCSV2_trackEtaRel_0/F");
  trOUT->Branch("TagVarCSV2_trackEtaRel_1 ",   &TagVarCSV2_trackEtaRel_1, "TagVarCSV2_trackEtaRel_1/F");
  trOUT->Branch("TagVarCSV2_trackEtaRel_2 ",   &TagVarCSV2_trackEtaRel_2, "TagVarCSV2_trackEtaRel_2/F");


  trOUT->Branch("SubJet1_pt ", &SubJet1_pt, "SubJet1_pt/F");
  trOUT->Branch("SubJet1_eta",   &SubJet1_eta, "SubJet1_eta/F");
  trOUT->Branch("SubJet1_phi      ",   &SubJet1_phi, "SubJet1_phi/F");
  trOUT->Branch("SubJet1_mass     ",   &SubJet1_mass, "SubJet1_mass/F");
  trOUT->Branch("SubJet1_flavour  ",   &SubJet1_flavour, "SubJet1_flavour/I");
  trOUT->Branch("SubJet1_nbHadrons",   &SubJet1_nbHadrons, "SubJet1_nbHadrons/I");
  trOUT->Branch("SubJet1_JP       ",   &SubJet1_JP, "SubJet1_JP/F");
  trOUT->Branch("SubJet1_JBP      ",   &SubJet1_JBP, "SubJet1_JBP/F");
  trOUT->Branch("SubJet1_CSV      ",   &SubJet1_CSV, "SubJet1_CSV/F");
  trOUT->Branch("SubJet1_CSVIVF   ",   &SubJet1_CSVIVF, "SubJet1_CSVIVF/F");
  trOUT->Branch("TagVarCSV1_jetNTracks             ",   &TagVarCSV1_jetNTracks, "TagVarCSV1_jetNTracks/I");
  trOUT->Branch("TagVarCSV1_jetNTracksEtaRel       ",   &TagVarCSV1_jetNTracksEtaRel, "TagVarCSV1_jetNTracksEtaRel/I");
  trOUT->Branch("TagVarCSV1_trackSumJetEtRatio     ",   &TagVarCSV1_trackSumJetEtRatio, "TagVarCSV1_trackSumJetEtRatio/F");
  trOUT->Branch("TagVarCSV1_trackSumJetDeltaR      ",   &TagVarCSV1_trackSumJetDeltaR, "TagVarCSV1_trackSumJetDeltaR/F");
  trOUT->Branch("TagVarCSV1_trackSip2dValAboveCharm",   &TagVarCSV1_trackSip2dValAboveCharm, "TagVarCSV1_trackSip2dValAboveCharm/F");
  trOUT->Branch("TagVarCSV1_trackSip2dSigAboveCharm",   &TagVarCSV1_trackSip2dSigAboveCharm, "TagVarCSV1_trackSip2dSigAboveCharm/F");
  trOUT->Branch("TagVarCSV1_trackSip3dValAboveCharm",   &TagVarCSV1_trackSip3dValAboveCharm, "TagVarCSV1_trackSip3dValAboveCharm/F");
  trOUT->Branch("TagVarCSV1_trackSip3dSigAboveCharm",   &TagVarCSV1_trackSip3dSigAboveCharm, "TagVarCSV1_trackSip3dSigAboveCharm/F");
  trOUT->Branch("TagVarCSV1_vertexCategory         ",   &TagVarCSV1_vertexCategory, "TagVarCSV1_vertexCategory/i");
  trOUT->Branch("TagVarCSV1_jetNSecondaryVertices  ",   &TagVarCSV1_jetNSecondaryVertices, "TagVarCSV1_jetNSecondaryVertices/I");
  trOUT->Branch("TagVarCSV1_vertexMass             ",   &TagVarCSV1_vertexMass, "TagVarCSV1_vertexMass/F");
  trOUT->Branch("TagVarCSV1_vertexNTracks          ",   &TagVarCSV1_vertexNTracks, "TagVarCSV1_vertexNTracks/I");
  trOUT->Branch("TagVarCSV1_vertexEnergyRatio      ",   &TagVarCSV1_vertexEnergyRatio, "TagVarCSV1_vertexEnergyRatio/F");
  trOUT->Branch("TagVarCSV1_vertexJetDeltaR        ",   &TagVarCSV1_vertexJetDeltaR, "TagVarCSV1_vertexJetDeltaR/F");
  trOUT->Branch("TagVarCSV1_flightDistance2dVal    ",   &TagVarCSV1_flightDistance2dVal, "TagVarCSV1_flightDistance2dVal/F");
  trOUT->Branch("TagVarCSV1_flightDistance2dSig    ",   &TagVarCSV1_flightDistance2dSig, "TagVarCSV1_flightDistance2dSig/F");
  trOUT->Branch("TagVarCSV1_flightDistance3dVal    ",   &TagVarCSV1_flightDistance3dVal, "TagVarCSV1_flightDistance3dVal/F");
  trOUT->Branch("TagVarCSV1_flightDistance3dSig    ",   &TagVarCSV1_flightDistance3dSig, "TagVarCSV1_flightDistance3dSig/F");
  trOUT->Branch("TagVarCSV1_trackEtaRel_0 ",   &TagVarCSV1_trackEtaRel_0, "TagVarCSV1_trackEtaRel_0/F");
  trOUT->Branch("TagVarCSV1_trackEtaRel_1 ",   &TagVarCSV1_trackEtaRel_1, "TagVarCSV1_trackEtaRel_1/F");
  trOUT->Branch("TagVarCSV1_trackEtaRel_2 ",   &TagVarCSV1_trackEtaRel_2, "TagVarCSV1_trackEtaRel_2/F");

  trOUT->Branch("PFLepton_deltaR",&PFLepton_deltaR,"PFLepton_deltaR/F");
  trOUT->Branch("PFLepton_ptrel",&PFLepton_ptrel,"PFLepton_ptrel/F");
  trOUT->Branch("PFLepton_ratioRel",&PFLepton_ratioRel,"PFLepton_ratioRel/F");
  trOUT->Branch("PFLepton_IP2D",&PFLepton_IP2D,"PFLepton_IP2D/F");
  trOUT->Branch("PFLepton_IP",&PFLepton_IP,"PFLepton_IP/F");
  trOUT->Branch("PFLepton_ratio",&PFLepton_ratio,"PFLepton_ratio/F");
  trOUT->Branch("PFLepton_pt",&PFLepton_pt,"PFLepton_pt/F");


  trOUT->Branch("PFLepton_deltaR_0",&PFLepton_deltaR_0,"PFLepton_deltaR_0/F");
  trOUT->Branch("PFLepton_ptrel_0",&PFLepton_ptrel_0,"PFLepton_ptrel_0/F");
  trOUT->Branch("PFLepton_ratioRel_0",&PFLepton_ratioRel_0,"PFLepton_ratioRel_0/F");
  trOUT->Branch("PFLepton_IP2D_0",&PFLepton_IP2D_0,"PFLepton_IP2D_0/F");
  trOUT->Branch("PFLepton_deltaR_1",&PFLepton_deltaR_1,"PFLepton_deltaR_1/F");
  trOUT->Branch("PFLepton_ptrel_1",&PFLepton_ptrel_1,"PFLepton_ptrel_1/F");
  trOUT->Branch("PFLepton_ratioRel_1",&PFLepton_ratioRel_1,"PFLepton_ratioRel_1/F");
  trOUT->Branch("PFLepton_IP2D_1",&PFLepton_IP2D_1,"PFLepton_IP2D_1/F");



  trOUT->Branch("SV_axis_deltaR",  &SV_axis_deltaR ,"SV_axis_deltaR/F");
  trOUT->Branch("SV_axis_deltaPhi", &SV_axis_deltaPhi, "SV_axis_deltaPhi/F");
  trOUT->Branch("SV_axis_dot", &SV_axis_dot , "SV_axis_dot/F");
  trOUT->Branch("SV_axis_cross", &SV_axis_cross, "SV_axis_cross/F"); 
  trOUT->Branch("z_ratio",  &z_ratio ,"z_ratio/F");
  trOUT->Branch("z_ratio1",  &z_ratio1 ,"z_ratio1/F");
  trOUT->Branch("z_ratio2",  &z_ratio2 ,"z_ratio2/F");
  trOUT->Branch("dR_SV", &dR_SV, "dR_SV/F");
  trOUT->Branch("nSL",&nSL ,"nSL/I");
  trOUT->Branch("nSL_1",&nSL_1,"nSL_1/I");
  trOUT->Branch("nSL_2",&nSL_2,"nSL_2/I");
  trOUT->Branch("nSL_3",&nSL_3,"nSL_3/I");
  trOUT->Branch("nSE",&nSE ,"nSE/I");
  trOUT->Branch("nSM",&nSM ,"nSM/I");
  trOUT->Branch("PFElectron_deltaR",&PFElectron_deltaR,"PFElectron_deltaR/F");
  trOUT->Branch("PFElectron_ptrel",&PFElectron_ptrel,"PFElectron_ptrel/F");
  trOUT->Branch("PFElectron_ratioRel",&PFElectron_ratioRel,"PFElectron_ratioRel/F");
  trOUT->Branch("PFElectron_IP2D",&PFElectron_IP2D,"PFElectron_IP2D/F");
  trOUT->Branch("PFMuon_deltaR",&PFMuon_deltaR,"PFMuon_deltaR/F");
  trOUT->Branch("PFMuon_ptrel",&PFMuon_ptrel,"PFMuon_ptrel/F");
  trOUT->Branch("PFMuon_ratioRel",&PFMuon_ratioRel,"PFMuon_ratioRel/F");
  trOUT->Branch("PFMuon_IP2D",&PFMuon_IP2D,"PFMuon_IP2D/F");


  trOUT->Branch("massPruned"  ,&massPruned, "massPruned/F");
  trOUT->Branch("ptPruned"  ,&ptPruned, "ptPruned/F");
  trOUT->Branch("etaPruned"  ,&etaPruned, "etaPruned/F");
  trOUT->Branch("phiPruned"  ,&phiPruned, "phiPruned/F");
  trOUT->Branch("flavour", &flavour, "flavour/I");
  trOUT->Branch("nbHadrons", &nbHadrons, "nbHadrons/I");
  trOUT->Branch("tau2"  ,&tau2, "tau2/F");
  trOUT->Branch("tau1"  ,&tau1, "tau1/F");
  trOUT->Branch("tau2IVF"  ,&tau2IVF, "tau2IVF/F");
  trOUT->Branch("tau1IVF"  ,&tau1IVF, "tau1IVF/F");

  trOUT->Branch("nSV"  ,&nSV, "nSV/I");
  trOUT->Branch("SV_nTrk_0",&SV_nTrk_0, "SV_nTrk_0/I");
  trOUT->Branch("SV_nTrk_1",&SV_nTrk_1, "SV_nTrk_1/I"); 
  trOUT->Branch("SV_flight_0", &SV_flight_0, "SV_flight_0/F");
  trOUT->Branch("SV_flight_1", &SV_flight_1, "SV_flight_1/F");
  trOUT->Branch("SV_flightErr_0", &SV_flightErr_0, "SV_flightErr_0/F");
  trOUT->Branch("SV_flightErr_1", &SV_flightErr_1, "SV_flightErr_1/F");
  trOUT->Branch("SV_flight2D_0", &SV_flight2D_0, "SV_flight2D_0/F");
  trOUT->Branch("SV_flight2D_1", &SV_flight2D_1, "SV_flight2D_1/F");
  trOUT->Branch("SV_flight2DErr_0", &SV_flight2DErr_0, "SV_flight2DErr_0/F");
  trOUT->Branch("SV_flight2DErr_1", &SV_flight2DErr_1, "SV_flight2DErr_1/F");
  trOUT->Branch("SV_totCharge_0", &SV_totCharge_0, "SV_totCharge_0/F");
  trOUT->Branch("SV_totCharge_1", &SV_totCharge_1, "SV_totCharge_1/F");
  trOUT->Branch("SV_mass_0", &SV_mass_0, "SV_mass_0/F");
  trOUT->Branch("SV_mass_correctedT_0", &SV_mass_correctedT_0, "SV_mass_correctedT_0/F");
  trOUT->Branch("SV_mass_corrected_0", &SV_mass_corrected_0, "SV_mass_corrected_0/F");
  trOUT->Branch("SV_mass_correctedT_1", &SV_mass_correctedT_1, "SV_mass_correctedT_1/F");
  trOUT->Branch("SV_mass_corrected_1", &SV_mass_corrected_1, "SV_mass_corrected_1/F");
  trOUT->Branch("SV_mass_1", &SV_mass_1, "SV_mass_1/F");
  trOUT->Branch("SV_vtx_pt_0", &SV_vtx_pt_0, "SV_vtx_pt_0/F");
  trOUT->Branch("SV_vtx_pt_1", &SV_vtx_pt_1, "SV_vtx_pt_1/F");
  trOUT->Branch("SV_vtx_EnergyRatio_1", &SV_vtx_EnergyRatio_1, "SV_vtx_EnergyRatio_1/F");//SV_vtx_deltaR_1
  trOUT->Branch("SV_vtx_EnergyRatio_pseudo_1", &SV_vtx_EnergyRatio_pseudo_1, "SV_vtx_EnergyRatio_pseudo_1/F");
  trOUT->Branch("SV_vtx_EnergyRatio_tighter_1", &SV_vtx_EnergyRatio_tighter_1, "SV_vtx_EnergyRatio_tighter_1/F");
  trOUT->Branch("SV_vtx_EnergyRatio_pseudo_0", &SV_vtx_EnergyRatio_pseudo_0, "SV_vtx_EnergyRatio_pseudo_0/F");
  trOUT->Branch("SV_vtx_EnergyRatio_tighter_0", &SV_vtx_EnergyRatio_tighter_0, "SV_vtx_EnergyRatio_tighter_0/F");

  trOUT->Branch("SV_vtx_EnergyRatio_0", &SV_vtx_EnergyRatio_0, "SV_vtx_EnergyRatio_0/F");
  trOUT->Branch("SV_vtx_deltaR_1", &SV_vtx_deltaR_1, "SV_vtx_deltaR_1/F");
  trOUT->Branch("SV_vtx_deltaR_0", &SV_vtx_deltaR_0, "SV_vtx_deltaR_0/F");
  trOUT->Branch("jetNTracks"  ,&jetNTracks, "jetNTracks/I");
  trOUT->Branch("jetNTracks_ip"  ,&jetNTracks_ip, "jetNTracks_ip/I");

  trOUT->Branch("jetNTracksEtaRel"  ,&jetNTracksEtaRel, "jetNTracksEtaRel/I");
  trOUT->Branch("jetNSecondaryVertices"  ,&jetNSecondaryVertices, "jetNSecondaryVertices/I");
  trOUT->Branch("vertexNTracks"  ,&vertexNTracks, "vertexNTracks/I");
  trOUT->Branch("vertexCategory"  ,&vertexCategory, "vertexCategory/I");
  trOUT->Branch("trackSumJetEtRatio", &trackSumJetEtRatio, "trackSumJetEtRatio/F");
  trOUT->Branch("trackSumJetDeltaR", &trackSumJetDeltaR, "trackSumJetDeltaR/F");
  trOUT->Branch("trackSip2dValAboveCharm", &trackSip2dValAboveCharm, "trackSip2dValAboveCharm/F");
  trOUT->Branch("trackSip2dSigAboveCharm", &trackSip2dSigAboveCharm, "trackSip2dSigAboveCharm/F");
  trOUT->Branch("trackSip3dValAboveCharm", &trackSip3dValAboveCharm, "trackSip3dValAboveCharm/F");
  trOUT->Branch("trackSip3dSigAboveCharm", &trackSip3dSigAboveCharm, "trackSip3dSigAboveCharm/F");
  trOUT->Branch("vertexMass", &vertexMass, "vertexMass/F");
  trOUT->Branch("vertexEnergyRatio", &vertexEnergyRatio, "vertexEnergyRatio/F");
  trOUT->Branch("vertexJetDeltaR", &vertexJetDeltaR, "vertexJetDeltaR/F");
  trOUT->Branch("flightDistance2dVal", &flightDistance2dVal, "flightDistance2dVal/F");
  trOUT->Branch("flightDistance2dSig", &flightDistance2dSig, "flightDistance2dSig/F");
  trOUT->Branch("flightDistance3dVal", &flightDistance3dVal, "flightDistance3dVal/F");
  trOUT->Branch("flightDistance3dSig", &flightDistance3dSig, "flightDistance3dSig/F");
  trOUT->Branch("trackSip2dSig_0", &trackSip2dSig_0, "trackSip2dSig_0/F");
  trOUT->Branch("trackSip2dSig_1", &trackSip2dSig_1, "trackSip2dSig_1/F");
  trOUT->Branch("trackSip2dSig_2", &trackSip2dSig_2, "trackSip2dSig_2/F");
  trOUT->Branch("trackSip2dSig_3", &trackSip2dSig_3, "trackSip2dSig_3/F");
  trOUT->Branch("trackSipdSig_0", &trackSipdSig_0, "trackSipdSig_0/F");
  trOUT->Branch("trackSipdSig_1", &trackSipdSig_1, "trackSipdSig_1/F");
  trOUT->Branch("trackSipdSig_2", &trackSipdSig_2, "trackSipdSig_2/F");
  trOUT->Branch("trackSipdSig_3", &trackSipdSig_3, "trackSipdSig_3/F");
  trOUT->Branch("trackSipdSig_1_0", &trackSipdSig_1_0, "trackSipdSig_1_0/F");
  trOUT->Branch("trackSipdSig_1_1", &trackSipdSig_1_1, "trackSipdSig_1_1/F");
  trOUT->Branch("trackSipdSig_1_2", &trackSipdSig_1_2, "trackSipdSig_1_2/F");
  trOUT->Branch("trackSipdSig_1_3", &trackSipdSig_1_3, "trackSipdSig_1_3/F");
  trOUT->Branch("trackSipdSig_0_0", &trackSipdSig_0_0, "trackSipdSig_0_0/F");
  trOUT->Branch("trackSipdSig_0_1", &trackSipdSig_0_1, "trackSipdSig_0_1/F");
  trOUT->Branch("trackSipdSig_0_2", &trackSipdSig_0_2, "trackSipdSig_0_2/F");
  trOUT->Branch("trackSipdSig_0_3", &trackSipdSig_0_3, "trackSipdSig_0_3/F");


  trOUT->Branch("trackSip3dSig_0", &trackSip3dSig_0, "trackSip3dSig_0/F");
  trOUT->Branch("trackSip3dSig_1", &trackSip3dSig_1, "trackSip3dSig_1/F");
  trOUT->Branch("trackSip3dSig_2", &trackSip3dSig_2, "trackSip3dSig_2/F");
  trOUT->Branch("trackSip3dSig_3", &trackSip3dSig_3, "trackSip3dSig_3/F");
  trOUT->Branch("trackPtRel_0", &trackPtRel_0, "trackPtRel_0/F");
  trOUT->Branch("trackPtRel_1", &trackPtRel_1, "trackPtRel_1/F");
  trOUT->Branch("trackPtRel_2", &trackPtRel_2, "trackPtRel_2/F");
  trOUT->Branch("trackPtRel_3", &trackPtRel_3, "trackPtRel_3/F");
  trOUT->Branch("trackEtaRel_0", &trackEtaRel_0, "trackEtaRel_0/F");
  trOUT->Branch("trackEtaRel_1", &trackEtaRel_1, "trackEtaRel_1/F");
  trOUT->Branch("trackEtaRel_2", &trackEtaRel_2, "trackEtaRel_2/F");
  trOUT->Branch("trackEtaRel_3", &trackEtaRel_3, "trackEtaRel_3/F");

  trOUT->Branch("tau1_trackEtaRel_0", &tau1_trackEtaRel_0, "tau1_trackEtaRel_0/F");
  trOUT->Branch("tau1_trackEtaRel_1", &tau1_trackEtaRel_1, "tau1_trackEtaRel_1/F");
  trOUT->Branch("tau1_trackEtaRel_2", &tau1_trackEtaRel_2, "tau1_trackEtaRel_2/F");

  trOUT->Branch("tau0_trackEtaRel_0", &tau0_trackEtaRel_0, "tau0_trackEtaRel_0/F");
  trOUT->Branch("tau0_trackEtaRel_1", &tau0_trackEtaRel_1, "tau0_trackEtaRel_1/F");
  trOUT->Branch("tau0_trackEtaRel_2", &tau0_trackEtaRel_2, "tau0_trackEtaRel_2/F");


  trOUT->Branch("tau01_trackEtaRel_0", &tau01_trackEtaRel_0, "tau01_trackEtaRel_0/F");
  trOUT->Branch("tau01_trackEtaRel_1", &tau01_trackEtaRel_1, "tau01_trackEtaRel_1/F");
  trOUT->Branch("tau01_trackEtaRel_2", &tau01_trackEtaRel_2, "tau01_trackEtaRel_2/F");

  trOUT->Branch("tau00_trackEtaRel_0", &tau00_trackEtaRel_0, "tau00_trackEtaRel_0/F");
  trOUT->Branch("tau00_trackEtaRel_1", &tau00_trackEtaRel_1, "tau00_trackEtaRel_1/F");
  trOUT->Branch("tau00_trackEtaRel_2", &tau00_trackEtaRel_2, "tau00_trackEtaRel_2/F");

  trOUT->Branch("trackSip2dSigAboveCharm_0",   &trackSip2dSigAboveCharm_0, "trackSip2dSigAboveCharm_0/F");
  trOUT->Branch("trackSip2dSigAboveCharm_1",   &trackSip2dSigAboveCharm_1, "trackSip2dSigAboveCharm_1/F");
  trOUT->Branch("trackSip2dSigAboveBottom_0",   &trackSip2dSigAboveBottom_0, "trackSip2dSigAboveBottom_0/F");
  trOUT->Branch("trackSip2dSigAboveBottom_1",   &trackSip2dSigAboveBottom_1, "trackSip2dSigAboveBottom_1/F");

  trOUT->Branch("DoubleB", &DoubleB, "DoubleB/F");
  trOUT->Branch("BDTG",&BDTG,"BDTG/F");	

  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================

   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         PVez;
   Float_t         GenPVz;
   Float_t         nPUtrue;
   Int_t           nPU;
   Int_t           PU_bunch[68];   //[nPU]
   Float_t         PU_z[68];   //[nPU]
   Float_t         PU_sumpT_low[68];   //[nPU]
   Float_t         PU_sumpT_high[68];   //[nPU]
   Int_t           PU_ntrks_low[68];   //[nPU]
   Int_t           PU_ntrks_high[68];   //[nPU]  
   Int_t           nGenPruned;
   Float_t         GenPruned_pT[17];   //[nGenPruned]
   Float_t         GenPruned_eta[17];   //[nGenPruned]
   Float_t         GenPruned_phi[17];   //[nGenPruned]
   Float_t         GenPruned_mass[17];   //[nGenPruned]
   Int_t           GenPruned_pdgID[17];   //[nGenPruned]
   Int_t           GenPruned_status[17];   //[nGenPruned]
   Int_t           GenPruned_mother[17];   //[nGenPruned]
   Float_t         PV_x[56];   //[nPV]
   Float_t         PV_y[56];   //[nPV]
   Float_t         PV_z[56];   //[nPV]
   Float_t         PV_ex[56];   //[nPV]
   Float_t         PV_ey[56];   //[nPV]
   Float_t         PV_ez[56];   //[nPV]
   Float_t         PV_chi2[56];   //[nPV]
   Float_t         PV_ndf[56];   //[nPV]
   Int_t           PV_isgood[56];   //[nPV]
   Int_t           PV_isfake[56];   //[nPV]
   Int_t           FatJetInfo_nJet;
   Float_t         FatJetInfo_Jet_pt[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_genpt[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_residual[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_jes[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_eta[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_phi[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_mass[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_ntracks[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nseltracks[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_flavour[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nbHadrons[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_ncHadrons[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_ProbaN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_ProbaP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_Proba[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_BprobN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_BprobP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_Bprob[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SvxN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_Svx[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SvxNHP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SvxHP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombSvxN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombSvxP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombSvx[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombIVF[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombIVF_P[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_CombIVF_N[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftMuN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftMuP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftMu[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftElN[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftElP[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_SoftEl[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_DoubleSV[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_cMVA[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_SV_multi[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nSM[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nSE[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_looseID[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_tightID[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nFirstSV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastSV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nSV;
   Float_t         FatJetInfo_SV_x[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_y[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_z[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_ex[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_ey[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_ez[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_chi2[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_ndf[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_flight[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_flightErr[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_deltaR_jet[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_deltaR_sum_jet[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_deltaR_sum_dir[12];   //[FatJetInfo.nSV]
    Float_t         FatJetInfo_SV_vtx_pt[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_flight2D[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_flight2DErr[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_totCharge[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_vtxDistJetAxis[12];   //[FatJetInfo.nSV]
   Int_t           FatJetInfo_SV_nTrk[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_mass[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_vtx_eta[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_vtx_phi[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_EnergyRatio[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_dir_x[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_dir_y[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_SV_dir_z[12];   //[FatJetInfo.nSV]
   Float_t         FatJetInfo_Jet_ptSoftDrop[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_etaSoftDrop[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_phiSoftDrop[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_massSoftDrop[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_jecF0SoftDrop[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_ptPruned[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_etaPruned[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_phiPruned[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_massPruned[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_jecF0Pruned[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1IVF[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2IVF[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis1_px[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis1_py[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis1_pz[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis2_px[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis2_py[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tauAxis2_pz[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_z_ratio[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip3dSig_3[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip3dSig_2[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip3dSig_1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip3dSig_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackEtaRel_2[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackEtaRel_1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackEtaRel_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_2[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_2[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_nSecondaryVertices[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_nSecondaryVertices[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_flightDistance2dSig[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_flightDistance2dSig[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_vertexDeltaR[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_vertexDeltaR[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_vertexEnergyRatio[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_vertexEnergyRatio[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_vertexMass[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_vertexMass[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_vertexMass_corrected[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_vertexMass_corrected[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau1_vertexNTracks[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_tau2_vertexNTracks[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_BDTG_SV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nTrack;
   Float_t         FatJetInfo_Track_lengthTau[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_distTau[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Jet_trackSip2dSigAboveCharm_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip2dSigAboveCharm_1[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip2dSigAboveBottom_0[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Jet_trackSip2dSigAboveBottom_1[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nFirstTrack[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastTrack[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_Track_dxy[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_dz[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_length[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_dist[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IP2D[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IP2Dsig[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IP2Derr[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IP[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IPsig[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_IPerr[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_Proba[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_p[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_pt[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_eta[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_phi[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_chi2[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_charge[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitStrip[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitPixel[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitAll[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitTIB[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitTID[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitTOB[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitTEC[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitPXB[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_nHitPXF[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_isHitL1[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_PV[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_SV[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_PVweight[167];   //[FatJetInfo.nTrack]
   Float_t         FatJetInfo_Track_SVweight[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_isfromSV[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Track_isfromV0[167];   //[FatJetInfo.nTrack]
   Int_t           FatJetInfo_Jet_nFirstSE[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastSE[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nPFElectron;
   Int_t           FatJetInfo_PFElectron_IdxJet[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_pt[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_eta[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_phi[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_ptrel[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_deltaR[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_ratio[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_ratioRel[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_IP[8];   //[FatJetInfo.nPFElectron]
   Float_t         FatJetInfo_PFElectron_IP2D[8];   //[FatJetInfo.nPFElectron]
   Int_t           FatJetInfo_Jet_nFirstSM[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastSM[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nPFMuon;
   Int_t           FatJetInfo_PFMuon_IdxJet[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nMuHit[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nTkHit[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nPixHit[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nOutHit[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nTkLwM[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nPixLwM[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_nMatched[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_chi2[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_chi2Tk[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_PFMuon_isGlobal[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_pt[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_eta[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_phi[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_ptrel[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_deltaR[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_ratio[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_ratioRel[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_IP[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_IP2D[12];   //[FatJetInfo.nPFMuon]
   Float_t         FatJetInfo_PFMuon_dz[12];   //[FatJetInfo.nPFMuon]
   Int_t           FatJetInfo_Jet_nFirstTrkTagVarCSV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastTrkTagVarCSV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nFirstTrkEtaRelTagVarCSV[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastTrkEtaRelTagVarCSV[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackJetPt[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_jetNTracks[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_jetNTracksEtaRel[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSumJetEtRatio[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSumJetDeltaR[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSip2dValAboveCharm[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSip3dValAboveCharm[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_vertexCategory[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_jetNSecondaryVertices[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_vertexMass[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_vertexNTracks[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_vertexEnergyRatio[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_vertexJetDeltaR[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_flightDistance2dVal[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_flightDistance2dSig[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_flightDistance3dVal[12];   //[FatJetInfo.nJet]
   Float_t         FatJetInfo_TagVarCSV_flightDistance3dSig[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nTrkTagVarCSV;
   Int_t           FatJetInfo_nTrkEtaRelTagVarCSV;
   Float_t         FatJetInfo_TagVarCSV_trackMomentum[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackEta[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackPhi[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackPtRel[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackPPar[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackDeltaR[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackPtRatio[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackPParRatio[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackSip2dVal[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackSip2dSig[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackSip3dVal[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackSip3dSig[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackDecayLenVal[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackDecayLenSig[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackJetDistVal[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackJetDistSig[87];   //[FatJetInfo.nTrkTagVarCSV]
   Float_t         FatJetInfo_TagVarCSV_trackEtaRel[42];   //[FatJetInfo.nTrkEtaRelTagVarCSV]	
    Int_t           FatJetInfo_Jet_nSubJets_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nFirstSJ_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nLastSJ_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nsharedtracks_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nsubjettracks_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_Jet_nsharedsubjettracks_Pruned[12];   //[FatJetInfo.nJet]
   Int_t           FatJetInfo_nSubJet_Pruned;
   Int_t           FatJetInfo_SubJetIdx_Pruned[12];   //[FatJetInfo.nSubJet_Pruned]
   Int_t           PrunedSubJetInfo_nJet;
   Float_t         PrunedSubJetInfo_Jet_pt[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_genpt[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_residual[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_jes[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_eta[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_phi[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_mass[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_ntracks[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nseltracks[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_flavour[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nbHadrons[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_ncHadrons[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_ProbaN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_ProbaP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_Proba[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_BprobN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_BprobP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_Bprob[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SvxN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_Svx[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SvxNHP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SvxHP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombSvxN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombSvxP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombSvx[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombIVF[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombIVF_P[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_CombIVF_N[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftMuN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftMuP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftMu[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftElN[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftElP[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_SoftEl[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_DoubleSV[16];   //[PrunedSubJetInfo.nJet]
   Float_t         PrunedSubJetInfo_Jet_cMVA[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_SV_multi[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nSM[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nSE[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_looseID[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_tightID[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nFirstSV[16];   //[PrunedSubJetInfo.nJet]
   Int_t           PrunedSubJetInfo_Jet_nLastSV[16];   //[PrunedSubJetInfo.nJet]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVez;   //!
   TBranch        *b_GenPVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_PU_bunch;   //!
   TBranch        *b_PU_z;   //!
   TBranch        *b_PU_sumpT_low;   //! 
   TBranch        *b_PU_sumpT_high;   //!  
   TBranch        *b_PU_ntrks_low;   //!  
   TBranch        *b_PU_ntrks_high;   //!
   TBranch        *b_nGenPruned;   //!
   TBranch        *b_GenPruned_pT;   //!
   TBranch        *b_GenPruned_eta;   //!
   TBranch        *b_GenPruned_phi;   //!
   TBranch        *b_GenPruned_mass;   //!   
   TBranch        *b_GenPruned_pdgID;   //!
   TBranch        *b_GenPruned_status;   //!
   TBranch        *b_GenPruned_mother;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_ex;   //!
   TBranch        *b_PV_ey;   //!
   TBranch        *b_PV_ez;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ndf;   //!
   TBranch        *b_PV_isgood;   //!
   TBranch        *b_PV_isfake;   //!
   TBranch        *b_FatJetInfo_nJet;   //!
   TBranch        *b_FatJetInfo_Jet_pt;   //!
   TBranch        *b_FatJetInfo_Jet_genpt;   //!
   TBranch        *b_FatJetInfo_Jet_residual;   //!  
   TBranch        *b_FatJetInfo_Jet_jes;   //!  
   TBranch        *b_FatJetInfo_Jet_eta;   //!   
   TBranch        *b_FatJetInfo_Jet_phi;   //!   
   TBranch        *b_FatJetInfo_Jet_mass;   //!   
   TBranch        *b_FatJetInfo_Jet_ntracks;   //! 
   TBranch        *b_FatJetInfo_Jet_nseltracks;   //!
   TBranch        *b_FatJetInfo_Jet_flavour;   //!
   TBranch        *b_FatJetInfo_Jet_nbHadrons;   //!
   TBranch        *b_FatJetInfo_Jet_ncHadrons;   //!
   TBranch        *b_FatJetInfo_Jet_ProbaN;   //!
   TBranch        *b_FatJetInfo_Jet_ProbaP;   //!  
   TBranch        *b_FatJetInfo_Jet_Proba;   //!   
   TBranch        *b_FatJetInfo_Jet_BprobN;   //!  
   TBranch        *b_FatJetInfo_Jet_BprobP;   //!  
   TBranch        *b_FatJetInfo_Jet_Bprob;   //!
   TBranch        *b_FatJetInfo_Jet_SvxN;   //!
   TBranch        *b_FatJetInfo_Jet_Svx;   //!
   TBranch        *b_FatJetInfo_Jet_SvxNHP;   //!
   TBranch        *b_FatJetInfo_Jet_SvxHP;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvxN;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvxP;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvx;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF_P;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF_N;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMuN;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMuP;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMu;   //!
   TBranch        *b_FatJetInfo_Jet_SoftElN;   //!
   TBranch        *b_FatJetInfo_Jet_SoftElP;   //!
   TBranch        *b_FatJetInfo_Jet_SoftEl;   //!
   TBranch        *b_FatJetInfo_Jet_DoubleSV;   //!
   TBranch        *b_FatJetInfo_Jet_cMVA;   //!
   TBranch        *b_FatJetInfo_Jet_SV_multi;   //!
   TBranch        *b_FatJetInfo_Jet_nSM;   //!
   TBranch        *b_FatJetInfo_Jet_nSE;   //!
   TBranch        *b_FatJetInfo_Jet_looseID;   //!
   TBranch        *b_FatJetInfo_Jet_tightID;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSV;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSV;   //!
   TBranch        *b_FatJetInfo_nSV;   //!
   TBranch        *b_FatJetInfo_SV_x;   //!
   TBranch        *b_FatJetInfo_SV_y;   //!
   TBranch        *b_FatJetInfo_SV_z;   //!
   TBranch        *b_FatJetInfo_SV_ex;   //!
   TBranch        *b_FatJetInfo_SV_ey;   //!
   TBranch        *b_FatJetInfo_SV_ez;   //!
   TBranch        *b_FatJetInfo_SV_chi2;   //!
   TBranch        *b_FatJetInfo_SV_ndf;   //!
   TBranch        *b_FatJetInfo_SV_flight;   //!
   TBranch        *b_FatJetInfo_SV_flightErr;   //!
   TBranch        *b_FatJetInfo_SV_deltaR_jet;   //!
   TBranch        *b_FatJetInfo_SV_deltaR_sum_jet;   //!
   TBranch        *b_FatJetInfo_SV_deltaR_sum_dir;   //!
   TBranch        *b_FatJetInfo_SV_vtx_pt;   //!
   TBranch        *b_FatJetInfo_SV_flight2D;   //!
   TBranch        *b_FatJetInfo_SV_flight2DErr;   //!
   TBranch        *b_FatJetInfo_SV_totCharge;   //!
   TBranch        *b_FatJetInfo_SV_vtxDistJetAxis;   //!
   TBranch        *b_FatJetInfo_SV_nTrk;   //!
   TBranch        *b_FatJetInfo_SV_mass;   //!
   TBranch        *b_FatJetInfo_SV_vtx_eta;   //!
   TBranch        *b_FatJetInfo_SV_vtx_phi;   //!
   TBranch        *b_FatJetInfo_SV_EnergyRatio;   //!
   TBranch        *b_FatJetInfo_SV_dir_x;   //!
   TBranch        *b_FatJetInfo_SV_dir_y;   //!
   TBranch        *b_FatJetInfo_SV_dir_z;   //!
   TBranch        *b_FatJetInfo_Jet_ptSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_etaSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_phiSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_massSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_jecF0SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_ptPruned;   //!
   TBranch        *b_FatJetInfo_Jet_etaPruned;   //!
   TBranch        *b_FatJetInfo_Jet_phiPruned;   //!
   TBranch        *b_FatJetInfo_Jet_massPruned;   //!
   TBranch        *b_FatJetInfo_Jet_jecF0Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_tau1;   //!
   TBranch        *b_FatJetInfo_Jet_tau2;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_px;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_py;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_pz;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_px;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_py;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_pz;   //!
   TBranch        *b_FatJetInfo_Jet_z_ratio;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip3dSig_3;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip3dSig_2;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip3dSig_1;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip3dSig_0;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_nSecondaryVertices;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_nSecondaryVertices;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_flightDistance2dSig;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_flightDistance2dSig;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexDeltaR;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexDeltaR;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexEnergyRatio;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexEnergyRatio;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexMass;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexMass;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexMass_corrected;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexMass_corrected;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexNTracks;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexNTracks;   //!
   TBranch        *b_FatJetInfo_Jet_BDTG_SV;   //!
   TBranch        *b_FatJetInfo_nTrack;   //!
   TBranch        *b_FatJetInfo_Track_lengthTau;   //!
   TBranch        *b_FatJetInfo_Track_distTau;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip2dSigAboveCharm_0;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip2dSigAboveCharm_1;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip2dSigAboveBottom_0;   //!
   TBranch        *b_FatJetInfo_Jet_trackSip2dSigAboveBottom_1;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstTrack;   //!
   TBranch        *b_FatJetInfo_Jet_nLastTrack;   //!
   TBranch        *b_FatJetInfo_Track_dxy;   //!
   TBranch        *b_FatJetInfo_Track_dz;   //!
   TBranch        *b_FatJetInfo_Track_length;   //!
   TBranch        *b_FatJetInfo_Track_dist;   //!
   TBranch        *b_FatJetInfo_Track_IP2D;   //!
   TBranch        *b_FatJetInfo_Track_IP2Dsig;   //!
   TBranch        *b_FatJetInfo_Track_IP2Derr;   //!
   TBranch        *b_FatJetInfo_Track_IP;   //!
   TBranch        *b_FatJetInfo_Track_IPsig;   //!
   TBranch        *b_FatJetInfo_Track_IPerr;   //!
   TBranch        *b_FatJetInfo_Track_Proba;   //!
   TBranch        *b_FatJetInfo_Track_p;   //!
   TBranch        *b_FatJetInfo_Track_pt;   //!
   TBranch        *b_FatJetInfo_Track_eta;   //!
   TBranch        *b_FatJetInfo_Track_phi;   //!
   TBranch        *b_FatJetInfo_Track_chi2;   //!
   TBranch        *b_FatJetInfo_Track_charge;   //!
   TBranch        *b_FatJetInfo_Track_nHitStrip;   //!
   TBranch        *b_FatJetInfo_Track_nHitPixel;   //!
   TBranch        *b_FatJetInfo_Track_nHitAll;   //!
   TBranch        *b_FatJetInfo_Track_nHitTIB;   //!
   TBranch        *b_FatJetInfo_Track_nHitTID;   //!
   TBranch        *b_FatJetInfo_Track_nHitTOB;   //!
   TBranch        *b_FatJetInfo_Track_nHitTEC;   //!
   TBranch        *b_FatJetInfo_Track_nHitPXB;   //!
    TBranch        *b_FatJetInfo_Track_nHitPXF;   //!
   TBranch        *b_FatJetInfo_Track_isHitL1;   //!
   TBranch        *b_FatJetInfo_Track_PV;   //!
   TBranch        *b_FatJetInfo_Track_SV;   //!
   TBranch        *b_FatJetInfo_Track_PVweight;   //!
   TBranch        *b_FatJetInfo_Track_SVweight;   //!
   TBranch        *b_FatJetInfo_Track_isfromSV;   //!
   TBranch        *b_FatJetInfo_Track_isfromV0;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSE;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSE;   //!
   TBranch        *b_FatJetInfo_nPFElectron;   //!
   TBranch        *b_FatJetInfo_PFElectron_IdxJet;   //!
   TBranch        *b_FatJetInfo_PFElectron_pt;   //!
   TBranch        *b_FatJetInfo_PFElectron_eta;   //!
   TBranch        *b_FatJetInfo_PFElectron_phi;   //!
   TBranch        *b_FatJetInfo_PFElectron_ptrel;   //!
   TBranch        *b_FatJetInfo_PFElectron_deltaR;   //!
   TBranch        *b_FatJetInfo_PFElectron_ratio;   //!
   TBranch        *b_FatJetInfo_PFElectron_ratioRel;   //!
   TBranch        *b_FatJetInfo_PFElectron_IP;   //!
   TBranch        *b_FatJetInfo_PFElectron_IP2D;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSM;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSM;   //!
   TBranch        *b_FatJetInfo_nPFMuon;   //!
   TBranch        *b_FatJetInfo_PFMuon_IdxJet;   //!
   TBranch        *b_FatJetInfo_PFMuon_nMuHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nTkHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nPixHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nOutHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nTkLwM;   //!
   TBranch        *b_FatJetInfo_PFMuon_nPixLwM;   //!
   TBranch        *b_FatJetInfo_PFMuon_nMatched;   //!
   TBranch        *b_FatJetInfo_PFMuon_chi2;   //!
   TBranch        *b_FatJetInfo_PFMuon_chi2Tk;   //!
   TBranch        *b_FatJetInfo_PFMuon_isGlobal;   //!
   TBranch        *b_FatJetInfo_PFMuon_pt;   //!
   TBranch        *b_FatJetInfo_PFMuon_eta;   //!
   TBranch        *b_FatJetInfo_PFMuon_phi;   //!
   TBranch        *b_FatJetInfo_PFMuon_ptrel;   //!
   TBranch        *b_FatJetInfo_PFMuon_deltaR;   //!
   TBranch        *b_FatJetInfo_PFMuon_ratio;   //!
   TBranch        *b_FatJetInfo_PFMuon_ratioRel;   //!
   TBranch        *b_FatJetInfo_PFMuon_IP;   //!
   TBranch        *b_FatJetInfo_PFMuon_IP2D;   //!
   TBranch        *b_FatJetInfo_PFMuon_dz;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstTrkTagVarCSV;   //!
   TBranch        *b_FatJetInfo_Jet_nLastTrkTagVarCSV;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstTrkEtaRelTagVarCSV;   //!
   TBranch        *b_FatJetInfo_Jet_nLastTrkEtaRelTagVarCSV;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackJetPt;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_jetNTracks;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_jetNTracksEtaRel;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSumJetEtRatio;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSumJetDeltaR;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip2dValAboveCharm;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip3dValAboveCharm;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_vertexCategory;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_jetNSecondaryVertices;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_vertexMass;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_vertexNTracks;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_vertexEnergyRatio;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_vertexJetDeltaR;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_flightDistance2dVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_flightDistance2dSig;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_flightDistance3dVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_flightDistance3dSig;   //!
   TBranch        *b_FatJetInfo_nTrkTagVarCSV;   //!
   TBranch        *b_FatJetInfo_nTrkEtaRelTagVarCSV;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackMomentum;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackEta;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackPhi;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackPtRel;   //!
    TBranch        *b_FatJetInfo_TagVarCSV_trackPPar;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackDeltaR;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackPtRatio;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackPParRatio;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip2dVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip2dSig;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip3dVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackSip3dSig;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackDecayLenVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackDecayLenSig;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackJetDistVal;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackJetDistSig;   //!
   TBranch        *b_FatJetInfo_TagVarCSV_trackEtaRel;   //!
     TBranch        *b_FatJetInfo_Jet_nSubJets_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSJ_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSJ_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedtracks_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsubjettracks_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedsubjettracks_Pruned;   //!
   TBranch        *b_FatJetInfo_nSubJet_Pruned;   //!
   TBranch        *b_FatJetInfo_SubJetIdx_Pruned;   //!
   TBranch        *b_PrunedSubJetInfo_nJet;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_pt;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_genpt;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_residual;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_jes;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_eta;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_phi;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_mass;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ntracks;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nseltracks;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_flavour;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nbHadrons;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ncHadrons;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ProbaN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ProbaP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Proba;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_BprobN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_BprobP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Bprob;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Svx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxNHP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxHP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvxN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvxP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF_P;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF_N;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMuN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMuP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMu;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftElN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftElP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftEl;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_DoubleSV;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_cMVA;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SV_multi;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nSM;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nSE;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_looseID;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_tightID;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nFirstSV;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nLastSV;   //!
	
	
  std::cout << "Processing " << infilename << "..." << std::endl;    
  TFile *infile    = TFile::Open(infilename.c_str()); assert(infile);
  //  TTree *fChain = (TTree*)infile->Get("btaganaSubJets");   assert(fChain);
  TDirectory *d = (TDirectory*)infile->Get("btaganaFatJets"); assert(d);
  //  d->cd();
  TTree *fChain = (TTree*)d->Get("ttree");   assert(fChain);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVez", &PVez, &b_PVez);
   fChain->SetBranchAddress("GenPVz", &GenPVz, &b_GenPVz);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("PU_bunch", PU_bunch, &b_PU_bunch);
   fChain->SetBranchAddress("PU_z", PU_z, &b_PU_z);
   fChain->SetBranchAddress("PU_sumpT_low", PU_sumpT_low, &b_PU_sumpT_low);
   fChain->SetBranchAddress("PU_sumpT_high", PU_sumpT_high, &b_PU_sumpT_high);
   fChain->SetBranchAddress("PU_ntrks_low", PU_ntrks_low, &b_PU_ntrks_low);
   fChain->SetBranchAddress("PU_ntrks_high", PU_ntrks_high, &b_PU_ntrks_high);
   fChain->SetBranchAddress("nGenPruned", &nGenPruned, &b_nGenPruned);
   fChain->SetBranchAddress("GenPruned_pT", GenPruned_pT, &b_GenPruned_pT);
   fChain->SetBranchAddress("GenPruned_eta", GenPruned_eta, &b_GenPruned_eta);
   fChain->SetBranchAddress("GenPruned_phi", GenPruned_phi, &b_GenPruned_phi);
   fChain->SetBranchAddress("GenPruned_mass", GenPruned_mass, &b_GenPruned_mass);
   fChain->SetBranchAddress("GenPruned_pdgID", GenPruned_pdgID, &b_GenPruned_pdgID);
   fChain->SetBranchAddress("GenPruned_status", GenPruned_status, &b_GenPruned_status);
   fChain->SetBranchAddress("GenPruned_mother", GenPruned_mother, &b_GenPruned_mother);
   fChain->SetBranchAddress("PV_x", PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_ex", PV_ex, &b_PV_ex);
   fChain->SetBranchAddress("PV_ey", PV_ey, &b_PV_ey);
   fChain->SetBranchAddress("PV_ez", PV_ez, &b_PV_ez);
   fChain->SetBranchAddress("PV_chi2", PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ndf", PV_ndf, &b_PV_ndf);
   fChain->SetBranchAddress("PV_isgood", PV_isgood, &b_PV_isgood);
   fChain->SetBranchAddress("PV_isfake", PV_isfake, &b_PV_isfake);
   fChain->SetBranchAddress("FatJetInfo.nJet", &FatJetInfo_nJet, &b_FatJetInfo_nJet);
   fChain->SetBranchAddress("FatJetInfo.Jet_pt", FatJetInfo_Jet_pt, &b_FatJetInfo_Jet_pt);
   fChain->SetBranchAddress("FatJetInfo.Jet_genpt", FatJetInfo_Jet_genpt, &b_FatJetInfo_Jet_genpt);
   fChain->SetBranchAddress("FatJetInfo.Jet_residual", FatJetInfo_Jet_residual, &b_FatJetInfo_Jet_residual);
   fChain->SetBranchAddress("FatJetInfo.Jet_jes", FatJetInfo_Jet_jes, &b_FatJetInfo_Jet_jes);
   fChain->SetBranchAddress("FatJetInfo.Jet_eta", FatJetInfo_Jet_eta, &b_FatJetInfo_Jet_eta);
   fChain->SetBranchAddress("FatJetInfo.Jet_phi", FatJetInfo_Jet_phi, &b_FatJetInfo_Jet_phi);
   fChain->SetBranchAddress("FatJetInfo.Jet_mass", FatJetInfo_Jet_mass, &b_FatJetInfo_Jet_mass);
   fChain->SetBranchAddress("FatJetInfo.Jet_ntracks", FatJetInfo_Jet_ntracks, &b_FatJetInfo_Jet_ntracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_nseltracks", FatJetInfo_Jet_nseltracks, &b_FatJetInfo_Jet_nseltracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_flavour", FatJetInfo_Jet_flavour, &b_FatJetInfo_Jet_flavour);
   fChain->SetBranchAddress("FatJetInfo.Jet_nbHadrons", FatJetInfo_Jet_nbHadrons, &b_FatJetInfo_Jet_nbHadrons);
   fChain->SetBranchAddress("FatJetInfo.Jet_ncHadrons", FatJetInfo_Jet_ncHadrons, &b_FatJetInfo_Jet_ncHadrons);
   fChain->SetBranchAddress("FatJetInfo.Jet_ProbaN", FatJetInfo_Jet_ProbaN, &b_FatJetInfo_Jet_ProbaN);
   fChain->SetBranchAddress("FatJetInfo.Jet_ProbaP", FatJetInfo_Jet_ProbaP, &b_FatJetInfo_Jet_ProbaP);
   fChain->SetBranchAddress("FatJetInfo.Jet_Proba", FatJetInfo_Jet_Proba, &b_FatJetInfo_Jet_Proba);
   fChain->SetBranchAddress("FatJetInfo.Jet_BprobN", FatJetInfo_Jet_BprobN, &b_FatJetInfo_Jet_BprobN);
   fChain->SetBranchAddress("FatJetInfo.Jet_BprobP", FatJetInfo_Jet_BprobP, &b_FatJetInfo_Jet_BprobP);
   fChain->SetBranchAddress("FatJetInfo.Jet_Bprob", FatJetInfo_Jet_Bprob, &b_FatJetInfo_Jet_Bprob);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxN", FatJetInfo_Jet_SvxN, &b_FatJetInfo_Jet_SvxN);
   fChain->SetBranchAddress("FatJetInfo.Jet_Svx", FatJetInfo_Jet_Svx, &b_FatJetInfo_Jet_Svx);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxNHP", FatJetInfo_Jet_SvxNHP, &b_FatJetInfo_Jet_SvxNHP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxHP", FatJetInfo_Jet_SvxHP, &b_FatJetInfo_Jet_SvxHP);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvxN", FatJetInfo_Jet_CombSvxN, &b_FatJetInfo_Jet_CombSvxN);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvxP", FatJetInfo_Jet_CombSvxP, &b_FatJetInfo_Jet_CombSvxP);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvx", FatJetInfo_Jet_CombSvx, &b_FatJetInfo_Jet_CombSvx);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF", FatJetInfo_Jet_CombIVF, &b_FatJetInfo_Jet_CombIVF);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF_P", FatJetInfo_Jet_CombIVF_P, &b_FatJetInfo_Jet_CombIVF_P);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF_N", FatJetInfo_Jet_CombIVF_N, &b_FatJetInfo_Jet_CombIVF_N);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMuN", FatJetInfo_Jet_SoftMuN, &b_FatJetInfo_Jet_SoftMuN);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMuP", FatJetInfo_Jet_SoftMuP, &b_FatJetInfo_Jet_SoftMuP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMu", FatJetInfo_Jet_SoftMu, &b_FatJetInfo_Jet_SoftMu);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftElN", FatJetInfo_Jet_SoftElN, &b_FatJetInfo_Jet_SoftElN);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftElP", FatJetInfo_Jet_SoftElP, &b_FatJetInfo_Jet_SoftElP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftEl", FatJetInfo_Jet_SoftEl, &b_FatJetInfo_Jet_SoftEl);
   fChain->SetBranchAddress("FatJetInfo.Jet_DoubleSV", FatJetInfo_Jet_DoubleSV, &b_FatJetInfo_Jet_DoubleSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_cMVA", FatJetInfo_Jet_cMVA, &b_FatJetInfo_Jet_cMVA);
   fChain->SetBranchAddress("FatJetInfo.Jet_SV_multi", FatJetInfo_Jet_SV_multi, &b_FatJetInfo_Jet_SV_multi);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSM", FatJetInfo_Jet_nSM, &b_FatJetInfo_Jet_nSM);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSE", FatJetInfo_Jet_nSE, &b_FatJetInfo_Jet_nSE);
   fChain->SetBranchAddress("FatJetInfo.Jet_looseID", FatJetInfo_Jet_looseID, &b_FatJetInfo_Jet_looseID);
   fChain->SetBranchAddress("FatJetInfo.Jet_tightID", FatJetInfo_Jet_tightID, &b_FatJetInfo_Jet_tightID);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSV", FatJetInfo_Jet_nFirstSV, &b_FatJetInfo_Jet_nFirstSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSV", FatJetInfo_Jet_nLastSV, &b_FatJetInfo_Jet_nLastSV);
   fChain->SetBranchAddress("FatJetInfo.nSV", &FatJetInfo_nSV, &b_FatJetInfo_nSV);
   fChain->SetBranchAddress("FatJetInfo.SV_x", FatJetInfo_SV_x, &b_FatJetInfo_SV_x);
   fChain->SetBranchAddress("FatJetInfo.SV_y", FatJetInfo_SV_y, &b_FatJetInfo_SV_y);
   fChain->SetBranchAddress("FatJetInfo.SV_z", FatJetInfo_SV_z, &b_FatJetInfo_SV_z);
   fChain->SetBranchAddress("FatJetInfo.SV_ex", FatJetInfo_SV_ex, &b_FatJetInfo_SV_ex);
   fChain->SetBranchAddress("FatJetInfo.SV_ey", FatJetInfo_SV_ey, &b_FatJetInfo_SV_ey);
   fChain->SetBranchAddress("FatJetInfo.SV_ez", FatJetInfo_SV_ez, &b_FatJetInfo_SV_ez);
   fChain->SetBranchAddress("FatJetInfo.SV_chi2", FatJetInfo_SV_chi2, &b_FatJetInfo_SV_chi2);
   fChain->SetBranchAddress("FatJetInfo.SV_ndf", FatJetInfo_SV_ndf, &b_FatJetInfo_SV_ndf);
   fChain->SetBranchAddress("FatJetInfo.SV_flight", FatJetInfo_SV_flight, &b_FatJetInfo_SV_flight);
   fChain->SetBranchAddress("FatJetInfo.SV_flightErr", FatJetInfo_SV_flightErr, &b_FatJetInfo_SV_flightErr);
   fChain->SetBranchAddress("FatJetInfo.SV_deltaR_jet", FatJetInfo_SV_deltaR_jet, &b_FatJetInfo_SV_deltaR_jet);
   fChain->SetBranchAddress("FatJetInfo.SV_deltaR_sum_jet", FatJetInfo_SV_deltaR_sum_jet, &b_FatJetInfo_SV_deltaR_sum_jet);
   fChain->SetBranchAddress("FatJetInfo.SV_deltaR_sum_dir", FatJetInfo_SV_deltaR_sum_dir, &b_FatJetInfo_SV_deltaR_sum_dir);
   fChain->SetBranchAddress("FatJetInfo.SV_vtx_pt", FatJetInfo_SV_vtx_pt, &b_FatJetInfo_SV_vtx_pt);
   fChain->SetBranchAddress("FatJetInfo.SV_flight2D", FatJetInfo_SV_flight2D, &b_FatJetInfo_SV_flight2D);
   fChain->SetBranchAddress("FatJetInfo.SV_flight2DErr", FatJetInfo_SV_flight2DErr, &b_FatJetInfo_SV_flight2DErr);
   fChain->SetBranchAddress("FatJetInfo.SV_totCharge", FatJetInfo_SV_totCharge, &b_FatJetInfo_SV_totCharge);
   fChain->SetBranchAddress("FatJetInfo.SV_vtxDistJetAxis", FatJetInfo_SV_vtxDistJetAxis, &b_FatJetInfo_SV_vtxDistJetAxis);
   fChain->SetBranchAddress("FatJetInfo.SV_nTrk", FatJetInfo_SV_nTrk, &b_FatJetInfo_SV_nTrk);
   fChain->SetBranchAddress("FatJetInfo.SV_mass", FatJetInfo_SV_mass, &b_FatJetInfo_SV_mass);
   fChain->SetBranchAddress("FatJetInfo.SV_vtx_eta", FatJetInfo_SV_vtx_eta, &b_FatJetInfo_SV_vtx_eta);
   fChain->SetBranchAddress("FatJetInfo.SV_vtx_phi", FatJetInfo_SV_vtx_phi, &b_FatJetInfo_SV_vtx_phi);
   fChain->SetBranchAddress("FatJetInfo.SV_EnergyRatio", FatJetInfo_SV_EnergyRatio, &b_FatJetInfo_SV_EnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.SV_dir_x", FatJetInfo_SV_dir_x, &b_FatJetInfo_SV_dir_x);
   fChain->SetBranchAddress("FatJetInfo.SV_dir_y", FatJetInfo_SV_dir_y, &b_FatJetInfo_SV_dir_y);
   fChain->SetBranchAddress("FatJetInfo.SV_dir_z", FatJetInfo_SV_dir_z, &b_FatJetInfo_SV_dir_z);
   fChain->SetBranchAddress("FatJetInfo.Jet_ptSoftDrop", FatJetInfo_Jet_ptSoftDrop, &b_FatJetInfo_Jet_ptSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_etaSoftDrop", FatJetInfo_Jet_etaSoftDrop, &b_FatJetInfo_Jet_etaSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_phiSoftDrop", FatJetInfo_Jet_phiSoftDrop, &b_FatJetInfo_Jet_phiSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_massSoftDrop", FatJetInfo_Jet_massSoftDrop, &b_FatJetInfo_Jet_massSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_jecF0SoftDrop", FatJetInfo_Jet_jecF0SoftDrop, &b_FatJetInfo_Jet_jecF0SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_ptPruned", FatJetInfo_Jet_ptPruned, &b_FatJetInfo_Jet_ptPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_etaPruned", FatJetInfo_Jet_etaPruned, &b_FatJetInfo_Jet_etaPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_phiPruned", FatJetInfo_Jet_phiPruned, &b_FatJetInfo_Jet_phiPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_massPruned", FatJetInfo_Jet_massPruned, &b_FatJetInfo_Jet_massPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_jecF0Pruned", FatJetInfo_Jet_jecF0Pruned, &b_FatJetInfo_Jet_jecF0Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1", FatJetInfo_Jet_tau1, &b_FatJetInfo_Jet_tau1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2", FatJetInfo_Jet_tau2, &b_FatJetInfo_Jet_tau2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_px", FatJetInfo_Jet_tauAxis1_px, &b_FatJetInfo_Jet_tauAxis1_px);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_py", FatJetInfo_Jet_tauAxis1_py, &b_FatJetInfo_Jet_tauAxis1_py);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_pz", FatJetInfo_Jet_tauAxis1_pz, &b_FatJetInfo_Jet_tauAxis1_pz);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_px", FatJetInfo_Jet_tauAxis2_px, &b_FatJetInfo_Jet_tauAxis2_px);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_py", FatJetInfo_Jet_tauAxis2_py, &b_FatJetInfo_Jet_tauAxis2_py);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_pz", FatJetInfo_Jet_tauAxis2_pz, &b_FatJetInfo_Jet_tauAxis2_pz);
   fChain->SetBranchAddress("FatJetInfo.Jet_z_ratio", FatJetInfo_Jet_z_ratio, &b_FatJetInfo_Jet_z_ratio);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip3dSig_3", FatJetInfo_Jet_trackSip3dSig_3, &b_FatJetInfo_Jet_trackSip3dSig_3);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip3dSig_2", FatJetInfo_Jet_trackSip3dSig_2, &b_FatJetInfo_Jet_trackSip3dSig_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip3dSig_1", FatJetInfo_Jet_trackSip3dSig_1, &b_FatJetInfo_Jet_trackSip3dSig_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip3dSig_0", FatJetInfo_Jet_trackSip3dSig_0, &b_FatJetInfo_Jet_trackSip3dSig_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_2", FatJetInfo_Jet_trackEtaRel_2, &b_FatJetInfo_Jet_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_1", FatJetInfo_Jet_trackEtaRel_1, &b_FatJetInfo_Jet_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_0", FatJetInfo_Jet_trackEtaRel_0, &b_FatJetInfo_Jet_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_0", FatJetInfo_Jet_tau1_trackEtaRel_0, &b_FatJetInfo_Jet_tau1_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_1", FatJetInfo_Jet_tau1_trackEtaRel_1, &b_FatJetInfo_Jet_tau1_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_2", FatJetInfo_Jet_tau1_trackEtaRel_2, &b_FatJetInfo_Jet_tau1_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_0", FatJetInfo_Jet_tau2_trackEtaRel_0, &b_FatJetInfo_Jet_tau2_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_1", FatJetInfo_Jet_tau2_trackEtaRel_1, &b_FatJetInfo_Jet_tau2_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_2", FatJetInfo_Jet_tau2_trackEtaRel_2, &b_FatJetInfo_Jet_tau2_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_nSecondaryVertices", FatJetInfo_Jet_tau1_nSecondaryVertices, &b_FatJetInfo_Jet_tau1_nSecondaryVertices);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_nSecondaryVertices", FatJetInfo_Jet_tau2_nSecondaryVertices, &b_FatJetInfo_Jet_tau2_nSecondaryVertices);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_flightDistance2dSig", FatJetInfo_Jet_tau1_flightDistance2dSig, &b_FatJetInfo_Jet_tau1_flightDistance2dSig);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_flightDistance2dSig", FatJetInfo_Jet_tau2_flightDistance2dSig, &b_FatJetInfo_Jet_tau2_flightDistance2dSig);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexDeltaR", FatJetInfo_Jet_tau1_vertexDeltaR, &b_FatJetInfo_Jet_tau1_vertexDeltaR);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexDeltaR", FatJetInfo_Jet_tau2_vertexDeltaR, &b_FatJetInfo_Jet_tau2_vertexDeltaR);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexEnergyRatio", FatJetInfo_Jet_tau1_vertexEnergyRatio, &b_FatJetInfo_Jet_tau1_vertexEnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexEnergyRatio", FatJetInfo_Jet_tau2_vertexEnergyRatio, &b_FatJetInfo_Jet_tau2_vertexEnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexMass", FatJetInfo_Jet_tau1_vertexMass, &b_FatJetInfo_Jet_tau1_vertexMass);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexMass", FatJetInfo_Jet_tau2_vertexMass, &b_FatJetInfo_Jet_tau2_vertexMass);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexMass_corrected", FatJetInfo_Jet_tau1_vertexMass_corrected, &b_FatJetInfo_Jet_tau1_vertexMass_corrected);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexMass_corrected", FatJetInfo_Jet_tau2_vertexMass_corrected, &b_FatJetInfo_Jet_tau2_vertexMass_corrected);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexNTracks", FatJetInfo_Jet_tau1_vertexNTracks, &b_FatJetInfo_Jet_tau1_vertexNTracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexNTracks", FatJetInfo_Jet_tau2_vertexNTracks, &b_FatJetInfo_Jet_tau2_vertexNTracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_BDTG_SV", FatJetInfo_Jet_BDTG_SV, &b_FatJetInfo_Jet_BDTG_SV);
   fChain->SetBranchAddress("FatJetInfo.nTrack", &FatJetInfo_nTrack, &b_FatJetInfo_nTrack);
   fChain->SetBranchAddress("FatJetInfo.Track_lengthTau", FatJetInfo_Track_lengthTau, &b_FatJetInfo_Track_lengthTau);
   fChain->SetBranchAddress("FatJetInfo.Track_distTau", FatJetInfo_Track_distTau, &b_FatJetInfo_Track_distTau);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip2dSigAboveCharm_0", FatJetInfo_Jet_trackSip2dSigAboveCharm_0, &b_FatJetInfo_Jet_trackSip2dSigAboveCharm_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip2dSigAboveCharm_1", FatJetInfo_Jet_trackSip2dSigAboveCharm_1, &b_FatJetInfo_Jet_trackSip2dSigAboveCharm_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip2dSigAboveBottom_0", FatJetInfo_Jet_trackSip2dSigAboveBottom_0, &b_FatJetInfo_Jet_trackSip2dSigAboveBottom_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackSip2dSigAboveBottom_1", FatJetInfo_Jet_trackSip2dSigAboveBottom_1, &b_FatJetInfo_Jet_trackSip2dSigAboveBottom_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstTrack", FatJetInfo_Jet_nFirstTrack, &b_FatJetInfo_Jet_nFirstTrack);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastTrack", FatJetInfo_Jet_nLastTrack, &b_FatJetInfo_Jet_nLastTrack);
   fChain->SetBranchAddress("FatJetInfo.Track_dxy", FatJetInfo_Track_dxy, &b_FatJetInfo_Track_dxy);
   fChain->SetBranchAddress("FatJetInfo.Track_dz", FatJetInfo_Track_dz, &b_FatJetInfo_Track_dz);
   fChain->SetBranchAddress("FatJetInfo.Track_length", FatJetInfo_Track_length, &b_FatJetInfo_Track_length);
   fChain->SetBranchAddress("FatJetInfo.Track_dist", FatJetInfo_Track_dist, &b_FatJetInfo_Track_dist);
   fChain->SetBranchAddress("FatJetInfo.Track_IP2D", FatJetInfo_Track_IP2D, &b_FatJetInfo_Track_IP2D);
   fChain->SetBranchAddress("FatJetInfo.Track_IP2Dsig", FatJetInfo_Track_IP2Dsig, &b_FatJetInfo_Track_IP2Dsig);
   fChain->SetBranchAddress("FatJetInfo.Track_IP2Derr", FatJetInfo_Track_IP2Derr, &b_FatJetInfo_Track_IP2Derr);
   fChain->SetBranchAddress("FatJetInfo.Track_IP", FatJetInfo_Track_IP, &b_FatJetInfo_Track_IP);
   fChain->SetBranchAddress("FatJetInfo.Track_IPsig", FatJetInfo_Track_IPsig, &b_FatJetInfo_Track_IPsig);
   fChain->SetBranchAddress("FatJetInfo.Track_IPerr", FatJetInfo_Track_IPerr, &b_FatJetInfo_Track_IPerr);
   fChain->SetBranchAddress("FatJetInfo.Track_Proba", FatJetInfo_Track_Proba, &b_FatJetInfo_Track_Proba);
   fChain->SetBranchAddress("FatJetInfo.Track_p", FatJetInfo_Track_p, &b_FatJetInfo_Track_p);
   fChain->SetBranchAddress("FatJetInfo.Track_pt", FatJetInfo_Track_pt, &b_FatJetInfo_Track_pt);
   fChain->SetBranchAddress("FatJetInfo.Track_eta", FatJetInfo_Track_eta, &b_FatJetInfo_Track_eta);
   fChain->SetBranchAddress("FatJetInfo.Track_phi", FatJetInfo_Track_phi, &b_FatJetInfo_Track_phi);
   fChain->SetBranchAddress("FatJetInfo.Track_chi2", FatJetInfo_Track_chi2, &b_FatJetInfo_Track_chi2);
   fChain->SetBranchAddress("FatJetInfo.Track_charge", FatJetInfo_Track_charge, &b_FatJetInfo_Track_charge);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitStrip", FatJetInfo_Track_nHitStrip, &b_FatJetInfo_Track_nHitStrip);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitPixel", FatJetInfo_Track_nHitPixel, &b_FatJetInfo_Track_nHitPixel);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitAll", FatJetInfo_Track_nHitAll, &b_FatJetInfo_Track_nHitAll);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitTIB", FatJetInfo_Track_nHitTIB, &b_FatJetInfo_Track_nHitTIB);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitTID", FatJetInfo_Track_nHitTID, &b_FatJetInfo_Track_nHitTID);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitTOB", FatJetInfo_Track_nHitTOB, &b_FatJetInfo_Track_nHitTOB);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitTEC", FatJetInfo_Track_nHitTEC, &b_FatJetInfo_Track_nHitTEC);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitPXB", FatJetInfo_Track_nHitPXB, &b_FatJetInfo_Track_nHitPXB);
   fChain->SetBranchAddress("FatJetInfo.Track_nHitPXF", FatJetInfo_Track_nHitPXF, &b_FatJetInfo_Track_nHitPXF);
   fChain->SetBranchAddress("FatJetInfo.Track_isHitL1", FatJetInfo_Track_isHitL1, &b_FatJetInfo_Track_isHitL1);
    fChain->SetBranchAddress("FatJetInfo.Track_PV", FatJetInfo_Track_PV, &b_FatJetInfo_Track_PV);
   fChain->SetBranchAddress("FatJetInfo.Track_SV", FatJetInfo_Track_SV, &b_FatJetInfo_Track_SV);
   fChain->SetBranchAddress("FatJetInfo.Track_PVweight", FatJetInfo_Track_PVweight, &b_FatJetInfo_Track_PVweight);
   fChain->SetBranchAddress("FatJetInfo.Track_SVweight", FatJetInfo_Track_SVweight, &b_FatJetInfo_Track_SVweight);
   fChain->SetBranchAddress("FatJetInfo.Track_isfromSV", FatJetInfo_Track_isfromSV, &b_FatJetInfo_Track_isfromSV);
   fChain->SetBranchAddress("FatJetInfo.Track_isfromV0", FatJetInfo_Track_isfromV0, &b_FatJetInfo_Track_isfromV0);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSE", FatJetInfo_Jet_nFirstSE, &b_FatJetInfo_Jet_nFirstSE);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSE", FatJetInfo_Jet_nLastSE, &b_FatJetInfo_Jet_nLastSE);
   fChain->SetBranchAddress("FatJetInfo.nPFElectron", &FatJetInfo_nPFElectron, &b_FatJetInfo_nPFElectron);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IdxJet", FatJetInfo_PFElectron_IdxJet, &b_FatJetInfo_PFElectron_IdxJet);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_pt", FatJetInfo_PFElectron_pt, &b_FatJetInfo_PFElectron_pt);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_eta", FatJetInfo_PFElectron_eta, &b_FatJetInfo_PFElectron_eta);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_phi", FatJetInfo_PFElectron_phi, &b_FatJetInfo_PFElectron_phi);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ptrel", FatJetInfo_PFElectron_ptrel, &b_FatJetInfo_PFElectron_ptrel);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_deltaR", FatJetInfo_PFElectron_deltaR, &b_FatJetInfo_PFElectron_deltaR);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ratio", FatJetInfo_PFElectron_ratio, &b_FatJetInfo_PFElectron_ratio);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ratioRel", FatJetInfo_PFElectron_ratioRel, &b_FatJetInfo_PFElectron_ratioRel);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IP", FatJetInfo_PFElectron_IP, &b_FatJetInfo_PFElectron_IP);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IP2D", FatJetInfo_PFElectron_IP2D, &b_FatJetInfo_PFElectron_IP2D);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSM", FatJetInfo_Jet_nFirstSM, &b_FatJetInfo_Jet_nFirstSM);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSM", FatJetInfo_Jet_nLastSM, &b_FatJetInfo_Jet_nLastSM);
   fChain->SetBranchAddress("FatJetInfo.nPFMuon", &FatJetInfo_nPFMuon, &b_FatJetInfo_nPFMuon);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IdxJet", FatJetInfo_PFMuon_IdxJet, &b_FatJetInfo_PFMuon_IdxJet);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nMuHit", FatJetInfo_PFMuon_nMuHit, &b_FatJetInfo_PFMuon_nMuHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nTkHit", FatJetInfo_PFMuon_nTkHit, &b_FatJetInfo_PFMuon_nTkHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nPixHit", FatJetInfo_PFMuon_nPixHit, &b_FatJetInfo_PFMuon_nPixHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nOutHit", FatJetInfo_PFMuon_nOutHit, &b_FatJetInfo_PFMuon_nOutHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nTkLwM", FatJetInfo_PFMuon_nTkLwM, &b_FatJetInfo_PFMuon_nTkLwM);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nPixLwM", FatJetInfo_PFMuon_nPixLwM, &b_FatJetInfo_PFMuon_nPixLwM);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nMatched", FatJetInfo_PFMuon_nMatched, &b_FatJetInfo_PFMuon_nMatched);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_chi2", FatJetInfo_PFMuon_chi2, &b_FatJetInfo_PFMuon_chi2);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_chi2Tk", FatJetInfo_PFMuon_chi2Tk, &b_FatJetInfo_PFMuon_chi2Tk);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_isGlobal", FatJetInfo_PFMuon_isGlobal, &b_FatJetInfo_PFMuon_isGlobal);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_pt", FatJetInfo_PFMuon_pt, &b_FatJetInfo_PFMuon_pt);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_eta", FatJetInfo_PFMuon_eta, &b_FatJetInfo_PFMuon_eta);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_phi", FatJetInfo_PFMuon_phi, &b_FatJetInfo_PFMuon_phi);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ptrel", FatJetInfo_PFMuon_ptrel, &b_FatJetInfo_PFMuon_ptrel);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_deltaR", FatJetInfo_PFMuon_deltaR, &b_FatJetInfo_PFMuon_deltaR);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ratio", FatJetInfo_PFMuon_ratio, &b_FatJetInfo_PFMuon_ratio);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ratioRel", FatJetInfo_PFMuon_ratioRel, &b_FatJetInfo_PFMuon_ratioRel);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IP", FatJetInfo_PFMuon_IP, &b_FatJetInfo_PFMuon_IP);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IP2D", FatJetInfo_PFMuon_IP2D, &b_FatJetInfo_PFMuon_IP2D);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_dz", FatJetInfo_PFMuon_dz, &b_FatJetInfo_PFMuon_dz);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstTrkTagVarCSV", FatJetInfo_Jet_nFirstTrkTagVarCSV, &b_FatJetInfo_Jet_nFirstTrkTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastTrkTagVarCSV", FatJetInfo_Jet_nLastTrkTagVarCSV, &b_FatJetInfo_Jet_nLastTrkTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstTrkEtaRelTagVarCSV", FatJetInfo_Jet_nFirstTrkEtaRelTagVarCSV, &b_FatJetInfo_Jet_nFirstTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastTrkEtaRelTagVarCSV", FatJetInfo_Jet_nLastTrkEtaRelTagVarCSV, &b_FatJetInfo_Jet_nLastTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackJetPt", FatJetInfo_TagVarCSV_trackJetPt, &b_FatJetInfo_TagVarCSV_trackJetPt);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_jetNTracks", FatJetInfo_TagVarCSV_jetNTracks, &b_FatJetInfo_TagVarCSV_jetNTracks);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_jetNTracksEtaRel", FatJetInfo_TagVarCSV_jetNTracksEtaRel, &b_FatJetInfo_TagVarCSV_jetNTracksEtaRel);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSumJetEtRatio", FatJetInfo_TagVarCSV_trackSumJetEtRatio, &b_FatJetInfo_TagVarCSV_trackSumJetEtRatio);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSumJetDeltaR", FatJetInfo_TagVarCSV_trackSumJetDeltaR, &b_FatJetInfo_TagVarCSV_trackSumJetDeltaR);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip2dValAboveCharm", FatJetInfo_TagVarCSV_trackSip2dValAboveCharm, &b_FatJetInfo_TagVarCSV_trackSip2dValAboveCharm);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip2dSigAboveCharm", FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm, &b_FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip3dValAboveCharm", FatJetInfo_TagVarCSV_trackSip3dValAboveCharm, &b_FatJetInfo_TagVarCSV_trackSip3dValAboveCharm);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip3dSigAboveCharm", FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm, &b_FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_vertexCategory", FatJetInfo_TagVarCSV_vertexCategory, &b_FatJetInfo_TagVarCSV_vertexCategory);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_jetNSecondaryVertices", FatJetInfo_TagVarCSV_jetNSecondaryVertices, &b_FatJetInfo_TagVarCSV_jetNSecondaryVertices);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_vertexMass", FatJetInfo_TagVarCSV_vertexMass, &b_FatJetInfo_TagVarCSV_vertexMass);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_vertexNTracks", FatJetInfo_TagVarCSV_vertexNTracks, &b_FatJetInfo_TagVarCSV_vertexNTracks);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_vertexEnergyRatio", FatJetInfo_TagVarCSV_vertexEnergyRatio, &b_FatJetInfo_TagVarCSV_vertexEnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_vertexJetDeltaR", FatJetInfo_TagVarCSV_vertexJetDeltaR, &b_FatJetInfo_TagVarCSV_vertexJetDeltaR);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_flightDistance2dVal", FatJetInfo_TagVarCSV_flightDistance2dVal, &b_FatJetInfo_TagVarCSV_flightDistance2dVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_flightDistance2dSig", FatJetInfo_TagVarCSV_flightDistance2dSig, &b_FatJetInfo_TagVarCSV_flightDistance2dSig);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_flightDistance3dVal", FatJetInfo_TagVarCSV_flightDistance3dVal, &b_FatJetInfo_TagVarCSV_flightDistance3dVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_flightDistance3dSig", FatJetInfo_TagVarCSV_flightDistance3dSig, &b_FatJetInfo_TagVarCSV_flightDistance3dSig);
   fChain->SetBranchAddress("FatJetInfo.nTrkTagVarCSV", &FatJetInfo_nTrkTagVarCSV, &b_FatJetInfo_nTrkTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.nTrkEtaRelTagVarCSV", &FatJetInfo_nTrkEtaRelTagVarCSV, &b_FatJetInfo_nTrkEtaRelTagVarCSV);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackMomentum", FatJetInfo_TagVarCSV_trackMomentum, &b_FatJetInfo_TagVarCSV_trackMomentum);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackEta", FatJetInfo_TagVarCSV_trackEta, &b_FatJetInfo_TagVarCSV_trackEta);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackPhi", FatJetInfo_TagVarCSV_trackPhi, &b_FatJetInfo_TagVarCSV_trackPhi);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackPtRel", FatJetInfo_TagVarCSV_trackPtRel, &b_FatJetInfo_TagVarCSV_trackPtRel);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackPPar", FatJetInfo_TagVarCSV_trackPPar, &b_FatJetInfo_TagVarCSV_trackPPar);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackDeltaR", FatJetInfo_TagVarCSV_trackDeltaR, &b_FatJetInfo_TagVarCSV_trackDeltaR);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackPtRatio", FatJetInfo_TagVarCSV_trackPtRatio, &b_FatJetInfo_TagVarCSV_trackPtRatio);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackPParRatio", FatJetInfo_TagVarCSV_trackPParRatio, &b_FatJetInfo_TagVarCSV_trackPParRatio);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip2dVal", FatJetInfo_TagVarCSV_trackSip2dVal, &b_FatJetInfo_TagVarCSV_trackSip2dVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip2dSig", FatJetInfo_TagVarCSV_trackSip2dSig, &b_FatJetInfo_TagVarCSV_trackSip2dSig);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip3dVal", FatJetInfo_TagVarCSV_trackSip3dVal, &b_FatJetInfo_TagVarCSV_trackSip3dVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackSip3dSig", FatJetInfo_TagVarCSV_trackSip3dSig, &b_FatJetInfo_TagVarCSV_trackSip3dSig);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackDecayLenVal", FatJetInfo_TagVarCSV_trackDecayLenVal, &b_FatJetInfo_TagVarCSV_trackDecayLenVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackDecayLenSig", FatJetInfo_TagVarCSV_trackDecayLenSig, &b_FatJetInfo_TagVarCSV_trackDecayLenSig);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackJetDistVal", FatJetInfo_TagVarCSV_trackJetDistVal, &b_FatJetInfo_TagVarCSV_trackJetDistVal);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackJetDistSig", FatJetInfo_TagVarCSV_trackJetDistSig, &b_FatJetInfo_TagVarCSV_trackJetDistSig);
   fChain->SetBranchAddress("FatJetInfo.TagVarCSV_trackEtaRel", FatJetInfo_TagVarCSV_trackEtaRel, &b_FatJetInfo_TagVarCSV_trackEtaRel);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSubJets_Pruned", FatJetInfo_Jet_nSubJets_Pruned, &b_FatJetInfo_Jet_nSubJets_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSJ_Pruned", FatJetInfo_Jet_nFirstSJ_Pruned, &b_FatJetInfo_Jet_nFirstSJ_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSJ_Pruned", FatJetInfo_Jet_nLastSJ_Pruned, &b_FatJetInfo_Jet_nLastSJ_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedtracks_Pruned", FatJetInfo_Jet_nsharedtracks_Pruned, &b_FatJetInfo_Jet_nsharedtracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsubjettracks_Pruned", FatJetInfo_Jet_nsubjettracks_Pruned, &b_FatJetInfo_Jet_nsubjettracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedsubjettracks_Pruned", FatJetInfo_Jet_nsharedsubjettracks_Pruned, &b_FatJetInfo_Jet_nsharedsubjettracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.nSubJet_Pruned", &FatJetInfo_nSubJet_Pruned, &b_FatJetInfo_nSubJet_Pruned); 
   fChain->SetBranchAddress("FatJetInfo.SubJetIdx_Pruned", FatJetInfo_SubJetIdx_Pruned, &b_FatJetInfo_SubJetIdx_Pruned);
   fChain->SetBranchAddress("PrunedSubJetInfo.nJet", &PrunedSubJetInfo_nJet, &b_PrunedSubJetInfo_nJet);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_pt", PrunedSubJetInfo_Jet_pt, &b_PrunedSubJetInfo_Jet_pt);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_genpt", PrunedSubJetInfo_Jet_genpt, &b_PrunedSubJetInfo_Jet_genpt);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_residual", PrunedSubJetInfo_Jet_residual, &b_PrunedSubJetInfo_Jet_residual);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_jes", PrunedSubJetInfo_Jet_jes, &b_PrunedSubJetInfo_Jet_jes);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_eta", PrunedSubJetInfo_Jet_eta, &b_PrunedSubJetInfo_Jet_eta);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_phi", PrunedSubJetInfo_Jet_phi, &b_PrunedSubJetInfo_Jet_phi);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_mass", PrunedSubJetInfo_Jet_mass, &b_PrunedSubJetInfo_Jet_mass);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ntracks", PrunedSubJetInfo_Jet_ntracks, &b_PrunedSubJetInfo_Jet_ntracks);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nseltracks", PrunedSubJetInfo_Jet_nseltracks, &b_PrunedSubJetInfo_Jet_nseltracks);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_flavour", PrunedSubJetInfo_Jet_flavour, &b_PrunedSubJetInfo_Jet_flavour);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nbHadrons", PrunedSubJetInfo_Jet_nbHadrons, &b_PrunedSubJetInfo_Jet_nbHadrons);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ncHadrons", PrunedSubJetInfo_Jet_ncHadrons, &b_PrunedSubJetInfo_Jet_ncHadrons);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ProbaN", PrunedSubJetInfo_Jet_ProbaN, &b_PrunedSubJetInfo_Jet_ProbaN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ProbaP", PrunedSubJetInfo_Jet_ProbaP, &b_PrunedSubJetInfo_Jet_ProbaP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Proba", PrunedSubJetInfo_Jet_Proba, &b_PrunedSubJetInfo_Jet_Proba);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_BprobN", PrunedSubJetInfo_Jet_BprobN, &b_PrunedSubJetInfo_Jet_BprobN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_BprobP", PrunedSubJetInfo_Jet_BprobP, &b_PrunedSubJetInfo_Jet_BprobP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Bprob", PrunedSubJetInfo_Jet_Bprob, &b_PrunedSubJetInfo_Jet_Bprob);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxN", PrunedSubJetInfo_Jet_SvxN, &b_PrunedSubJetInfo_Jet_SvxN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Svx", PrunedSubJetInfo_Jet_Svx, &b_PrunedSubJetInfo_Jet_Svx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxNHP", PrunedSubJetInfo_Jet_SvxNHP, &b_PrunedSubJetInfo_Jet_SvxNHP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxHP", PrunedSubJetInfo_Jet_SvxHP, &b_PrunedSubJetInfo_Jet_SvxHP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvxN", PrunedSubJetInfo_Jet_CombSvxN, &b_PrunedSubJetInfo_Jet_CombSvxN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvxP", PrunedSubJetInfo_Jet_CombSvxP, &b_PrunedSubJetInfo_Jet_CombSvxP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvx", PrunedSubJetInfo_Jet_CombSvx, &b_PrunedSubJetInfo_Jet_CombSvx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF", PrunedSubJetInfo_Jet_CombIVF, &b_PrunedSubJetInfo_Jet_CombIVF);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF_P", PrunedSubJetInfo_Jet_CombIVF_P, &b_PrunedSubJetInfo_Jet_CombIVF_P);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF_N", PrunedSubJetInfo_Jet_CombIVF_N, &b_PrunedSubJetInfo_Jet_CombIVF_N);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMuN", PrunedSubJetInfo_Jet_SoftMuN, &b_PrunedSubJetInfo_Jet_SoftMuN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMuP", PrunedSubJetInfo_Jet_SoftMuP, &b_PrunedSubJetInfo_Jet_SoftMuP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMu", PrunedSubJetInfo_Jet_SoftMu, &b_PrunedSubJetInfo_Jet_SoftMu);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftElN", PrunedSubJetInfo_Jet_SoftElN, &b_PrunedSubJetInfo_Jet_SoftElN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftElP", PrunedSubJetInfo_Jet_SoftElP, &b_PrunedSubJetInfo_Jet_SoftElP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftEl", PrunedSubJetInfo_Jet_SoftEl, &b_PrunedSubJetInfo_Jet_SoftEl);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_DoubleSV", PrunedSubJetInfo_Jet_DoubleSV, &b_PrunedSubJetInfo_Jet_DoubleSV);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_cMVA", PrunedSubJetInfo_Jet_cMVA, &b_PrunedSubJetInfo_Jet_cMVA);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SV_multi", PrunedSubJetInfo_Jet_SV_multi, &b_PrunedSubJetInfo_Jet_SV_multi);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nSM", PrunedSubJetInfo_Jet_nSM, &b_PrunedSubJetInfo_Jet_nSM);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nSE", PrunedSubJetInfo_Jet_nSE, &b_PrunedSubJetInfo_Jet_nSE);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_looseID", PrunedSubJetInfo_Jet_looseID, &b_PrunedSubJetInfo_Jet_looseID);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_tightID", PrunedSubJetInfo_Jet_tightID, &b_PrunedSubJetInfo_Jet_tightID);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nFirstSV", PrunedSubJetInfo_Jet_nFirstSV, &b_PrunedSubJetInfo_Jet_nFirstSV);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nLastSV", PrunedSubJetInfo_Jet_nLastSV, &b_PrunedSubJetInfo_Jet_nLastSV);


   Long64_t nentries = fChain->GetEntriesFast();
   std::cout<<nentries<<std::endl;	
   int decade(0);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   double progress = 10.0*jentry/(1.0*nentries);
	   int k = TMath::FloorNint(progress);
	   if (k > decade)
		   cout<<10*k<<" %"<<endl;
	   decade = k;
	   Long64_t ientry = fChain->LoadTree(jentry);
    	   if (ientry < 0) break;
           nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   int countJ=0;

	   if (FatJetInfo_nJet>3) continue;
	   for(int i=0; i<FatJetInfo_nJet ; i++){

		   if(FatJetInfo_Jet_pt[i]>300. && FatJetInfo_Jet_massPruned[i]>70 && FatJetInfo_Jet_massPruned[i]< 200 && countJ<2){

			   SV_vtx_EnergyRatio_1=-3.;SV_vtx_EnergyRatio_0=-3.;SV_mass_0=-3.; z_ratio1=-3.; tau_dot=0.; nSV=-1.; nSL=-3.; z_ratio2=-3; z_ratio=-3;SV_mass_1=-3.;nSL_3=0;
			   PFLepton_IP2D=-3.;PFLepton_ratioRel=-3.;PFLepton_ptrel=-3.;PFLepton_deltaR=-3.;PFLepton_IP=-3.;
			   PFLepton_IP2D_0=-3.;PFLepton_ratioRel_0=-3.;PFLepton_ptrel_0=-3.;PFLepton_deltaR_0=-3.;
			   PFLepton_IP2D_1=-3.;PFLepton_ratioRel_1=-3.;PFLepton_ptrel_1=-3.;PFLepton_deltaR_1=-3.;
			   PFLepton_ratio=-3.;PFLepton_pt=-3.;

			   DoubleB = FatJetInfo_Jet_DoubleSV[i];
			   BDTG = FatJetInfo_Jet_BDTG_SV[i];
			   trackSip2dSig_0=-50;
			   trackSip2dSig_1=-50;
			   trackSip2dSig_2=-50;
			   trackSip2dSig_3=-50.;
			   trackSipdSig_0=-50;
			   trackSipdSig_1=-50;
			   trackSipdSig_2=-50;
			   trackSipdSig_3=-50.;
			   trackSip3dSig_0=-50;
			   trackSip3dSig_1=-50;
			   trackSip3dSig_2=-50;
			   trackSip3dSig_3=-50.;
			   trackSipdSig_1_0 =-50;
			   trackSipdSig_1_1 =-50;
			   trackSipdSig_1_2 =-50;
			   trackSipdSig_1_3 =-50;
			   trackSipdSig_0_0 =-50;
			   trackSipdSig_0_1 =-50;
			   trackSipdSig_0_2 =-50;
			   trackSipdSig_0_3 =-50;




			   //std::cout<<" *************************   NEW  EVENT  ****************************** "<<jentry<<std::endl;

			   tau_vertexMass_0=-3.;tau_vertexMass_1=-3.;tau_vertexMass_corrected_0=-3.;tau_vertexMass_corrected_1=-3.;
			   tau_vertexEnergyRatio_1=-3.;
			   tau_vertexEnergyRatio_0=-3.;tau_vertexDeltaR_1=-3.;     tau_vertexDeltaR_0=-3.; 
			   tau_vertexNtrk_0=-3;	tau_vertexNtrk_1=-3;
			   tau_vertexFlight2DErr_0=-3;tau_vertexFlight2DErr_1=-3; 
			   tau_flightDistance2dSig_0=-3;tau_flightDistance2dSig_1=-3;	
			   tau_vertexNsv_0=-3;       tau_vertexNsv_1=-3;	
			   tau1_trackEtaRel_0=-3;
			   tau1_trackEtaRel_1=-3;
			   tau1_trackEtaRel_2=-3;

			   tau0_trackEtaRel_0=-3;
			   tau0_trackEtaRel_1=-3;
			   tau0_trackEtaRel_2=-3;

			   trackSip2dSigAboveBottom_0=FatJetInfo_Jet_trackSip2dSigAboveBottom_0[i];
			   trackSip2dSigAboveBottom_1=FatJetInfo_Jet_trackSip2dSigAboveBottom_1[i]; 	
			   trackSip2dSigAboveCharm_0=FatJetInfo_Jet_trackSip2dSigAboveCharm_0[i];
			   trackSip2dSigAboveCharm_1=FatJetInfo_Jet_trackSip2dSigAboveCharm_1[i];


			   //Track_distance_TwoHighest3DSig=FatJetInfo_Jet_TrackDistance_TwoHighest3DSig[i];
			   if( (FatJetInfo_Jet_tau1_vertexMass[i]<0 || FatJetInfo_Jet_tau2_vertexMass[i]<0)){
				   tau_vertexMass_0 = TMath::Max(FatJetInfo_Jet_tau1_vertexMass[i],FatJetInfo_Jet_tau2_vertexMass[i]);
				   tau_vertexMass_corrected_0 = TMath::Max(FatJetInfo_Jet_tau1_vertexMass_corrected[i],FatJetInfo_Jet_tau2_vertexMass_corrected[i]);	
				   tau_vertexEnergyRatio_0 = TMath::Max(FatJetInfo_Jet_tau1_vertexEnergyRatio[i],FatJetInfo_Jet_tau2_vertexEnergyRatio[i]);	
				   tau_vertexDeltaR_0 = TMath::Max(FatJetInfo_Jet_tau1_vertexDeltaR[i],FatJetInfo_Jet_tau2_vertexDeltaR[i]);	
				   tau_flightDistance2dSig_0 = TMath::Max(FatJetInfo_Jet_tau1_flightDistance2dSig[i],FatJetInfo_Jet_tau2_flightDistance2dSig[i]);
				   tau_vertexNtrk_0= TMath::Max(FatJetInfo_Jet_tau1_vertexNTracks[i],FatJetInfo_Jet_tau2_vertexNTracks[i]);
				   tau_vertexNsv_0= TMath::Max(FatJetInfo_Jet_tau2_nSecondaryVertices[i],FatJetInfo_Jet_tau1_nSecondaryVertices[i]); 	
				   tau0_trackEtaRel_0 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_0[i],FatJetInfo_Jet_tau2_trackEtaRel_0[i]);
				   tau0_trackEtaRel_1 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_1[i],FatJetInfo_Jet_tau2_trackEtaRel_1[i]);
				   tau0_trackEtaRel_2 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_2[i],FatJetInfo_Jet_tau2_trackEtaRel_2[i]);


				   tau_vertexMass_1 = TMath::Min(FatJetInfo_Jet_tau1_vertexMass[i],FatJetInfo_Jet_tau2_vertexMass[i]);
				   tau_vertexMass_corrected_1 = TMath::Min(FatJetInfo_Jet_tau1_vertexMass_corrected[i],FatJetInfo_Jet_tau2_vertexMass_corrected[i]);
				   tau_vertexEnergyRatio_1 = TMath::Min(FatJetInfo_Jet_tau1_vertexEnergyRatio[i],FatJetInfo_Jet_tau2_vertexEnergyRatio[i]);
				   tau_vertexDeltaR_1 = TMath::Min(FatJetInfo_Jet_tau1_vertexDeltaR[i],FatJetInfo_Jet_tau2_vertexDeltaR[i]);
				   tau_flightDistance2dSig_1 = TMath::Min(FatJetInfo_Jet_tau1_flightDistance2dSig[i],FatJetInfo_Jet_tau2_flightDistance2dSig[i]);
				   tau_vertexNtrk_1= TMath::Min(FatJetInfo_Jet_tau1_vertexNTracks[i],FatJetInfo_Jet_tau2_vertexNTracks[i]);
				   tau_vertexNsv_1= TMath::Min(FatJetInfo_Jet_tau1_nSecondaryVertices[i],FatJetInfo_Jet_tau2_nSecondaryVertices[i]);	
				   tau1_trackEtaRel_0 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_0[i],FatJetInfo_Jet_tau2_trackEtaRel_0[i]);
				   tau1_trackEtaRel_1 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_1[i],FatJetInfo_Jet_tau2_trackEtaRel_1[i]);
				   tau1_trackEtaRel_2 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_2[i],FatJetInfo_Jet_tau2_trackEtaRel_2[i]);


			   }else if( (FatJetInfo_Jet_tau1_vertexMass[i]>0 && FatJetInfo_Jet_tau2_vertexMass[i]>0)){
				   tau_vertexMass_0 = FatJetInfo_Jet_tau1_vertexMass[i];
				   tau_vertexMass_corrected_0 =FatJetInfo_Jet_tau1_vertexMass_corrected[i];
				   tau_vertexEnergyRatio_0 = FatJetInfo_Jet_tau1_vertexEnergyRatio[i];
				   tau_vertexDeltaR_0 = FatJetInfo_Jet_tau1_vertexDeltaR[i];
				   tau_flightDistance2dSig_0 = FatJetInfo_Jet_tau1_flightDistance2dSig[i];
				   tau_vertexNtrk_0= FatJetInfo_Jet_tau1_vertexNTracks[i];
				   tau_vertexNsv_0= FatJetInfo_Jet_tau2_nSecondaryVertices[i];
				   tau0_trackEtaRel_0 = FatJetInfo_Jet_tau1_trackEtaRel_0[i];
				   tau0_trackEtaRel_1 = FatJetInfo_Jet_tau1_trackEtaRel_1[i];
				   tau0_trackEtaRel_2 = FatJetInfo_Jet_tau1_trackEtaRel_2[i];
				   tau_vertexMass_1 = FatJetInfo_Jet_tau2_vertexMass[i];
				   tau_vertexMass_corrected_1 = FatJetInfo_Jet_tau2_vertexMass_corrected[i];
				   tau_vertexEnergyRatio_1 = FatJetInfo_Jet_tau2_vertexEnergyRatio[i];
				   tau_vertexDeltaR_1 =FatJetInfo_Jet_tau2_vertexDeltaR[i];
				   tau_flightDistance2dSig_1 = FatJetInfo_Jet_tau2_flightDistance2dSig[i];
				   tau_vertexNtrk_1= FatJetInfo_Jet_tau2_vertexNTracks[i];
				   tau_vertexNsv_1=FatJetInfo_Jet_tau2_nSecondaryVertices[i];
				   tau1_trackEtaRel_0 = FatJetInfo_Jet_tau2_trackEtaRel_0[i];
				   tau1_trackEtaRel_1 =FatJetInfo_Jet_tau2_trackEtaRel_1[i];
				   tau1_trackEtaRel_2 = FatJetInfo_Jet_tau2_trackEtaRel_2[i];




			   }



			   if(tau_vertexEnergyRatio_0 >50.) tau_vertexEnergyRatio_0=50.;
			   if(tau_vertexEnergyRatio_1 >50.) tau_vertexEnergyRatio_1=50.;	
			   nsel++;
			   TVector3 tau_axis1, tau_axis2;

			   countJ++;	
			   int contSV=0;	 
			   std::map<double, int> VTXpt;
			   //for (mjj_map::reverse_iterator i=a.rbegin(); i!=a.rend(); ++i)
			   //qvol = FatJetInfo_Jet_qvol[i];
			   //if(FatJetInfo_Jet_SD_chi[i]>0.) SDchi = TMath::Log(FatJetInfo_Jet_SD_chi[i]);	
			   //std::cout<< qvol << "  qvol "<<SDchi <<std::endl;
			   ptPruned = FatJetInfo_Jet_ptPruned[i];
			   etaPruned = FatJetInfo_Jet_etaPruned[i];
			   phiPruned = FatJetInfo_Jet_phiPruned[i];
			   tau1 = FatJetInfo_Jet_tau1[i];
			   tau2 = FatJetInfo_Jet_tau2[i];


			   TagVarCSV1_jetNTracks              = FatJetInfo_TagVarCSV_jetNTracks[i];
			   TagVarCSV1_jetNTracksEtaRel        = FatJetInfo_TagVarCSV_jetNTracksEtaRel[i];
			   TagVarCSV1_trackSumJetEtRatio      = (FatJetInfo_TagVarCSV_trackSumJetEtRatio[i] < -800. ? -1. : FatJetInfo_TagVarCSV_trackSumJetEtRatio[i]);
			   TagVarCSV1_trackSumJetDeltaR       = (FatJetInfo_TagVarCSV_trackSumJetDeltaR[i] < -800. ? -1. : FatJetInfo_TagVarCSV_trackSumJetDeltaR[i]);
			   TagVarCSV1_trackSip2dValAboveCharm = (FatJetInfo_TagVarCSV_trackSip2dValAboveCharm[i] < -800. ?  -9. : FatJetInfo_TagVarCSV_trackSip2dValAboveCharm[i]);
			   TagVarCSV1_trackSip2dSigAboveCharm = (FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm[i]);
			   TagVarCSV1_trackSip3dValAboveCharm = (FatJetInfo_TagVarCSV_trackSip3dValAboveCharm[i] < -800. ?  -9. : FatJetInfo_TagVarCSV_trackSip3dValAboveCharm[i]);
			   TagVarCSV1_trackSip3dSigAboveCharm = (FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm[i]);
			   TagVarCSV1_vertexCategory          = (FatJetInfo_TagVarCSV_vertexCategory[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexCategory[i]);
			   TagVarCSV1_jetNSecondaryVertices   = FatJetInfo_TagVarCSV_jetNSecondaryVertices[i];
			   TagVarCSV1_vertexMass              = (FatJetInfo_TagVarCSV_vertexMass[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexMass[i]);
			   TagVarCSV1_vertexNTracks           = FatJetInfo_TagVarCSV_vertexNTracks[i];
			   TagVarCSV1_vertexEnergyRatio       = (FatJetInfo_TagVarCSV_vertexEnergyRatio[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexEnergyRatio[i]);
			   TagVarCSV1_vertexJetDeltaR         = (FatJetInfo_TagVarCSV_vertexJetDeltaR[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexJetDeltaR[i]);
			   TagVarCSV1_flightDistance2dVal     = (FatJetInfo_TagVarCSV_flightDistance2dVal[i] = -9999. ? -1. : FatJetInfo_TagVarCSV_flightDistance2dVal[i]);
			   TagVarCSV1_flightDistance2dSig     = (FatJetInfo_TagVarCSV_flightDistance2dSig[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance2dSig[i]);
			   TagVarCSV1_flightDistance3dVal     = (FatJetInfo_TagVarCSV_flightDistance3dVal[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance3dVal[i]);
			   TagVarCSV1_flightDistance3dSig     = (FatJetInfo_TagVarCSV_flightDistance3dSig[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance3dSig[i]);


			   tau1IVF = FatJetInfo_Jet_tau1IVF[i];
			   tau2IVF = FatJetInfo_Jet_tau2IVF[i];
			   massPruned = FatJetInfo_Jet_massPruned[i];		
			   nbHadrons = FatJetInfo_Jet_nbHadrons[i];
			   flavour = FatJetInfo_Jet_flavour[i];
			   tau_axis1.SetXYZ(FatJetInfo_Jet_tauAxis1_px[i],FatJetInfo_Jet_tauAxis1_py[i],FatJetInfo_Jet_tauAxis1_pz[i]);			
			   tau_axis2.SetXYZ(FatJetInfo_Jet_tauAxis2_px[i],FatJetInfo_Jet_tauAxis2_py[i],FatJetInfo_Jet_tauAxis2_pz[i]);

			   int contSL = 0 , contSE=0, contSM=0; float minPFMuon_ptrel=-999 ;
			   int contSL_1=0, contSL_2=0,  contSL_3=0;
			   for(int l=0; l<FatJetInfo_nPFElectron ;l++){
				   if(FatJetInfo_PFElectron_IdxJet[l] == i ) {
					   contSL++;
					   //	if(FatJetInfo_PFMuon_ptrel[l]>35.) contSL_1++;
					   //	else if(FatJetInfo_PFMuon_ptrel[l]>20.) contSL_2++;
					   //	else 
					   if(FatJetInfo_PFMuon_ptrel[l]>5.) contSL_3++; 

					   contSE++;
					   if(FatJetInfo_PFElectron_ptrel[l]> PFLepton_ptrel){
						   PFLepton_deltaR = FatJetInfo_PFElectron_deltaR[l];
						   PFLepton_ptrel = FatJetInfo_PFElectron_ptrel[l];
						   PFLepton_ratioRel = FatJetInfo_PFElectron_ratioRel[l];
						   PFLepton_IP2D = FatJetInfo_PFElectron_IP2D[l];
						   PFLepton_IP = FatJetInfo_PFElectron_IP[l];	
						   PFLepton_ratio = FatJetInfo_PFElectron_ratio[l];
						   PFLepton_pt = FatJetInfo_PFElectron_pt[l];	
					   }
					   if (FatJetInfo_PFElectron_ptrel[l]> PFLepton_ptrel_0) {

						   PFLepton_ptrel_1 = PFLepton_ptrel_0;
						   PFLepton_deltaR_1 = PFLepton_deltaR_0;
						   PFLepton_ratioRel_1 = PFLepton_ratioRel_0;
						   PFLepton_IP2D_1 = PFLepton_IP2D_0;



						   PFLepton_ptrel_0 = FatJetInfo_PFElectron_ptrel[l];
						   PFLepton_deltaR_0 = FatJetInfo_PFElectron_deltaR[l];
						   PFLepton_ratioRel_0 = FatJetInfo_PFElectron_ratioRel[l];
						   PFLepton_IP2D_0 = FatJetInfo_PFElectron_IP2D[l];

					   }
					   if(FatJetInfo_PFElectron_ptrel[l]>PFLepton_ptrel_1&& FatJetInfo_PFElectron_ptrel[l] < PFLepton_ptrel_0){

						   PFLepton_ptrel_1 = FatJetInfo_PFElectron_ptrel[l];
						   PFLepton_deltaR_1 = FatJetInfo_PFElectron_deltaR[l];
						   PFLepton_ratioRel_1 = FatJetInfo_PFElectron_ratioRel[l];
						   PFLepton_IP2D_1 = FatJetInfo_PFElectron_IP2D[l];
					   }

					   if(FatJetInfo_PFElectron_ptrel[l]> PFElectron_ptrel ){
						   PFElectron_deltaR = FatJetInfo_PFElectron_deltaR[l];
						   PFElectron_ptrel = FatJetInfo_PFElectron_ptrel[l];
						   PFElectron_ratioRel = FatJetInfo_PFElectron_ratioRel[l];
						   PFElectron_IP2D = FatJetInfo_PFElectron_IP2D[l];
					   }
				   }
			   }
			   for(int l=0; l<FatJetInfo_nPFMuon ;l++){
				   if(FatJetInfo_PFMuon_IdxJet[l] == i){
					   contSL++;
					   if(FatJetInfo_PFMuon_ptrel[l]>=5.) contSL_3++;	
					   contSM++;
					   if(FatJetInfo_PFMuon_ptrel[l]> PFLepton_ptrel ){
						   PFLepton_deltaR = FatJetInfo_PFMuon_deltaR[l];
						   PFLepton_ptrel = FatJetInfo_PFMuon_ptrel[l];
						   PFLepton_ratioRel = FatJetInfo_PFMuon_ratioRel[l];
						   PFLepton_IP2D = FatJetInfo_PFMuon_IP2D[l];
						   PFLepton_IP = FatJetInfo_PFMuon_IP[l];	
						   PFLepton_ratio = FatJetInfo_PFMuon_ratio[l];
						   PFLepton_pt = FatJetInfo_PFMuon_pt[l];      
					   }	

					   if (FatJetInfo_PFMuon_ptrel[l] > PFLepton_ptrel_0) {
						   PFLepton_ptrel_1 = PFLepton_ptrel_0;
						   PFLepton_deltaR_1 = PFLepton_deltaR_0;
						   PFLepton_ratioRel_1 = PFLepton_ratioRel_0;
						   PFLepton_IP2D_1 = PFLepton_IP2D_0;

						   PFLepton_ptrel_0 = FatJetInfo_PFMuon_ptrel[l];
						   PFLepton_deltaR_0 = FatJetInfo_PFMuon_deltaR[l];
						   PFLepton_ratioRel_0 = FatJetInfo_PFMuon_ratioRel[l];
						   PFLepton_IP2D_0 = FatJetInfo_PFMuon_IP2D[l];

					   }
					   if(FatJetInfo_PFMuon_ptrel[l]>PFLepton_ptrel_1&& FatJetInfo_PFMuon_ptrel[l] < PFLepton_ptrel_0){
						   PFLepton_ptrel_1 = FatJetInfo_PFMuon_ptrel[l];
						   PFLepton_deltaR_1 = FatJetInfo_PFMuon_deltaR[l];
						   PFLepton_ratioRel_1 = FatJetInfo_PFMuon_ratioRel[l];
						   PFLepton_IP2D_1 = FatJetInfo_PFMuon_IP2D[l];


					   }


					   if(FatJetInfo_PFMuon_ptrel[l]> minPFMuon_ptrel ){
						   PFMuon_deltaR = FatJetInfo_PFMuon_deltaR[l];
						   PFMuon_ptrel = FatJetInfo_PFMuon_ptrel[l];
						   PFMuon_ratioRel = FatJetInfo_PFMuon_ratioRel[l];
						   minPFMuon_ptrel = FatJetInfo_PFMuon_ratioRel[l];
						   PFMuon_IP2D = FatJetInfo_PFMuon_IP2D[l];
					   }
				   }

			   }
			   if(PFLepton_ptrel<5.) {PFLepton_ptrel =-3.;PFLepton_ratio=-3; PFLepton_pt=-3; PFLepton_IP2D =-3; PFLepton_deltaR=-3; PFLepton_IP = -3;}
			   nSM = contSM;
			   nSL= contSL;
			   nSL_1 =contSL_1;
			   nSL_2 =contSL_2;
			   nSL_3 =contSL_3;
			   TVector3 SV_axis0, SV_axis1;
			   TVector3 FJet;
			   FJet.SetPtEtaPhi(FatJetInfo_Jet_ptPruned[i], FatJetInfo_Jet_etaPruned[i], FatJetInfo_Jet_phiPruned[i]);//, FatJetInfo_Jet_massPruned[i]);
			   for(int j =0 ; j< FatJetInfo_nSV; j++)
			   {
				   TVector3 SV_;
				   SV_.SetXYZ(FatJetInfo_SV_dir_x[j], FatJetInfo_SV_dir_y[j], FatJetInfo_SV_dir_z[j]);//, FatJetInfo_SV_mass[j]);
				   if(SV_.DeltaR(FJet)<0.7) {
					   contSV++;
					   VTXpt[FatJetInfo_SV_flightErr[j]]=j;

				   }
			   }   
			   nSV= contSV;
			   if(nSV==1) cc1++;
			   if(nSV==2) cc2++; 
			   if(nSV==0) cc0++;
			   if(nSV==3) cc3++;
			   if(nSV>3) cc4++;
			   float dummyVtx=-3.;

			   if(contSV==0){
				   SV_vtx_EnergyRatio_pseudo_1=-3;SV_vtx_EnergyRatio_tighter_1=-3;SV_vtx_EnergyRatio_pseudo_0=-3;SV_vtx_EnergyRatio_tighter_0=-3;
				   SV_mass_correctedT_0=-3;SV_mass_corrected_0=3;SV_mass_correctedT_1=-3;SV_mass_corrected_1=-3;

				   SV_flight_0 = dummyVtx;
				   SV_flightErr_0 = dummyVtx;
				   SV_vtx_pt_0 = dummyVtx;
				   SV_flight2D_0 =  dummyVtx;
				   SV_flight2DErr_0 =   dummyVtx;
				   SV_totCharge_0 =    dummyVtx;
				   SV_nTrk_0 =dummyVtx;
				   SV_mass_0 =   dummyVtx;
				   SV_vtx_EnergyRatio_1 = dummyVtx;
				   SV_vtx_EnergyRatio_0 = dummyVtx;
				   SV_vtx_deltaR_1 = dummyVtx;	
				   SV_vtx_deltaR_0= dummyVtx;
				   SV_dirX_0 = dummyVtx;
				   SV_dirY_0 = dummyVtx;
				   SV_dirZ_0 = dummyVtx;

				   SV_dirX_1 = dummyVtx;
				   SV_dirY_1 = dummyVtx;
				   SV_dirZ_1 = dummyVtx;


				   SV_flight_1 = dummyVtx;
				   SV_flightErr_1 = dummyVtx;
				   SV_vtx_pt_1 = dummyVtx;
				   SV_flight2D_1 =  dummyVtx;
				   SV_flight2DErr_1 =   dummyVtx;
				   SV_totCharge_1 =    dummyVtx;
				   SV_nTrk_1 = dummyVtx;
				   SV_mass_1 =   dummyVtx;
				   SV_deltaR_axis_1 =  dummyVtx;
				   SV_deltaR_axis_0 =  dummyVtx;
				   SV_axis_deltaR = dummyVtx;
				   SV_axis_deltaPhi = dummyVtx; 
				   SV_axis_dot = dummyVtx;  //SV_axis1.Dot(SV_axis0);
				   SV_axis_cross =dummyVtx;  // (SV_axis1.Cross(SV_axis0)).Perp();
				   z_ratio = dummyVtx;  // SV_vtx_pt_1/(SV_vtx_pt_1+SV_vtx_pt_0);
				   z_ratio1 = dummyVtx; 
				   z_ratio2 = dummyVtx;
				   deltaR_tau2SV =dummyVtx;
				   deltaR_tau1SV =dummyVtx;
				   tau_dot = dummyVtx;


			   }


			   //	std::sort( VTXpt.begin(),VTXpt.end(),std::greater<float>() );
			   //if(contSV>0) std::cout<< "  " << contSV  << "  "<<VTXpt.size() << "  "<< iVtx->second << "  "<<std::endl;
			   int cont=0;     TLorentzVector SVsys, SV1, SV2;
			   for ( std::map<double, int>::iterator iVtx=VTXpt.begin(); iVtx!=VTXpt.end(); ++iVtx)        {
				   //for ( std::map<double, int>::reverse_iterator iVtx=VTXpt.rbegin(); iVtx!=VTXpt.rend(); ++iVtx)	{
				   //std::cout<<cont<<"  "<<FatJetInfo_SV_mass[iVtx->second] <<std::endl;	
				   cont++;	
				   if(contSV==1){


					   //SV_vtx_EnergyRatio_pseudo_0=FatJetInfo_SV_EnergyRatio_pseudo[iVtx->second];SV_vtx_EnergyRatio_tighter_0=FatJetInfo_SV_EnergyRatio_tighter[iVtx->second];
					   //SV_mass_correctedT_0=FatJetInfo_SV_mass_corrected_tighter[iVtx->second];SV_mass_corrected_0=FatJetInfo_SV_mass_corrected[iVtx->second];
					   //SV_vtx_EnergyRatio_pseudo_1=-3;SV_vtx_EnergyRatio_tighter_1=-3;
					   SV_mass_correctedT_1=-3;SV_mass_corrected_1=3;
					   SV_flight_0 = FatJetInfo_SV_flight[iVtx->second];
					   SV_flightErr_0 = FatJetInfo_SV_flightErr[iVtx->second];
					   SV_vtx_pt_0 = FatJetInfo_SV_vtx_pt[iVtx->second];
					   SV_flight2D_0 =  FatJetInfo_SV_flight2D[iVtx->second];
					   SV_flight2DErr_0 =   FatJetInfo_SV_flight2DErr[iVtx->second];
					   SV_totCharge_0 =    FatJetInfo_SV_totCharge[iVtx->second];
					   SV_nTrk_0 = FatJetInfo_SV_nTrk[iVtx->second];
					   SV_mass_0 =   FatJetInfo_SV_mass[iVtx->second];
					   SV_dirX_0 = FatJetInfo_SV_dir_x[iVtx->second];
					   SV_dirY_0 =  FatJetInfo_SV_dir_y[iVtx->second];
					   SV_dirZ_0 = FatJetInfo_SV_dir_z[iVtx->second];
					   SV_vtx_EnergyRatio_1 = dummyVtx;
					   SV_vtx_EnergyRatio_0 =  TMath::Min(100.,(double)FatJetInfo_SV_EnergyRatio[iVtx->second]);
					   SV_vtx_deltaR_0 =FatJetInfo_SV_deltaR_jet[iVtx->second];
					   SV_vtx_deltaR_1= dummyVtx;
					   SV_flight_1 = dummyVtx;
					   SV_flightErr_1 = dummyVtx;
					   SV_vtx_pt_1 = dummyVtx;
					   SV_flight2D_1 =  dummyVtx;
					   SV_flight2DErr_1 =   dummyVtx;
					   SV_totCharge_1 =    dummyVtx;
					   SV_nTrk_1 = dummyVtx;
					   SV_mass_1 =   dummyVtx;
					   SV_dirX_1 = dummyVtx;	
					   SV_dirY_1 = dummyVtx;
					   SV_dirZ_1 = dummyVtx;	
					   SV_axis_deltaR = dummyVtx;
					   SV_axis_deltaPhi = dummyVtx;
					   SV_axis_dot =  dummyVtx;
					   SV_deltaR_axis_1 =  dummyVtx;
					   z_ratio = tau_axis2.DeltaR(tau_axis1)*SV_vtx_pt_0/SV_mass_0;
					   z_ratio1 = tau_axis2.DeltaR(tau_axis1)*SV_vtx_pt_0/SV_mass_0;	
					   SV_axis0.SetXYZ(FatJetInfo_SV_dir_x[iVtx->second], FatJetInfo_SV_dir_y[iVtx->second], FatJetInfo_SV_dir_z[iVtx->second]);
					   TLorentzVector sub2;
					   sub2.SetPtEtaPhiM(SubJet2_pt,SubJet2_eta,SubJet2_phi,SubJet2_mass);
					   TLorentzVector sub1;
					   sub1.SetPtEtaPhiM(SubJet1_pt,SubJet1_eta,SubJet1_phi,SubJet1_mass);
					   TVector3 SJ_axis0 = sub1.Vect();
					   TVector3 SJ_axis1= sub2.Vect();
					   if(SV_axis0.DeltaR(tau_axis2)<SV_axis0.DeltaR(tau_axis1)) {

						   tau_dot = tau_axis2.Dot(SV_axis0)/(tau_axis2.Mag()*SV_axis0.Mag());
						   SV_deltaR_axis_0 = tau_axis2.DeltaR(SV_axis0);
						   //if(tau_axis1.DeltaR(SJ_axis0)<tau_axis1.DeltaR(SJ_axis1)) tau_dot = tau_dot+tau_axis1.Dot(SJ_axis0)/(tau_axis1.Mag()*SJ_axis0.Mag());
						   //else tau_dot = tau_dot+tau_axis1.Dot(SJ_axis1)/(tau_axis1.Mag()*SJ_axis1.Mag());

					   }else {
						   tau_dot = tau_axis1.Dot(SV_axis0)/(tau_axis1.Mag()*SV_axis0.Mag());
						   SV_deltaR_axis_0 = tau_axis1.DeltaR(SV_axis0);

						   //if(tau_axis2.DeltaR(SJ_axis0)<tau_axis2.DeltaR(SJ_axis1)) tau_dot = tau_dot+tau_axis2.Dot(SJ_axis0)/(tau_axis2.Mag()*SJ_axis0.Mag());
						   //else tau_dot = tau_dot+tau_axis2.Dot(SJ_axis1)/(tau_axis2.Mag()*SJ_axis1.Mag());
					   }	



				   }

				   if(contSV>1){

					   if(cont==1){
						   SV_flight_0 = FatJetInfo_SV_flight[iVtx->second];
						   SV_flightErr_0 = FatJetInfo_SV_flightErr[iVtx->second];
						   SV_vtx_pt_0 = FatJetInfo_SV_vtx_pt[iVtx->second];
						   SV_flight2D_0 =  FatJetInfo_SV_flight2D[iVtx->second];
						   SV_flight2DErr_0 =   FatJetInfo_SV_flight2DErr[iVtx->second];
						   SV_totCharge_0 =    FatJetInfo_SV_totCharge[iVtx->second];
						   SV_nTrk_0 = FatJetInfo_SV_nTrk[iVtx->second];
						   SV_mass_0 =   FatJetInfo_SV_mass[iVtx->second];
						   SV_dirX_0 = FatJetInfo_SV_dir_x[iVtx->second];
						   SV_dirY_0 =  FatJetInfo_SV_dir_y[iVtx->second];
						   SV_dirZ_0 = FatJetInfo_SV_dir_z[iVtx->second];	
						   //SV_vtx_EnergyRatio_pseudo_0=FatJetInfo_SV_EnergyRatio_pseudo[iVtx->second];SV_vtx_EnergyRatio_tighter_0=FatJetInfo_SV_EnergyRatio_tighter[iVtx->second];
						   //SV_mass_correctedT_0=FatJetInfo_SV_mass_corrected_tighter[iVtx->second];SV_mass_corrected_0=FatJetInfo_SV_mass_corrected[iVtx->second];  	
						   SV1.SetPtEtaPhiM(FatJetInfo_SV_vtx_pt[iVtx->second], FatJetInfo_SV_vtx_eta[iVtx->second], FatJetInfo_SV_vtx_phi[iVtx->second], FatJetInfo_SV_mass[iVtx->second]);


						   //std::cout<<FatJetInfo_SV_EnergyRatio[iVtx->second]<<std::endl;
						   SV_vtx_EnergyRatio_0 =  (double)TMath::Min(100.,(double)FatJetInfo_SV_EnergyRatio[iVtx->second]);
						   SV_vtx_deltaR_0 =FatJetInfo_SV_deltaR_jet[iVtx->second];
						   SV_axis0.SetXYZ(FatJetInfo_SV_dir_x[iVtx->second], FatJetInfo_SV_dir_y[iVtx->second], FatJetInfo_SV_dir_z[iVtx->second]);
						   //std::cout<<SV_axis0.Eta()<< " 1 "<< SV_axis0.Phi()<< ""<< SV_axis0.X() <<std::endl;
						   // std::cout<< FatJetInfo_SV_dir_x[iVtx->second] << " FatJetInfo_SV_dir_x[iVtx->second] "<<std::endl;	
					   }else if(cont==2){
						   SV_flight_1 = FatJetInfo_SV_flight[iVtx->second];
						   //SV_vtx_EnergyRatio_pseudo_1=FatJetInfo_SV_EnergyRatio_pseudo[iVtx->second];SV_vtx_EnergyRatio_tighter_1=FatJetInfo_SV_EnergyRatio_tighter[iVtx->second];
						   //SV_mass_correctedT_1=FatJetInfo_SV_mass_corrected_tighter[iVtx->second];SV_mass_corrected_1=FatJetInfo_SV_mass_corrected[iVtx->second];
						   SV_dirX_1 = FatJetInfo_SV_dir_x[iVtx->second];
						   SV_dirY_1 =  FatJetInfo_SV_dir_y[iVtx->second];
						   SV_dirZ_1 = FatJetInfo_SV_dir_z[iVtx->second];
						   SV2.SetPtEtaPhiM(FatJetInfo_SV_vtx_pt[iVtx->second], FatJetInfo_SV_vtx_eta[iVtx->second], FatJetInfo_SV_vtx_phi[iVtx->second], FatJetInfo_SV_mass[iVtx->second]);

						   SV_vtx_deltaR_1 =FatJetInfo_SV_deltaR_jet[iVtx->second];
						   SV_flightErr_1 = FatJetInfo_SV_flightErr[iVtx->second];
						   SV_vtx_pt_1 = FatJetInfo_SV_vtx_pt[iVtx->second];
						   SV_flight2D_1 =  FatJetInfo_SV_flight2D[iVtx->second];
						   SV_flight2DErr_1 =   FatJetInfo_SV_flight2DErr[iVtx->second];
						   SV_totCharge_1 =    FatJetInfo_SV_totCharge[iVtx->second];
						   SV_nTrk_1 = FatJetInfo_SV_nTrk[iVtx->second];
						   SV_vtx_EnergyRatio_1 = TMath::Min(100.,(double) FatJetInfo_SV_EnergyRatio[iVtx->second]);	
						   SV_mass_1 =   FatJetInfo_SV_mass[iVtx->second];
						   SV_axis1.SetXYZ(FatJetInfo_SV_dir_x[iVtx->second], FatJetInfo_SV_dir_y[iVtx->second], FatJetInfo_SV_dir_z[iVtx->second]);
						   SV_axis_deltaR = SV_axis1.DeltaR(SV_axis0);
						   SV_axis_deltaPhi = SV_axis1.DeltaPhi(SV_axis0);

						   SV_axis_dot = SV_axis1.Dot(SV_axis0);
						   SV_axis_cross = (SV_axis1.Cross(SV_axis0)).Mag();
						   SVsys = SV1+SV2;
						   z_ratio = SV_axis_deltaR*((SV_vtx_pt_1+SV_vtx_pt_0)/((SV2+SV1).M()));
						   z_ratio = SV_axis_deltaR*((SV_vtx_pt_1)/((SV2+SV1).M()));
						   if(contSV==2)z_ratio2 = SV_axis_deltaR*(SV_vtx_pt_0/(SV1+SV2).M());


						   deltaR_tau1SV = tau_axis2.DeltaR(SV_axis0);
						   deltaR_tau2SV = tau_axis1.DeltaR(SV_axis0);
						   if(deltaR_tau1SV<deltaR_tau2SV) {
							   tau_dot = tau_axis1.Dot(SV_axis0)/(tau_axis1.Mag()*SV_axis0.Mag());//+tau_axis1.Dot(SV_axis1)/(tau_axis1.Mag()*SV_axis1.Mag());
							   SV_deltaR_axis_0 = tau_axis1.DeltaR(SV_axis0);
						   }
						   else {
							   tau_dot = tau_axis2.Dot(SV_axis0)/(tau_axis2.Mag()*SV_axis0.Mag());//+tau_axis1.Dot(SV_axis0)/(tau_axis1.Mag()*SV_axis0.Mag()); 
							   SV_deltaR_axis_0 = tau_axis2.DeltaR(SV_axis0);
						   }
						   deltaR_tau1SV = tau_axis2.DeltaR(SV_axis1);
						   deltaR_tau2SV = tau_axis1.DeltaR(SV_axis1);
						   if(deltaR_tau1SV<deltaR_tau2SV) {
							   SV_deltaR_axis_1 = tau_axis1.DeltaR(SV_axis1);
						   }
						   else {
							   SV_deltaR_axis_1 = tau_axis2.DeltaR(SV_axis1);
						   }



					   }
				   }

			   }


			   h_multiplicity->Fill(countJ);

			   jetNTracks              = FatJetInfo_TagVarCSV_jetNTracks[i];
			   jetNTracks_ip              = FatJetInfo_Jet_ntracks[i];	
			   jetNTracksEtaRel        = FatJetInfo_TagVarCSV_jetNTracksEtaRel[i];


			   trackSumJetEtRatio      = (FatJetInfo_TagVarCSV_trackSumJetEtRatio[i] < -800. ? -1. : FatJetInfo_TagVarCSV_trackSumJetEtRatio[i]);
			   trackSumJetDeltaR       = (FatJetInfo_TagVarCSV_trackSumJetDeltaR[i] < -800. ? -1. : FatJetInfo_TagVarCSV_trackSumJetDeltaR[i]);
			   trackSip2dValAboveCharm = (FatJetInfo_TagVarCSV_trackSip2dValAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip2dValAboveCharm[i]);
			   trackSip2dSigAboveCharm = (FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip2dSigAboveCharm[i]);
			   trackSip3dValAboveCharm = (FatJetInfo_TagVarCSV_trackSip3dValAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip3dValAboveCharm[i]);
			   trackSip3dSigAboveCharm = (FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm[i] < -800. ? -99. : FatJetInfo_TagVarCSV_trackSip3dSigAboveCharm[i]);

			   vertexCategory          = (FatJetInfo_TagVarCSV_vertexCategory[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexCategory[i]);
			   jetNSecondaryVertices   = FatJetInfo_TagVarCSV_jetNSecondaryVertices[i];
			   vertexMass              = (FatJetInfo_TagVarCSV_vertexMass[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexMass[i]);
			   vertexNTracks           = FatJetInfo_TagVarCSV_vertexNTracks[i];
			   vertexEnergyRatio       = (FatJetInfo_TagVarCSV_vertexEnergyRatio[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexEnergyRatio[i]);
			   vertexJetDeltaR         = (FatJetInfo_TagVarCSV_vertexJetDeltaR[i] < -800. ? -1. : FatJetInfo_TagVarCSV_vertexJetDeltaR[i]);

			   flightDistance2dVal     = (FatJetInfo_TagVarCSV_flightDistance2dVal[i] = -9999. ? -1. : FatJetInfo_TagVarCSV_flightDistance2dVal[i]);
			   flightDistance2dSig     = (FatJetInfo_TagVarCSV_flightDistance2dSig[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance2dSig[i]);
			   flightDistance3dVal     = (FatJetInfo_TagVarCSV_flightDistance3dVal[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance3dVal[i]);
			   flightDistance3dSig     = (FatJetInfo_TagVarCSV_flightDistance3dSig[i] < -800. ? -1. : FatJetInfo_TagVarCSV_flightDistance3dSig[i]);


			   std::vector<float> IP2Ds, IP2Ds_1, IP2Ds_2;
			   std::vector<float> IP3Ds;
			   std::vector<float> IPD;	
			   std::vector<float> PtRel;
			   std::vector<float> etaRels_1, etaRels_2;
			   int contTrk=0;	
			   int contTrk_tight=0;	
			   for (int iTrk = 0; iTrk < jetNTracks;  ++iTrk){
				   IPD.push_back(FatJetInfo_Track_IPsig[iTrk]);
				   TVector3 track_;
				   track_.SetPtEtaPhi(FatJetInfo_Track_pt[iTrk],FatJetInfo_Track_eta[iTrk],FatJetInfo_Track_phi[iTrk]);
				   if(track_.DeltaR(tau_axis1)<track_.DeltaR(tau_axis2) ) {
					   etaRels_1.push_back(TMath::Abs(etaRel(tau_axis1,track_)));
				   }
				   else {
					   etaRels_2.push_back(TMath::Abs(etaRel(tau_axis2,track_)));
				   }	

				   if(FatJetInfo_Track_isfromV0[iTrk]!=1 && FatJetInfo_Track_lengthTau[iTrk]<5. && FatJetInfo_Track_distTau[iTrk]<0.07) {	
					   float temp = FatJetInfo_Track_IPsig[iTrk];	 
					   if(FatJetInfo_Track_IPsig[iTrk]<-50.) temp =-50.;	
					   IP2Ds.push_back(temp);	
					   if(track_.DeltaR(tau_axis1)<track_.DeltaR(tau_axis2) ) {
						   IP2Ds_1.push_back(temp);
					   }
					   else {	
						   IP2Ds_2.push_back(temp);
					   }

					   contTrk++;
				   }
				   if(FatJetInfo_Track_isfromV0[iTrk]!=1 && FatJetInfo_Track_lengthTau[iTrk]<5. && FatJetInfo_Track_distTau[iTrk]<0.055) {
					   IPD.push_back(FatJetInfo_Track_IPsig[iTrk]);
					   contTrk_tight++;
				   }	

			   }		
			   for (int iTrk = FatJetInfo_Jet_nFirstTrkTagVarCSV[i]; iTrk < FatJetInfo_Jet_nLastTrkTagVarCSV[i]; ++iTrk){
				   //IP2Ds.push_back( FatJetInfo_TagVarCSV_trackSip2dSig[iTrk] );
				   IP3Ds.push_back( FatJetInfo_TagVarCSV_trackSip3dSig[iTrk] );
				   PtRel.push_back( FatJetInfo_TagVarCSV_trackPtRel[iTrk] );
			   }

			   std::vector<float> etaRels; // stores |trackEtaRel|!
			   for (int iTrk = FatJetInfo_Jet_nFirstTrkEtaRelTagVarCSV[i]; iTrk < FatJetInfo_Jet_nLastTrkEtaRelTagVarCSV[i]; ++iTrk)
				   etaRels.push_back( fabs(FatJetInfo_TagVarCSV_trackEtaRel[iTrk]) );

			   std::sort( IP2Ds.begin(),IP2Ds.end(),std::greater<float>() );
			   std::sort( IP3Ds.begin(),IP3Ds.end(),std::greater<float>() );
			   std::sort( PtRel.begin(),PtRel.end(),std::greater<float>() );
			   std::sort(IPD.begin(),IPD.end(),std::greater<float>() );
			   std::sort( IP2Ds_1.begin(),IP2Ds_1.end(),std::greater<float>() );
			   std::sort( IP2Ds_2.begin(),IP2Ds_2.end(),std::greater<float>() );
			   std::sort( etaRels.begin(),etaRels.end() ); //std::sort sorts in ascending order by default
			   std::sort( etaRels_1.begin(),etaRels_1.end() ); 		
			   std::sort( etaRels_2.begin(),etaRels_2.end() );
			   int num_1 = IP2Ds_1.size();
			   int num_2 = IP2Ds_2.size();
			   int numTracks = FatJetInfo_TagVarCSV_jetNTracks[i];
			   int numEtaRelTracks = FatJetInfo_TagVarCSV_jetNTracksEtaRel[i];
			   float dummyTrack = -50.;
			   float dummyPtRel = 0.;
			   float dummyEtaRel = -50.;
			   switch(numTracks){
				   case 0:

					   trackSip2dSig_0 = dummyTrack;
					   trackSip2dSig_1 = dummyTrack;
					   trackSip2dSig_2 = dummyTrack;
					   trackSip2dSig_3 = dummyTrack;

					   trackSip3dSig_0 = dummyTrack;
					   trackSip3dSig_1 = dummyTrack;
					   trackSip3dSig_2 = dummyTrack;
					   trackSip3dSig_3 = dummyTrack;

					   trackSipdSig_0 = dummyTrack;
					   trackSipdSig_1 = dummyTrack;
					   trackSipdSig_2 = dummyTrack;
					   trackSipdSig_3 = dummyTrack;

					   trackPtRel_0 = dummyPtRel;
					   trackPtRel_1 = dummyPtRel;
					   trackPtRel_2 = dummyPtRel;
					   trackPtRel_3 = dummyPtRel;
					   break;

				   case 1:


					   trackSip3dSig_0 = IP3Ds.at(0);
					   trackSip3dSig_1 = dummyTrack;
					   trackSip3dSig_2 = dummyTrack;
					   trackSip3dSig_3 = dummyTrack;

					   trackPtRel_0 = PtRel.at(0);
					   trackPtRel_1 = dummyPtRel;
					   trackPtRel_2 = dummyPtRel;
					   trackPtRel_3 = dummyPtRel;
					   break;

				   case 2:


					   trackSip3dSig_0 = IP3Ds.at(0);
					   trackSip3dSig_1 = IP3Ds.at(1);
					   trackSip3dSig_2 = dummyTrack;
					   trackSip3dSig_3 = dummyTrack;


					   trackPtRel_0 = PtRel.at(0);
					   trackPtRel_1 = PtRel.at(1);
					   trackPtRel_2 = dummyPtRel;
					   trackPtRel_3 = dummyPtRel;
					   break;

				   case 3:



					   trackSip3dSig_0 = IP3Ds.at(0);
					   trackSip3dSig_1 = IP3Ds.at(1);
					   trackSip3dSig_2 = IP3Ds.at(2);
					   trackSip3dSig_3 = dummyTrack;

					   trackPtRel_0 = PtRel.at(0);
					   trackPtRel_1 = PtRel.at(1);
					   trackPtRel_2 = PtRel.at(2);
					   trackPtRel_3 = dummyPtRel;
					   break;

				   default:


					   trackSip3dSig_0 = IP3Ds.at(0);
					   trackSip3dSig_1 = IP3Ds.at(1);
					   trackSip3dSig_2 = IP3Ds.at(2);
					   trackSip3dSig_3 = IP3Ds.at(3);

					   trackPtRel_0 = PtRel.at(0);
					   trackPtRel_1 = PtRel.at(1);
					   trackPtRel_2 = PtRel.at(2);
					   trackPtRel_3 = PtRel.at(3);

			   } // end switch on number of tracks for IP

			   switch(contTrk){
				   case 0:


					   trackSipdSig_0 = dummyTrack;
					   trackSipdSig_1 = dummyTrack;
					   trackSipdSig_2 = dummyTrack;
					   trackSipdSig_3 = dummyTrack;	

					   break;

				   case 1:


					   trackSipdSig_0 = IP2Ds.at(0);
					   trackSipdSig_1 = dummyTrack;
					   trackSipdSig_2 = dummyTrack;
					   trackSipdSig_3 = dummyTrack;	

					   break;

				   case 2:

					   trackSipdSig_0 = IP2Ds.at(0);
					   trackSipdSig_1 = IP2Ds.at(1);
					   trackSipdSig_2 = dummyTrack;
					   trackSipdSig_3 = dummyTrack;

					   break;

				   case 3:

					   trackSipdSig_0 = IP2Ds.at(0);
					   trackSipdSig_1 = IP2Ds.at(1);
					   trackSipdSig_2 = IP2Ds.at(2);
					   trackSipdSig_3 = dummyTrack;	

					   break;

				   default:


					   trackSipdSig_0 = IP2Ds.at(0);
					   trackSipdSig_1 = IP2Ds.at(1);
					   trackSipdSig_2 = IP2Ds.at(2);
					   trackSipdSig_3 = IP2Ds.at(3);	


			   } // end switch on number of tracks for IP
			   switch(num_1){
				   case 0:


					   trackSipdSig_0_0 = dummyTrack;
					   trackSipdSig_0_1 = dummyTrack;
					   trackSipdSig_0_2 = dummyTrack;
					   trackSipdSig_0_3 = dummyTrack;

					   trackEtaRel_0_0 = dummyEtaRel;
					   trackEtaRel_0_1 = dummyEtaRel;
					   trackEtaRel_0_2 = dummyEtaRel;	 

					   break;

				   case 1:


					   trackSipdSig_0_0 = IP2Ds_1.at(0);
					   trackSipdSig_0_1 = dummyTrack;
					   trackSipdSig_0_2 = dummyTrack;
					   trackSipdSig_0_3 = dummyTrack;

					   trackEtaRel_0_0 = etaRels_1.at(0);
					   trackEtaRel_0_1 = dummyEtaRel;
					   trackEtaRel_0_2 = dummyEtaRel;

					   break;

				   case 2:

					   trackSipdSig_0_0 = IP2Ds_1.at(0);
					   trackSipdSig_0_1 = IP2Ds_1.at(1);
					   trackSipdSig_0_2 = dummyTrack;
					   trackSipdSig_0_3 = dummyTrack;

					   trackEtaRel_0_0 = etaRels_1.at(0);
					   trackEtaRel_0_1 = etaRels_1.at(1);
					   trackEtaRel_0_2 = dummyEtaRel;	

					   break;

				   case 3:

					   trackSipdSig_0_0 = IP2Ds_1.at(0);
					   trackSipdSig_0_1 = IP2Ds_1.at(1);
					   trackSipdSig_0_2 = IP2Ds_1.at(2);
					   trackSipdSig_0_3 = dummyTrack;



					   trackEtaRel_0_0 = etaRels_1.at(0);
					   trackEtaRel_0_1 = etaRels_1.at(1);
					   trackEtaRel_0_2 = etaRels_1.at(2);

					   break;

				   default:


					   trackSipdSig_0_0 = IP2Ds_1.at(0);
					   trackSipdSig_0_1 = IP2Ds_1.at(1);
					   trackSipdSig_0_2 = IP2Ds_1.at(2);
					   trackSipdSig_0_3 = IP2Ds_1.at(3);


					   trackEtaRel_0_0 = etaRels_1.at(0);
					   trackEtaRel_0_1 = etaRels_1.at(1);
					   trackEtaRel_0_2 = etaRels_1.at(2);	


			   }	

			   switch(num_2){
				   case 0:


					   trackSipdSig_1_0 = dummyTrack;
					   trackSipdSig_1_1 = dummyTrack;
					   trackSipdSig_1_2 = dummyTrack;
					   trackSipdSig_1_3 = dummyTrack;

					   trackEtaRel_1_0 = dummyEtaRel;
					   trackEtaRel_1_1 = dummyEtaRel;
					   trackEtaRel_1_2 = dummyEtaRel;

					   break;

				   case 1:


					   trackSipdSig_1_0 = IP2Ds_2.at(0);
					   trackSipdSig_1_1 = dummyTrack;
					   trackSipdSig_1_2 = dummyTrack;
					   trackSipdSig_1_3 = dummyTrack;


					   trackEtaRel_1_0 = etaRels_2.at(0);
					   trackEtaRel_1_1 = dummyEtaRel;
					   trackEtaRel_1_2 = dummyEtaRel;


					   break;

				   case 2:

					   trackSipdSig_1_0 = IP2Ds_2.at(0);
					   trackSipdSig_1_1 = IP2Ds_2.at(1);
					   trackSipdSig_1_2 = dummyTrack;
					   trackSipdSig_1_3 = dummyTrack;

					   trackEtaRel_1_0 = etaRels_2.at(0);
					   trackEtaRel_1_1 = etaRels_2.at(1);
					   trackEtaRel_1_2 = dummyEtaRel;	



					   break;

				   case 3:

					   trackSipdSig_1_0 = IP2Ds_2.at(0);
					   trackSipdSig_1_1 = IP2Ds_2.at(1);
					   trackSipdSig_1_2 = IP2Ds_2.at(2);
					   trackSipdSig_1_3 = dummyTrack;

					   trackEtaRel_1_0 = etaRels_2.at(0);
					   trackEtaRel_1_1 = etaRels_2.at(1);
					   trackEtaRel_1_2 = etaRels_2.at(2);


					   break;

				   default:


					   trackSipdSig_1_0 = IP2Ds_2.at(0);
					   trackSipdSig_1_1 = IP2Ds_2.at(1);
					   trackSipdSig_1_2 = IP2Ds_2.at(2);
					   trackSipdSig_1_3 = IP2Ds_2.at(3);

					   trackEtaRel_1_0 = etaRels_2.at(0);
					   trackEtaRel_1_1 = etaRels_2.at(1);
					   trackEtaRel_1_2 = etaRels_2.at(2);	


			   }


			   switch(contTrk_tight){
				   case 0:


					   trackSip2dSig_0 = dummyTrack;
					   trackSip2dSig_1 = dummyTrack;
					   trackSip2dSig_2 = dummyTrack;
					   trackSip2dSig_3 = dummyTrack;

					   break;

				   case 1:


					   trackSip2dSig_0 = IPD.at(0);
					   trackSip2dSig_1 = dummyTrack;
					   trackSip2dSig_2 = dummyTrack;
					   trackSip2dSig_3 = dummyTrack;

					   break;

				   case 2:

					   trackSip2dSig_0 = IPD.at(0);
					   trackSip2dSig_1 = IPD.at(1);
					   trackSip2dSig_2 = dummyTrack;
					   trackSip2dSig_3 = dummyTrack;

					   break;

				   case 3:

					   trackSip2dSig_0 = IPD.at(0);
					   trackSip2dSig_1 = IPD.at(1);
					   trackSip2dSig_2 = IPD.at(2);
					   trackSip2dSig_3 = dummyTrack;

					   break;

				   default:


					   trackSip2dSig_0 = IPD.at(0);
					   trackSip2dSig_1 = IPD.at(1);
					   trackSip2dSig_2 = IPD.at(2);
					   trackSip2dSig_3 = IPD.at(3);


			   }	
			   switch(numEtaRelTracks){
				   case 0:

					   trackEtaRel_0 = dummyEtaRel;
					   trackEtaRel_1 = dummyEtaRel;
					   trackEtaRel_2 = dummyEtaRel;
					   break;

				   case 1:

					   trackEtaRel_0 = etaRels.at(0);
					   trackEtaRel_1 = dummyEtaRel;
					   trackEtaRel_2 = dummyEtaRel;
					   break;

				   case 2:

					   trackEtaRel_0 = etaRels.at(0);
					   trackEtaRel_1 = etaRels.at(1);
					   trackEtaRel_2 = dummyEtaRel;
					   break;

				   default:

					   trackEtaRel_0 = etaRels.at(0);
					   trackEtaRel_1 = etaRels.at(1);
					   trackEtaRel_2 = etaRels.at(2);

			   } //end switch on number of etarel tracks



			   if( (FatJetInfo_Jet_tau1_vertexMass[i]>0 && FatJetInfo_Jet_tau2_vertexMass[i]>0)){
				   tau00_trackEtaRel_0 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_0[i],FatJetInfo_Jet_tau2_trackEtaRel_0[i]);
				   tau00_trackEtaRel_1 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_1[i],FatJetInfo_Jet_tau2_trackEtaRel_1[i]);
				   tau00_trackEtaRel_2 = TMath::Max(FatJetInfo_Jet_tau1_trackEtaRel_2[i],FatJetInfo_Jet_tau2_trackEtaRel_2[i]);


				   tau01_trackEtaRel_0 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_0[i],FatJetInfo_Jet_tau2_trackEtaRel_0[i]);
				   tau01_trackEtaRel_1 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_1[i],FatJetInfo_Jet_tau2_trackEtaRel_1[i]);
				   tau01_trackEtaRel_2 = TMath::Min(FatJetInfo_Jet_tau1_trackEtaRel_2[i],FatJetInfo_Jet_tau2_trackEtaRel_2[i]);


			   }  else if(FatJetInfo_Jet_tau2_vertexMass[i]<0 && FatJetInfo_Jet_tau1_vertexMass[i]>0){

				   tau00_trackEtaRel_0 = FatJetInfo_Jet_tau1_trackEtaRel_0[i];
				   tau00_trackEtaRel_1 = FatJetInfo_Jet_tau1_trackEtaRel_1[i];
				   tau00_trackEtaRel_2 = FatJetInfo_Jet_tau1_trackEtaRel_2[i];
				   tau01_trackEtaRel_0 = trackEtaRel_1_0;
				   tau01_trackEtaRel_1 = trackEtaRel_1_1;
				   tau01_trackEtaRel_2= trackEtaRel_1_2;



			   }else if(FatJetInfo_Jet_tau1_vertexMass[i]<0 && FatJetInfo_Jet_tau2_vertexMass[i]>0){
				   tau00_trackEtaRel_0 = FatJetInfo_Jet_tau2_trackEtaRel_0[i];
				   tau00_trackEtaRel_1 = FatJetInfo_Jet_tau2_trackEtaRel_1[i];
				   tau00_trackEtaRel_2 = FatJetInfo_Jet_tau2_trackEtaRel_2[i];
				   tau01_trackEtaRel_0 = trackEtaRel_0_0;
				   tau01_trackEtaRel_1 = trackEtaRel_0_1;
				   tau01_trackEtaRel_2 =trackEtaRel_0_2;


			   }
			   else{
				   tau01_trackEtaRel_0 = trackEtaRel_1_0;
				   tau01_trackEtaRel_1 = trackEtaRel_1_1;
				   tau01_trackEtaRel_2= trackEtaRel_1_2;
				   tau00_trackEtaRel_0 = trackEtaRel_0_0;
				   tau00_trackEtaRel_1 = trackEtaRel_0_1;
				   tau00_trackEtaRel_2 =trackEtaRel_0_2;



			   }


			   int iSubJet1 = FatJetInfo_SubJetIdx_Pruned[FatJetInfo_Jet_nFirstSJ_Pruned[i]];
			   int iSubJet2 = FatJetInfo_SubJetIdx_Pruned[FatJetInfo_Jet_nFirstSJ_Pruned[i]+1];
			   if( PrunedSubJetInfo_Jet_pt[iSubJet1] < PrunedSubJetInfo_Jet_pt[iSubJet2] )
			   {
				   iSubJet1 = iSubJet1 + iSubJet2;
				   iSubJet2 = iSubJet1 - iSubJet2; // iSubJet2 = iSubJet1
				   iSubJet1 = iSubJet1 - iSubJet2;
			   }

			   SubJet1_pt        = PrunedSubJetInfo_Jet_pt[iSubJet1];
			   SubJet1_eta       = PrunedSubJetInfo_Jet_eta[iSubJet1];
			   SubJet1_phi       = PrunedSubJetInfo_Jet_phi[iSubJet1];
			   SubJet1_mass      = PrunedSubJetInfo_Jet_mass[iSubJet1];
			   SubJet1_flavour   = PrunedSubJetInfo_Jet_flavour[iSubJet1];
			   SubJet1_nbHadrons = PrunedSubJetInfo_Jet_nbHadrons[iSubJet1];
			   SubJet1_JP        = PrunedSubJetInfo_Jet_Proba[iSubJet1];
			   SubJet1_JBP       = PrunedSubJetInfo_Jet_Bprob[iSubJet1];
			   SubJet1_CSV       = PrunedSubJetInfo_Jet_CombSvx[iSubJet1];
			   SubJet1_CSVIVF    = PrunedSubJetInfo_Jet_CombIVF[iSubJet1];
			   std::vector<float> etaRels1; // stores |trackEtaRel|!
			   /*for (int iTrk = PrunedSubJetInfo_Jet_nFirstTrkEtaRelTagVarCSV[iSubJet1]; iTrk < PrunedSubJetInfo_Jet_nLastTrkEtaRelTagVarCSV[iSubJet1]; ++iTrk)
				   etaRels1.push_back( fabs(PrunedSubJetInfo_TagVarCSV_trackEtaRel[iTrk]) );

			   std::sort( etaRels1.begin(),etaRels1.end() ); //std::sort sorts in ascending order by default

			  int numEtaRelTracks1 = PrunedSubJetInfo_TagVarCSV_jetNTracksEtaRel[iSubJet1];

			     switch(numEtaRelTracks1){
			     case 0:

			     TagVarCSV1_trackEtaRel_0 = dummyEtaRel;
			     TagVarCSV1_trackEtaRel_1 = dummyEtaRel;
			     TagVarCSV1_trackEtaRel_2 = dummyEtaRel;
			     break;

			     case 1:

			     TagVarCSV1_trackEtaRel_0 = etaRels1.at(0);
			     TagVarCSV1_trackEtaRel_1 = dummyEtaRel;
			     TagVarCSV1_trackEtaRel_2 = dummyEtaRel;
			     break;

			     case 2:

			     TagVarCSV1_trackEtaRel_0 = etaRels1.at(0);
			     TagVarCSV1_trackEtaRel_1 = etaRels1.at(1);
			     TagVarCSV1_trackEtaRel_2 = dummyEtaRel;
			     break;

			     default:

			     TagVarCSV1_trackEtaRel_0 = etaRels1.at(0);
			     TagVarCSV1_trackEtaRel_1 = etaRels1.at(1);
			     TagVarCSV1_trackEtaRel_2 = etaRels1.at(2);

			     } 
			    */
			   SubJet2_pt        = PrunedSubJetInfo_Jet_pt[iSubJet2];
			   SubJet2_eta       = PrunedSubJetInfo_Jet_eta[iSubJet2];
			   SubJet2_phi       = PrunedSubJetInfo_Jet_phi[iSubJet2];
			   SubJet2_mass      = PrunedSubJetInfo_Jet_mass[iSubJet2];
			   SubJet2_flavour   = PrunedSubJetInfo_Jet_flavour[iSubJet2];
			   SubJet2_nbHadrons = PrunedSubJetInfo_Jet_nbHadrons[iSubJet2];
			   SubJet2_JP        = PrunedSubJetInfo_Jet_Proba[iSubJet2];
			   SubJet2_JBP       = PrunedSubJetInfo_Jet_Bprob[iSubJet2];
			   SubJet2_CSV       = PrunedSubJetInfo_Jet_CombSvx[iSubJet2];
			   SubJet2_CSVIVF    = PrunedSubJetInfo_Jet_CombIVF[iSubJet2];
			   SubJet_csv=-1;

			   if(SubJet1_pt>0) SubJet_csv =PrunedSubJetInfo_Jet_CombIVF[iSubJet1];		
			   if(SubJet2_pt>0 && SubJet1_pt>0) SubJet_csv = TMath::Min(PrunedSubJetInfo_Jet_CombIVF[iSubJet2],PrunedSubJetInfo_Jet_CombIVF[iSubJet1]);	
			   if(SubJet_csv <-1 || SubJet_csv> 1.) SubJet_csv =-1;
			   FatJet_csv  =  FatJetInfo_Jet_CombIVF[i]; 
			   if(FatJet_csv <-1 || FatJet_csv> 1.) FatJet_csv =-1;

/*			   if(FatJetInfo_Jet_z_ratio[i]!=z_ratio) std::cout<<FatJetInfo_Jet_z_ratio[i]<< " z_Ratio "<<z_ratio<<std::endl; 
			   if(FatJetInfo_Jet_tau1_trackEtaRel_0[i]!=tau1_trackEtaRel_0) std::cout<<FatJetInfo_Jet_tau1_trackEtaRel_0[i]<<" trackEtaRel_0 "<<tau0_trackEtaRel_0<<std::endl; 
			   if(FatJetInfo_Jet_tau1_vertexMass[i]!=tau_vertexMass_0) std::cout<<FatJetInfo_Jet_tau1_vertexMass[i]<<" vertexMass_0 "<<tau_vertexMass_0<<std::endl;
			   if(FatJetInfo_Jet_trackSip2dSigAboveBottom_0[i]!=trackSip2dSigAboveBottom_0) std::cout<<FatJetInfo_Jet_trackSip2dSigAboveBottom_0[i]<<" trackSip2dSigAboveBottom_0 "<<trackSip2dSigAboveBottom_0<<std::endl;
			   if(FatJetInfo_Jet_tau1_vertexEnergyRatio[i]!=tau_vertexEnergyRatio_0) std::cout<<FatJetInfo_Jet_tau1_vertexEnergyRatio[i]<<" WARNING "<<tau_vertexEnergyRatio_0<<std::endl;

*/
			   trOUT->Fill();
			   }
		   }// jet loop  
	   }
	   std::cout<<" sV =0 " <<(double) cc0/nsel <<std::endl;
	   std::cout<<" sV =1 " <<(double) cc1/nsel <<std::endl;
	   std::cout<<" sV =2 " <<(double) cc2/nsel <<std::endl;
	   std::cout<<" sV =3 " <<(double) cc3/nsel <<std::endl;
	   std::cout<<" sV =4 " <<(double) cc4/nsel <<std::endl;




	   outf->cd();     
	   trOUT->Write(); 
	   h_multiplicity->Write();	
	   outf->Close();  
   }     

