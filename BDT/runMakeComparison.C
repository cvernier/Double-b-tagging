{


	gROOT->ProcessLine(".L makeComparisonBackground.C");

	//makeComparisonBackground("z_ratio", 0., 40., "z variable");
/*	makeComparisonBackground("PFLepton_ptrel"); 

	makeComparisonBackground("PFLepton_ratio"); 
*/

	
//        makeComparisonBackground("tau_vertexMass_corrected_0", 0.,6., "SV_{0} mass (GeV)"); 
//        makeComparisonBackground("tau_vertexMass_corrected_1", 0., 6., "SV_{1} mass (GeV)");
        //makeComparisonBackground("tau_vertexEnergyRatio_0", 0., 4., "SV_{0} p_{T} ratio"); 
/*	makeComparisonBackground("tau_vertexEnergyRatio_1", 0., 4., "SV_{1} Energy ratio"); 
	makeComparisonBackground("tau_vertexDeltaR_0", 0., 1., "dR(SV_{0}, #tau)"); 
	//makeComparisonBackground("tau_vertexDeltaR_1", 0., 1.);  
     

	makeComparisonBackground("tau_flightDistance2dSig_1", 0., 60, "SV_{1} flight distance sig"); 
	makeComparisonBackground("tau_flightDistance2dSig_0", 0., 120, "SV_{0} flight distance sig");


	makeComparisonBackground("trackSip2dSig_3", 0., 20, "IP Sig 4th track"); 
	makeComparisonBackground("trackSip2dSig_1", 0., 50, "IP Sig 2nd track"); 
	makeComparisonBackground("trackSip2dSig_2", 0., 50, "IP Sig 3rd track"); 
	makeComparisonBackground("trackSip2dSig_0", 0. ,50., "IP Sig 1st track"); 
      */  
       // makeComparisonBackground("trackSip2dSigAboveBottom_0", 0., 20.,"2D SIP for first track above bottom");  
      /*  makeComparisonBackground("trackSip2dSigAboveCharm_0", 0., 50, "IP Sig First Track Above Charm"); 
	makeComparisonBackground("trackSip2dSigAboveBottom_1", 0., 30, "IP Sig Second Track Above Bottom");
	makeComparisonBackground("tau1_trackEtaRel_0", 0., 8.,  "EtaRel 1st track for SV_{1}");
	makeComparisonBackground("tau1_trackEtaRel_1", 0., 12.,  "EtaRel 2nd track for SV_{1}");
	makeComparisonBackground("tau1_trackEtaRel_2", 0.,12.,  "EtaRel 3rd track for SV_{1}");
	makeComparisonBackground("tau0_trackEtaRel_0", 0., 8., "EtaRel 1st track for SV_{0}");
        makeComparisonBackground("tau0_trackEtaRel_1", 0., 12.,"EtaRel 2nd track for SV_{0}");
	makeComparisonBackground("tau0_trackEtaRel_2", 0., 12.,"EtaRel 3rd track for SV_{0}");
	makeComparisonBackground("massPruned", 50., 200.,"Pruned mass (GeV)");
	makeComparisonBackground("tau2/tau1", 0., 1.,"#tau_2/#tau_1");
	


*/
	//makeComparisonBackground("ptPruned", 0., 2500.,"p$_T$ (GeV)");
	makeComparisonBackground("nSV", 0., 8., "number of SV"); 
//	makeComparisonBackground("jetNTracks", 0., 50., "Number of tracks");
//i
//
	/*makeComparisonBackground("trackSipdSig_0_1", 0., 20, "IP Sig 1st track associated to #tau_{2}"); 
        makeComparisonBackground("trackSipdSig_1_1", 0., 50, "IP Sig 2nd track associated to #tau_{2}"); 
        makeComparisonBackground("trackSipdSig_0_0", 0., 50, "IP Sig 1st track associated to #tau_{1}"); 
        makeComparisonBackground("trackSipdSig_1_0", 0. ,50., "IP Sig 2nd track associated to #tau_{1}"); 
*/

	gROOT->ProcessLine(".q");

} 
