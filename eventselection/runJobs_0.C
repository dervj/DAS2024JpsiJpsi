void runJobs_0()
{
	gSystem->Load("myntuple_C.so");
	TChain * chain = new TChain("/mkcands/X_data","");

	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016B/Charmonium/crab_2016BReRecoskim/210122_001306/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016C/Charmonium/crab_2016CReRecoskim/210122_001412/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016D/Charmonium/crab_2016DReRecoskim/210122_001527/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016E/Charmonium/crab_2016EReRecoskim/210122_001648/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016F/Charmonium/crab_2016FReRecoskim/210122_001728/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016G/Charmonium/crab_2016GReRecoskim/210122_001820/0000/*root");
	chain->Add("/eos/uscms/store/user/lpcbphy/yik/Charmoniumskimv6/2016/2016H/Charmonium/crab_2016HReRecoskim/210122_001848/0000/*root");

	myntuple a(chain);
	a.Loop();

}
