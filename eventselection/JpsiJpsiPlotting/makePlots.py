from CMSPLOTS.myFunction import DrawHistos
import ROOT

ROOT.gROOT.SetBatch(True)

# change this input address
finput = "../myhbk.root"

f = ROOT.TFile(finput)

tolPalette = [
	"#332288", 		#persian indigo
	"#117733", 		#deep green/camarone
	"#44AA99", 		#sea green/gossamer
	"#88CCEE", 		#sky blue
	"#DDCC77", 		#soft yellow
	"#CC6677", 		#contessa red
	"#AA4499", 		#dark magenta
	"#882255"		#rose bud cherry
	]

hname = "myFourMuonmass_mass"
h = f.Get(hname)
h.SetLineColor(1)
h.SetMarkerStyle(20)

histList = [h]
histLabels = ["Data"]
binWidth = histList[0].GetBinWidth()
binUnit = "MeV"
# binUnit = "GeV"
DrawHistos(histList, histLabels, 6, 10, "m(#mu^{+}#mu^{-}#mu^{+}#mu^{-}) [%s/c^{2}]"%binUnit, 0, 50, "Events / (%s %s/c^{2})"%(binWidth, binUnit), outputname = "fourmuon_mass", dology=False, legendoptions=["ep"], nMaxDigits=3)

hname1 = "myDiMuon1mass_mass"
hname2 = "myDiMuon2mass_mass"

h1 = f.Get(hname1)
h2 = f.Get(hname2)
h1.SetLineColor(1)
h1.SetMarkerColor(1)
h1.SetMarkerStyle(20)
h2.SetLineColor(2)
h2.SetMarkerColor(2)
h2.SetMarkerStyle(20)
DrawHistos([h1, h2], ["1st J/psi", "2nd J/psi"], 2.5, 3.5, "m(#mu^{+}#mu^{-}) [GeV/c^{2}]", 0, 2000, "Events / (0.1 GeV/c^{2})", outputname = "dimuon_mass", dology=False, legendoptions=["ep", "ep"], addOverflow=True, addUnderflow=True, nMaxDigits=3)

histList_diMuon = [h1, h2]

for hist, color in 

h1.SetLineColor(TColor.GetColor())