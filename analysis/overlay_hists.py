import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFileDefault', type=str, default='./tau_gun/MAIA/TauFinderOutputs/tau_ana_default.root')
parser.add_argument('--inputFileLoose', type=str, default='./tau_gun/MAIA/TauFinderOutputs/tau_ana_loose.root')

args = parser.parse_args()

root_file_default = ROOT.TFile(args.inputFileDefault, 'READ')
root_file_loose = ROOT.TFile(args.inputFileLoose, 'READ')
 
hDefault = root_file_default.Get('reco_1p_pt_eff')
#hDefault3P = root_file_default.Get('reco_1p_phi_eff')
hLoose = root_file_loose.Get('reco_1p_pt_eff')
hPion = root_file_loose.Get('1p_pi_pt_eff')

hDefault.SetTitle('1-Prong Reconstruction Efficiencies vs p_{T}')
hDefault.GetXaxis().SetTitle('True Visible #tau^{-} p_{T} [GeV/c]')
hDefault.GetYaxis().SetRangeUser(0.15, 1.2)

hDefault.SetLineColor(ROOT.kRed)
#hDefault3P.SetLineColor(ROOT.kOrange)
hLoose.SetLineColor(ROOT.kBlue)
hPion.SetLineColor(ROOT.kGreen)

hDefault.SetLineWidth(2)
#hDefault3P.SetLineWidth(2)
hLoose.SetLineWidth(2)
hPion.SetLineWidth(2)

hDefault.SetMarkerStyle(ROOT.kFullDotLarge)
hLoose.SetMarkerStyle(ROOT.kMultiply)

hDefault.SetMarkerColor(ROOT.kRed)
hLoose.SetMarkerColor(ROOT.kBlue)

hDefault.SetMarkerSize(1.5)
hLoose.SetMarkerSize(1.5)

hDefault.SetStats(0)
#hDefault3P.SetStats(0)
hLoose.SetStats(0)
hPion.SetStats(0)

c = ROOT.TCanvas('c', 'overlay', 800, 600)

hDefault.Draw()
#hDefault3P.Draw('SAME')
hLoose.Draw('SAME')
hPion.Draw('SAME')

text = ROOT.TLatex()
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.04)
text.SetTextAlign(13)
text.DrawLatex(0.15, 0.89, "MAIA Detector Concept")
text.DrawLatex(0.15, 0.84, "Simulated #tau^{-} Events")

c.Update()

legend = ROOT.TLegend(0.70, 0.75, 0.90, 0.90)
#legend = ROOT.TLegend(0.75, 0.85, 1.05, 1.0)
legend.AddEntry(hPion, 'Pion Reco', 'lp')
legend.AddEntry(hLoose, 'Loose Tau Reco', 'lp')
legend.AddEntry(hDefault, 'Default Tau Reco', 'lp')
#legend.AddEntry(hDefault3P, '3-Prong', 'l')
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.Draw()

c.Update()

filename = hDefault.GetTitle() + '.png'
c.SaveAs(filename)
