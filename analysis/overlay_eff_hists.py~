import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFileDefault', type=str, default='./tau_gun/MAIA/TauFinderOutputs/tau_ana_default.root')
parser.add_argument('--inputFileLoose', type=str, default='./tau_gun/MAIA/TauFinderOutputs/tau_ana_loose.root')

args = parser.parse_args()

root_file_default = ROOT.TFile(args.inputFileDefault, 'READ')
root_file_loose = ROOT.TFile(args.inputFileLoose, 'READ')



nums = ['1', '3']
names = ['pt', 'phi', 'theta']
variables = ['p_{T}', '#phi', '#theta']
units = ['GeV/c', 'rad', 'rad']

for num in nums:
    for i in range(3):

        hDefault = root_file_default.Get('reco_' + num + 'p_' + names[i] + '_eff')
        hLoose = root_file_loose.Get('reco_' + num + 'p_' + names[i] + '_eff')
        hPion = root_file_loose.Get(num + 'p_pi_' + names[i] + '_eff')

        hDefault.SetTitle(num + '-Prong Reconstruction Efficiencies vs ' + variables[i])
        hDefault.GetXaxis().SetTitle('True Visible #tau^{-} ' + variables[i] + ' [' + units[i] + ']')
        hDefault.GetYaxis().SetRangeUser(0.0, 1.2)

        hDefault.SetLineColor(ROOT.kRed)
        hLoose.SetLineColor(ROOT.kBlue)
        hPion.SetLineColor(ROOT.kGreen)

        hDefault.SetLineWidth(2)
        hLoose.SetLineWidth(2)
        hPion.SetLineWidth(2)

        hDefault.SetMarkerStyle(ROOT.kFullDotLarge)
        hLoose.SetMarkerStyle(ROOT.kMultiply)

        hDefault.SetMarkerColor(ROOT.kRed)
        hLoose.SetMarkerColor(ROOT.kBlue)

        hDefault.SetMarkerSize(1.5)
        hLoose.SetMarkerSize(1.5)

        hDefault.SetStats(0)
        hLoose.SetStats(0)
        hPion.SetStats(0)

        c = ROOT.TCanvas('c', 'overlay', 800, 600)

        hDefault.Draw()
        hLoose.Draw('SAME')
        hPion.Draw('SAME')

        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.04)
        text.SetTextAlign(13)
        text.DrawLatex(0.12, 0.85, "#bf{#it{MAIA Detector Concept}}")
        text.DrawLatex(0.12, 0.80, "Simulated #tau^{-} Events")

        c.Update()

        legend = ROOT.TLegend(0.75, 0.75, 0.95, 0.90)
        legend.AddEntry(hPion, 'Pion Reco', 'lp')
        legend.AddEntry(hLoose, 'Loose Tau Reco', 'lp')
        legend.AddEntry(hDefault, 'Default Tau Reco', 'lp')
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        c.Update()

        filename = hDefault.GetTitle() + '.png'
        c.SaveAs(filename)
