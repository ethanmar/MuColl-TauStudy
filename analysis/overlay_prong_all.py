import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFileTau', type=str, default='tau_eff_custom.root')
parser.add_argument('--inputFilePi', type=str, default='pi_eff_custom.root')
args = parser.parse_args()

root_file_tau = ROOT.TFile(args.inputFileTau, 'READ')
root_file_pi = ROOT.TFile(args.inputFilePi, 'READ')

names = ['pt', 'phi', 'theta']
variables = ['#it{p}_{T}^{vis}', '#phi^{vis}', '#theta^{vis}']
units = ['GeV/c', 'rad', 'rad']

for i in range(3):
    hTau1P0N = root_file_tau.Get('1p0n_tau_' + names[i] + '_eff')
    hPion1P0N = root_file_pi.Get('1p0n_pi_' + names[i] + '_eff')
    hTau3P0N = root_file_tau.Get('3p0n_tau_' + names[i] + '_eff')
    hPion3P0N = root_file_pi.Get('3p0n_pi_' + names[i] + '_eff')

    hPion1P0N.SetTitle('')
    hPion1P0N.SetName('1p0n_3p0n_' + names[i] + '_eff')
    hPion1P0N.GetXaxis().SetTitle('Tau Truth ' + variables[i] + ' [' + units[i] + ']')
    hPion1P0N.GetYaxis().SetRangeUser(0.0, 1.3)
    
    hPion1P0N.SetLineColor(ROOT.kGreen)
    hTau1P0N.SetLineColor(ROOT.kBlue)
    hPion3P0N.SetLineColor(ROOT.kCyan)
    hTau3P0N.SetLineColor(ROOT.kMagenta)
    
    hPion1P0N.SetLineWidth(2)
    hTau1P0N.SetLineWidth(2)
    hPion3P0N.SetLineWidth(2)
    hTau3P0N.SetLineWidth(2)
    
    hPion1P0N.SetMarkerStyle(47)
    hTau1P0N.SetMarkerStyle(43)
    hPion3P0N.SetMarkerStyle(34)
    hTau3P0N.SetMarkerStyle(29)
    
    hPion1P0N.SetMarkerColor(ROOT.kGreen)
    hTau1P0N.SetMarkerColor(ROOT.kBlue)
    hPion3P0N.SetMarkerColor(ROOT.kCyan)
    hTau3P0N.SetMarkerColor(ROOT.kMagenta)
    
    hPion1P0N.SetMarkerSize(1.5)
    hTau1P0N.SetMarkerSize(1.5)
    hPion3P0N.SetMarkerSize(1.5)
    hTau3P0N.SetMarkerSize(1.5)
    
    hPion1P0N.SetStats(0)
    hTau1P0N.SetStats(0)
    hPion3P0N.SetStats(0)
    hTau3P0N.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)
    c.SetGridy()
    c.SetTicky()
        
    hPion1P0N.Draw()
    hTau1P0N.Draw('SAME')
    hPion3P0N.Draw('SAME')
    hTau3P0N.Draw('SAME')
    
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.SetTextAlign(13)
    text.DrawLatex(0.13, 0.87, "#bf{#it{MAIA}} Work in Progress")
    text.DrawLatex(0.13, 0.83, "Charged Tau Gun (No BIB)")
    c.Update()
        
    legend = ROOT.TLegend(0.70, 0.75, 0.92, 0.90)
    legend.AddEntry(hPion1P0N, '1P0N Pion Reco', 'lp')
    legend.AddEntry(hTau1P0N, '1P0N Tau Reco', 'lp')
    legend.AddEntry(hPion3P0N, '3P0N Pion Reco', 'lp')
    legend.AddEntry(hTau3P0N, '3P0N Tau Reco', 'lp')        
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()
    
    c.Update()
    
    filename = hPion1P0N.GetName() + '.png'
    c.SaveAs(filename)
