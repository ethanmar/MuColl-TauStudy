import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFile', type=str, default='tau_eff.root')
args = parser.parse_args()

root_file = ROOT.TFile(args.inputFile, 'READ')

names = ['pt', 'phi']
variables = ['#it{p_{T}^{vis}}', '#phi']
units = ['GeV/c', 'rad']
regions = ['centbarrel', 'transition', 'endcap']

for i in range(3):
    for j in range(2):
        hTau1P0N = root_file.Get('1p0n_tau_' + names[j] + '_eff_' + regions[i])
        hTau1P1N = root_file.Get('1p1n_tau_' + names[j] + '_eff_' + regions[i])
        hTau1P2N = root_file.Get('1p2n_tau_' + names[j] + '_eff_' + regions[i])
        hTau1PXN = root_file.Get('1pXn_tau_' + names[j] + '_eff_' + regions[i])
        hTau3P0N = root_file.Get('3p0n_tau_' + names[j] + '_eff_' + regions[i])

        hTau1P0N.SetTitle('')
        hTau1P0N.SetName('tau_' + names[j] + '_eff_' + regions[i])
        if j == 0:
            hTau1P0N.GetXaxis().SetTitle("#font[42]{Visible True Tau}#kern[-0.5]{ }#it{p_{T}}#kern[-0.5]{ }#font[42]{[GeV/c]}")
        else:
            hTau1P0N.GetXaxis().SetTitle("#font[42]{Visible True Tau}#kern[-0.5]{ }#phi#kern[-0.5]{ }#font[42]{[rad]}")
        hTau1P0N.GetYaxis().SetRangeUser(0.0, 1.3)
    
        hTau1P0N.SetLineColor(ROOT.kGreen)
        hTau1P1N.SetLineColor(ROOT.kBlue)
        hTau1P2N.SetLineColor(ROOT.kRed)
        hTau1PXN.SetLineColor(ROOT.kCyan)
        hTau3P0N.SetLineColor(ROOT.kMagenta)
        
        hTau1P0N.SetLineWidth(2)
        hTau1P1N.SetLineWidth(2)
        hTau1P2N.SetLineWidth(2)
        hTau1PXN.SetLineWidth(2)
        hTau3P0N.SetLineWidth(2)
    
        hTau1P0N.SetMarkerStyle(47)
        hTau1P1N.SetMarkerStyle(43)
        hTau1P2N.SetMarkerStyle(34)
        hTau1PXN.SetMarkerStyle(29)
        hTau3P0N.SetMarkerStyle(22)

        hTau1P0N.SetMarkerColor(ROOT.kGreen)
        hTau1P1N.SetMarkerColor(ROOT.kBlue)
        hTau1P2N.SetMarkerColor(ROOT.kRed)
        hTau1PXN.SetMarkerColor(ROOT.kCyan)
        hTau3P0N.SetMarkerColor(ROOT.kMagenta)
        
        hTau1P0N.SetMarkerSize(1.5)
        hTau1P1N.SetMarkerSize(1.5)
        hTau1P2N.SetMarkerSize(1.5)
        hTau1PXN.SetMarkerSize(1.5)
        hTau3P0N.SetMarkerSize(1.5)
        
        hTau1P0N.SetStats(0)
        hTau1P1N.SetStats(0)
        hTau1P2N.SetStats(0)
        hTau1PXN.SetStats(0)
        hTau3P0N.SetStats(0)
    
        c = ROOT.TCanvas('c', 'overlay', 800, 600)
        c.SetGridy()
        c.SetTicky()
    
        hTau1P0N.Draw()
        # hTau1P1N.Draw('SAME')
        # hTau1P2N.Draw('SAME')
        hTau1PXN.Draw('SAME')
        hTau3P0N.Draw('SAME')
    
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.04)
        text.SetTextAlign(13)
        text.DrawLatex(0.13, 0.87, "#bf{#it{MAIA}} Detector Concept")
        text.DrawLatex(0.13, 0.83, "Simulated Tau Gun (No BIB)")
        if i==0:
            text.DrawLatex(0.13, 0.79, "Central Barrel Region (1.0#kern[-0.6]{ }<#kern[-0.6]{ }#it{#theta}#kern[-0.6]{ }<#kern[-0.6]{ }2.0)")
        elif i==1:
            text.DrawLatex(0.13, 0.79, "Transition Region (0.577#kern[-0.8]{ }<#kern[-0.8]{ }#it{#theta}#kern[-0.8]{ }<#kern[-0.8]{ }1.0#kern[-0.8]{ }or#kern[-0.8]{ }2.0#kern[-0.8]{ }<#kern[-0.8]{ }#it{#theta}#kern[-0.8]{ }<#kern[-0.8]{ }2.56)")
        else:
            text.DrawLatex(0.13, 0.79, "Endcap Region (#it{#theta}#kern[-0.6]{ }<#kern[-0.6]{ }0.577#kern[-0.6]{ }or#kern[-0.6]{ }#it{#theta}#kern[-0.6]{ }>#kern[-0.6]{ }2.56)")
        c.Update()
        
        legend = ROOT.TLegend(0.75, 0.75, 0.92, 0.90)
        legend.AddEntry(hTau1P0N, '1P0N Tau Reco', 'lp')
        # legend.AddEntry(hTau1P1N, '1P1N Tau Reco', 'lp')
        # legend.AddEntry(hTau1P2N, '1P2N Tau Reco', 'lp')
        legend.AddEntry(hTau1PXN, '1PXN Tau Reco', 'lp')
        legend.AddEntry(hTau3P0N, '3P0N Tau Reco', 'lp')
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
    
        c.Update()
    
        filename = hTau1P0N.GetName() + '.png'
        c.SaveAs(filename)
