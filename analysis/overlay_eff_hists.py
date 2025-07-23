import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

# parser.add_argument('--inputFileOld', type=str, default='./tau_gun/MAIA/TauFinderOutputs/tau_ana_cone_10_pt_01.root')
parser.add_argument('--inputFileTau', type=str, default='tau_eff.root')
parser.add_argument('--inputFilePi', type=str, default='pi_eff.root')
args = parser.parse_args()

# root_file_old = ROOT.TFile(args.inputFileOld, 'READ')
root_file_tau = ROOT.TFile(args.inputFileTau, 'READ')
root_file_pi = ROOT.TFile(args.inputFilePi, 'READ')



nums = ['1', '3']
names = ['pt', 'phi', 'theta']
variables = ['#it{p}_{T}^{vis}', '#phi^{vis}', '#theta^{vis}']
units = ['GeV/c', 'rad', 'rad']
# regions = ['_barrel', '_centbarrel']

for num in nums:
    for i in range(3):
        # for j in range(2):
        # hOld = root_file_old.Get('reco_' + num + 'p_' + names[i] + '_eff' + regions[j])
        hTau = root_file_tau.Get(num + 'p0n_tau_' + names[i] + '_eff')
        hPion = root_file_pi.Get(num + 'p0n_pi_' + names[i] + '_eff')
        
        hPion.SetTitle('')
        hPion.SetName(num + 'p0n_' + names[i] + '_eff')
        hPion.GetXaxis().SetTitle('Tau Truth ' + variables[i] + ' [' + units[i] + ']')
        hPion.GetYaxis().SetRangeUser(0.0, 1.3)
            
        # hOld.SetLineColor(ROOT.kRed)
        hTau.SetLineColor(ROOT.kBlue)
        hPion.SetLineColor(ROOT.kGreen)
            
        # hOld.SetLineWidth(2)
        hTau.SetLineWidth(2)
        hPion.SetLineWidth(2)
            
        # hOld.SetMarkerStyle(34)
        hTau.SetMarkerStyle(43)
        hPion.SetMarkerStyle(47)

        # hOld.SetMarkerColor(ROOT.kRed)
        hTau.SetMarkerColor(ROOT.kBlue)
        hPion.SetMarkerColor(ROOT.kGreen)

        # hOld.SetMarkerSize(1.5)
        hTau.SetMarkerSize(1.5)
        hPion.SetMarkerSize(1.5)
        
        # hOld.SetStats(0)
        hTau.SetStats(0)
        hPion.SetStats(0)

        c = ROOT.TCanvas('c', 'overlay', 800, 600)

        # hOld.Draw()
        hPion.Draw()
        hTau.Draw('SAME')

        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.04)
        text.SetTextAlign(13)
        text.DrawLatex(0.15, 0.87, "#it{MAIA} #bf{Work in Progress}")
        text.DrawLatex(0.15, 0.83, "#bf{Charged Tau Gun (No BIB)}")
        # if regions[j] == '_centbarrel':
            # text.DrawLatex(0.15, 0.79, "#bf{Central Barrel Region (1.0 < #theta < 2.0 rad)}")
        # if regions[j] == '_barrel':
            # text.DrawLatex(0.15, 0.79, "#bf{Barrel Region (0.70 < #theta < 2.45 rad)}")
        c.Update()

        legend = ROOT.TLegend(0.70, 0.75, 0.90, 0.90)
        legend.AddEntry(hPion, num + 'P0N Pion Reco', 'lp')
        legend.AddEntry(hTau, num + 'P0N Tau Reco', 'lp')
        # legend.AddEntry(hOld, num + 'P0N Tau Reco (Old)', 'lp')
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        c.Update()

        filename = hPion.GetName() + '.png'
        c.SaveAs(filename)
