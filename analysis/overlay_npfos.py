import ROOT

prongs = ['1', '3']

for prong in prongs:
    root_file_25 = ROOT.TFile('./tau_gun/MAIA/TauFinderOutputs/' + prong + 'P0N_ana_cone_25_pt_00.root', 'READ')
    root_file_15 = ROOT.TFile('./tau_gun/MAIA/TauFinderOutputs/' + prong + 'P0N_ana_cone_15_pt_1000.root', 'READ')

    h15 = root_file_15.Get('nPFOs_' + prong + 'P0N')
    h25 = root_file_25.Get('nPFOs_' + prong + 'P0N')

    h15.SetTitle('')
    h15.SetName('npfos_' + prong + 'P0N_comparison')
    h15.GetXaxis().SetTitle('Number of Reconstructed PFOs')

    h15.SetLineColor(ROOT.kBlue)
    h25.SetLineColor(ROOT.kRed)

    h15.SetLineWidth(2)
    h25.SetLineWidth(2)

    h15.SetStats(0)
    h25.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)

    h15.Draw()
    h25.Draw('SAME')

    legend = ROOT.TLegend(0.45, 0.70, 0.85, 0.90)
    legend.AddEntry(h25, '#theta_{cone} = 0.25 rad, #it{p}_{T}^{pfos} > 0.0 GeV/c (' + prong + 'P0N)', 'l')
    legend.AddEntry(h15, '#theta_{cone} = 0.15 rad, #it{p}_{T}^{pfos} > 1.0 GeV/c (' + prong + 'P0N)', 'l')
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()

    filename = h15.GetName() + '.png'
    c.SaveAs(filename)
