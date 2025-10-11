import ROOT

variables = ['pt', 'energy', 'energyRatio']
units = ['[GeV/c]', '[GeV]', '']
names = [' Reconstructed#kern[-0.5]{ }#it{p_{T}}', '#it{E}^{reco}', '#it{E}^{reco}/#it{E}_{#tau}^{reco}']

root_file = ROOT.TFile('pfo_ana.root', 'READ')

n_taus = 46268

for i, variable in enumerate(variables):
    hPhoton = root_file.Get(variable + 'PFOs_pdg22')
    hPion = root_file.Get(variable + 'PFOs_pdg211')
    hNeutron = root_file.Get(variable + 'PFOs_pdg2112')

    hPhoton.Scale(1/n_taus)
    hPion.Scale(1/n_taus)
    hNeutron.Scale(1/n_taus)
    
    hNeutron.SetTitle('')
    hNeutron.SetName('pfo_' + variable)
    hNeutron.GetYaxis().SetTitle('Number Reconstructed per Tau')
    hNeutron.GetXaxis().SetTitle(names[i] + '#kern[-0.5]{ }' + units[i])
    hNeutron.GetXaxis().SetRangeUser(0.0, 5.0)

    hPhoton.SetLineColor(ROOT.kCyan)
    hPion.SetLineColor(ROOT.kMagenta)
    hNeutron.SetLineColor(ROOT.kGreen)

    hPhoton.SetLineWidth(2)
    hPion.SetLineWidth(2)
    hNeutron.SetLineWidth(2)

    hPhoton.SetStats(0)
    hPion.SetStats(0)
    hNeutron.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)
    c.SetLeftMargin(0.15)
    
    hNeutron.Draw()
    hPion.Draw('SAME')
    hPhoton.Draw('SAME')

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.SetTextAlign(13)
    text.DrawLatex(0.22, 0.87, "#bf{#it{MAIA}} Detector Concept")
    text.DrawLatex(0.22, 0.83, "Simulated Tau Gun (No BIB)")
    
    c.Update()
        
    legend = ROOT.TLegend(0.75, 0.75, 0.90, 0.90)
    legend.AddEntry(hPhoton, 'Photon', 'l')
    legend.AddEntry(hPion, 'Pion', 'l')
    legend.AddEntry(hNeutron, 'Neutron', 'l')
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    
    filename = hNeutron.GetName() + '.png'
    c.SaveAs(filename)
