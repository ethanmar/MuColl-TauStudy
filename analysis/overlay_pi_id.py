import ROOT

names = ['centroid_r', 'e_over_e_plus_h']
labels = ['Centroid R [mm]', 'E_{ECAL}/E_{TOT}']


root_file = ROOT.TFile('pi_id.root', 'READ')

n_pis = 37771 #23542
n_elecs = 21606 #4178

for i, name in enumerate(names):
    hPion = root_file.Get('pi_' + name)
    hElectron = root_file.Get('elec_' + name)

    if (i == 0):
        hPion = hPion.Rebin(4)
        hElectron = hElectron.Rebin(4)

    hPion.Scale(1/n_pis)
    hElectron.Scale(1/n_elecs)
    
    hPion.SetTitle('')
    hPion.SetName(name)
    hPion.GetYaxis().SetTitle('Normalized Count')
    hPion.GetXaxis().SetTitle(labels[i])
    if (i == 0):
        hPion.GetYaxis().SetRangeUser(0, hElectron.GetMaximum()+0.001)
    else:
        hPion.GetYaxis().SetRangeUser(0, hElectron.GetMaximum()+0.01)

    hPion.SetLineColor(ROOT.kMagenta)
    hElectron.SetLineColor(ROOT.kGreen)

    hPion.SetLineWidth(2)
    hElectron.SetLineWidth(2)

    hPion.SetStats(0)
    hElectron.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)
    c.SetLeftMargin(0.15)
    c.SetGridy()
    
    hPion.Draw('HIST')
    hElectron.Draw('HIST SAME')
        
    legend = ROOT.TLegend(0.75, 0.75, 0.90, 0.90)
    legend.AddEntry(hPion, 'Matched Pion', 'l')
    legend.AddEntry(hElectron, 'Matched Electron', 'l')
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    
    filename = hPion.GetName() + '.png'
    c.SaveAs(filename)
