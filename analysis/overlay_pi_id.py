import ROOT

names = ['centroid_r', 'e_over_e_plus_h']
labels = ['Centroid R [mm]', 'E_{ECAL}/E_{TOT}']
regions = ['cent_barrel', 'transition', 'endcap']


root_file = ROOT.TFile('pi_id.root', 'READ')

for i, region in enumerate(regions):
    for j, name in enumerate(names):
        hPion = root_file.Get('pi_' + name + '_' + region)
        hElectron = root_file.Get('elec_' + name + '_' + region)

        if (j == 0):
            hPion = hPion.Rebin(4)
            hElectron = hElectron.Rebin(4)

        hPion.Scale(1/hPion.GetEntries())
        hElectron.Scale(1/hElectron.GetEntries())
        
        hPion.SetTitle('')
        hPion.SetName(name + '_' + region)
        hPion.GetYaxis().SetTitle('Normalized Count')
        hPion.GetXaxis().SetTitle(labels[j])
        y_max = max((hPion.GetMaximum(), hElectron.GetMaximum()))
        hPion.GetYaxis().SetRangeUser(0, y_max + y_max/4)
        '''if (j == 0):
            hPion.GetYaxis().SetRangeUser(0, 0.045)
        else:
            hPion.GetYaxis().SetRangeUser(0, 0.90)'''

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

        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.03)
        text.SetTextAlign(13)
        text.DrawLatex(0.18, 0.87, "#bf{#it{MAIA}} Detector Concept")
        text.DrawLatex(0.18, 0.83, "Simulated Tau Gun (No BIB)")
        if i==0:
            text.DrawLatex(0.18, 0.79, "Central Barrel Region (1.0#kern[-0.6]{ }<#kern[-0.6]{ }#it{#theta}#kern[-0.6]{ }<#kern[-0.6]{ }2.0)")
        elif i==1:
            text.DrawLatex(0.18, 0.79, "Transition Region (0.577#kern[-0.8]{ }<#kern[-0.8]{ }#it{#theta}#kern[-0.8]{ }<#kern[-0.8]{ }1.0#kern[-0.8]{ }or#kern[-0.8]{ }2.0#kern[-0.8]{ }<#kern[-0.8]{ }#it{#theta}#kern[-0.8]{ }<#kern[-0.8]{ }2.56)")
        else:
            text.DrawLatex(0.18, 0.79, "Endcap Region (#it{#theta}#kern[-0.6]{ }<#kern[-0.6]{ }0.577#kern[-0.6]{ }or#kern[-0.6]{ }#it{#theta}#kern[-0.6]{ }>#kern[-0.6]{ }2.56)")
        c.Update()
        
        legend = ROOT.TLegend(0.75, 0.75, 0.90, 0.90)
        legend.AddEntry(hPion, 'Matched Pion', 'l')
        legend.AddEntry(hElectron, 'Matched Electron', 'l')
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()

        c.Update()
        
        filename = hPion.GetName() + '.png'
        c.SaveAs(filename)
