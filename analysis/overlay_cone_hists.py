import ROOT

angles = ['05', '08', '10', '12', '15', '18', '20', '25']

root_files = []

for angle in angles:
    root_files.append(ROOT.TFile('./tau_gun/MAIA/TauFinderOutputs/tau_ana_loose_cone_' + angle + '.root', 'READ'))

nums = ['1', '3']
names = ['pt', 'phi', 'theta']
names_denom = ['pT', 'phi', 'theta']
variables = ['p_{T}', '#phi', '#theta']
units = ['GeV/c', 'rad', 'rad']

for num in nums:
    for i in range(3):

        hists = []
        hists_denom = []
        
        h05 = root_files[0].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h05)
        h08 = root_files[1].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h08)
        h10 = root_files[2].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h10)
        h12 = root_files[3].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h12)
        h15 = root_files[4].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h15)
        h18 = root_files[5].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h18)
        h20 = root_files[6].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h20)
        h25 = root_files[7].Get('reco_' + num + 'p_' + names[i] + '_eff')
        hists.append(h25)
        hPi = root_files[0].Get(num + 'p_pi_' + names[i] + '_eff')

        if num == '1':
            h05_denom = root_files[0].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h05_denom)
            h08_denom = root_files[1].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h08_denom)
            h10_denom = root_files[2].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h10_denom)
            h12_denom = root_files[3].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h12_denom)
            h15_denom = root_files[4].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h15_denom)
            h18_denom = root_files[5].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h18_denom)
            h20_denom = root_files[6].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h20_denom)
            h25_denom = root_files[7].Get(num + 'p_true_vis_' + names_denom[i])
            hists_denom.append(h25_denom)
        else:
            if names_denom[i] == 'phi':
                h05_denom = root_files[0].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h05_denom)
                h08_denom = root_files[1].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h08_denom)
                h10_denom = root_files[2].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h10_denom)
                h12_denom = root_files[3].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h12_denom)
                h15_denom = root_files[4].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h15_denom)
                h18_denom = root_files[5].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h18_denom)
                h20_denom = root_files[6].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h20_denom)
                h25_denom = root_files[7].Get(num + 'p_tau_' + names_denom[i] + '_true')
                hists_denom.append(h25_denom)
            else:
                h05_denom = root_files[0].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h05_denom)
                h08_denom = root_files[1].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h08_denom)
                h10_denom = root_files[2].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h10_denom)
                h12_denom = root_files[3].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h12_denom)
                h15_denom = root_files[4].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h15_denom)
                h18_denom = root_files[5].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h18_denom)
                h20_denom = root_files[6].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h20_denom)
                h25_denom = root_files[7].Get(num + 'p_tau_true_' + names_denom[i])
                hists_denom.append(h25_denom)

        angles_final = ['0.05', '0.08', '0.10', '0.12', '0.15', '0.18', '0.20', '0.25']
                
        for j in range(len(angles_final)):
            nbins = hists[j].GetNbinsX()
            weighted_sum = 0.0
            total_weight = 0.0
            variance_sum = 0.0

            for k in range(1, nbins + 1):
                eff = hists[j].GetBinContent(k)
                weight = hists_denom[j].GetBinContent(k)
                err = hists[j].GetBinError(k)
                weighted_sum += eff * weight
                total_weight += weight
                variance_sum += (weight**2) * (err**2)

            avg_eff = weighted_sum/total_weight if total_weight > 0 else 0
            avg_err = (variance_sum/(total_weight**2))**0.5 if total_weight > 0 else 0
                
            print(f'Mean Efficiency for {angles_final[j]} rad Search Cone ({names[i]}): {avg_eff} +/- {avg_err}')

        h05.SetTitle(num + '-Prong Reconstruction Efficiencies vs ' + variables[i])
        h05.GetXaxis().SetTitle('True Visible #tau^{-} ' + variables[i] + ' [' + units[i] + ']')
        if num == '1':
            h05.GetYaxis().SetRangeUser(0.15, 1.2)
        else:
            h05.GetYaxis().SetRangeUser(0.0, 1.0)
        
        h05.SetLineColor(ROOT.kRed)
        h08.SetLineColor(ROOT.kBlue)
        h10.SetLineColor(ROOT.kCyan)
        h12.SetLineColor(ROOT.kCyan)
        h15.SetLineColor(ROOT.kMagenta)
        h18.SetLineColor(ROOT.kOrange)
        h20.SetLineColor(ROOT.kOrange)
        h25.SetLineColor(ROOT.kViolet)
        hPi.SetLineColor(ROOT.kGreen)

        h05.SetLineWidth(2)
        h08.SetLineWidth(2)
        h10.SetLineWidth(2)
        h12.SetLineWidth(2)
        h15.SetLineWidth(2)
        h18.SetLineWidth(2)
        h20.SetLineWidth(2)
        h25.SetLineWidth(2)
        hPi.SetLineWidth(2)

        h05.SetMarkerStyle(8)
        h08.SetMarkerStyle(34)
        h10.SetMarkerStyle(33)
        h12.SetMarkerStyle(29)
        h15.SetMarkerStyle(23)
        h18.SetMarkerStyle(22)
        h20.SetMarkerStyle(21)
        h25.SetMarkerStyle(47)
        hPi.SetMarkerStyle(29)

        h05.SetMarkerColor(ROOT.kRed)
        h08.SetMarkerColor(ROOT.kBlue)
        h10.SetMarkerColor(ROOT.kCyan)
        h12.SetMarkerColor(ROOT.kCyan)
        h15.SetMarkerColor(ROOT.kMagenta)
        h18.SetMarkerColor(ROOT.kOrange)
        h20.SetMarkerColor(ROOT.kOrange)
        h25.SetMarkerColor(ROOT.kViolet)
        hPi.SetMarkerColor(ROOT.kGreen)
        
        h05.SetMarkerSize(1.5)
        h08.SetMarkerSize(1.5)
        h10.SetMarkerSize(1.5)
        h12.SetMarkerSize(1.5)
        h15.SetMarkerSize(1.5)
        h18.SetMarkerSize(1.5)
        h20.SetMarkerSize(1.5)
        h25.SetMarkerSize(1.5)
        hPi.SetMarkerSize(1.5)

        h05.SetStats(0)
        h08.SetStats(0)
        h10.SetStats(0)
        h12.SetStats(0)
        h15.SetStats(0)
        h18.SetStats(0)
        h20.SetStats(0)
        h25.SetStats(0)
        hPi.SetStats(0)

        c = ROOT.TCanvas('c', 'overlay', 800, 600)

        h05.Draw()
        # h08.Draw('SAME')
        h10.Draw('SAME')
        # h12.Draw('SAME')
        h15.Draw('SAME')
        # h18.Draw('SAME')
        h20.Draw('SAME')
        h25.Draw('SAME')
        hPi.Draw('SAME')
        
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.04)
        text.SetTextAlign(13)
        text.DrawLatex(0.12, 0.89, "#bf{#it{MAIA Detector Concept}}")
        text.DrawLatex(0.12, 0.84, "Simulated #tau^{-} Events")

        c.Update()

        legend = ROOT.TLegend(0.75, 0.70, 0.90, 0.90)
        legend.AddEntry(h05, '0.05 rad', 'lp')
        # legend.AddEntry(h08, '0.08 rad', 'lp')
        legend.AddEntry(h10, '0.10 rad', 'lp')
        # legend.AddEntry(h12, '0.12 rad', 'lp')
        legend.AddEntry(h15, '0.15 rad', 'lp')
        # legend.AddEntry(h18, '0.18 rad', 'lp')
        legend.AddEntry(h20, '0.20 rad', 'lp')
        legend.AddEntry(h25, '0.25 rad', 'lp')
        legend.AddEntry(hPi, 'Pion Eff', 'lp')
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        c.Update()

        filename = h05.GetTitle() + '.png'
        c.SaveAs(filename)
