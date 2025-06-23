from ROOT import TFile
import matplotlib.pyplot as plt

file = TFile.Open('./taufinder_cone_15_pt_1000.root')

tree = file.Get('anatree')

if not tree:
    print('TTree not found in file!')
    exit()

tau_pt = []
tau_isoE = []
tau_invM = []

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    
    for j in range(tree.ntau):
        tau_pt.append(tree.t_pt[j])
        tau_isoE.append(tree.t_isoE[j])
        tau_invM.append(tree.t_invMass[j])

plt.figure()
plt.scatter(tau_pt, tau_isoE, color='magenta', s=10)
# plt.axhline(y=5, linestyle='--', color='black', label='5 GeV')
plt.xlabel(r'$p_T^{reco}$ [GeV/c]', fontsize=12)
#plt.xlim(0,300)
#plt.ylim(0,30)
plt.ylabel(r'$E_{iso}$ [GeV]', fontsize=12)
plt.title('$\\tau^-$ Isolation Energy vs Transverse Momentum', fontsize=15)
# plt.text(0.01, 0.99, "MAIA Detector Concept\nSimulated $\\tau^-$ Events",
         # transform=plt.gca().transAxes, fontsize=12, va='top', ha='left')
# plt.legend()
plt.savefig('pt_iso_E.png')
plt.close()

plt.figure()
plt.scatter(tau_pt, tau_invM, color='magenta', s=10)
# plt.axhline(y=2, linestyle='--', color='black', label=r'2 $GeV/c^2$')
plt.xlabel(r'$p_T^{reco}$ [GeV/c]', fontsize=12)
#plt.xlim(0,300)
#plt.ylim(0,30)
plt.ylabel(r'$M_{inv}$ [$GeV/c^2$]', fontsize=12)
plt.title('$\\tau^-$ Invariant Mass vs Transverse Momentum', fontsize=15)
# plt.text(0.01, 0.99, "MAIA Detector Concept\nSimulated $\\tau^-$ Events",
         # transform=plt.gca().transAxes, fontsize=12, va='top', ha='left')
# plt.legend()
plt.savefig('pt_inv_M.png')
plt.close()

file.Close()
