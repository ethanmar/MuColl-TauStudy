from ROOT import TFile
import matplotlib.pyplot as plt

file = TFile.Open('./tau_gun/MAIA/TauFinderOutputs/Taus_loose.root')

tree = file.Get('anatree')

if not tree:
    print('TTree not found in file!')
    exit()

tau_pt = []
tau_isoE = []
tau_invM = []
isoE_100 = []
event_num = []
event_num_isoE_100 = []
nrej_isoE_tot = 0
n_isoE_more_5 = 0
n_invM_more_2 = 0


for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    nrej_isoE_tot += tree.nrej_isoE
    
    for j in range(tree.ntau):
        tau_pt.append(tree.t_pt[j])
        tau_isoE.append(tree.t_isoE[j])
        tau_invM.append(tree.t_minv[j])
        event_num.append(tree.event_num[j])
        
for isoE in tau_isoE:
    if isoE > 5:
        n_isoE_more_5 += 1

for invM in tau_invM:
    if invM > 2:
        n_invM_more_2 += 1
        
print(f'Number of taus that fail isoE: {nrej_isoE_tot}')
print(f'Fraction of taus with isoE > 5 GeV: {n_isoE_more_5/len(tau_isoE)}')
print(f'Fraction of taus with invM > 2 GeV/c^2: {n_invM_more_2/len(tau_invM)}')

for i in range(len(event_num)):
    if tau_isoE[i] > 100:
        event_num_isoE_100.append(event_num[i])
        isoE_100.append(tau_isoE[i])

print(f'Events with isoE > 100 GeV:')
for i in range(len(event_num_isoE_100)):
    print(f'Event number: {event_num_isoE_100[i]}     Isolation Energy: {isoE_100[i]}')
'''        
plt.scatter(tau_pt, tau_isoE, color='magenta', s=10)
plt.axhline(y=5, linestyle='--', color='black', label='5 GeV')
plt.xlabel(r'$p_T$ [GeV/c]', fontsize=12)
#plt.xlim(0,300)
#plt.ylim(0,30)
plt.ylabel(r'$E_{iso}$ [GeV]', fontsize=12)
plt.title('$\\tau^-$ Isolation Energy vs Transverse Momentum', fontsize=15)
plt.text(0.01, 0.99, "MAIA Detector Concept\nSimulated $\\tau^-$ Events",
         transform=plt.gca().transAxes, fontsize=12, va='top', ha='left')
plt.legend()
plt.savefig('pt_iso_E.png')

#plt.show()

'''
plt.scatter(tau_pt, tau_invM, color='magenta', s=10)
plt.axhline(y=2, linestyle='--', color='black', label=r'2 $GeV/c^2$')
plt.xlabel(r'$p_T$ [GeV/c]', fontsize=12)
#plt.xlim(0,300)
#plt.ylim(0,30)
plt.ylabel(r'$M_{inv}$ [$GeV/c^2$]', fontsize=12)
plt.title('$\\tau^-$ Invariant Mass vs Transverse Momentum', fontsize=15)
plt.text(0.01, 0.99, "MAIA Detector Concept\nSimulated $\\tau^-$ Events",
         transform=plt.gca().transAxes, fontsize=12, va='top', ha='left')
plt.legend()
plt.savefig('pt_inv_M.png')

#plt.show()


file.Close()
