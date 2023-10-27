from func.get import get_all_matrix, config, get_efield_phase
import numpy as np
import torch
import sys
sys.path.append("../plasmon-polariton")
from model import oscillator, sample_atoms 



def donut(r1, r2, npt):
    pt = torch.sqrt(torch.rand([npt]) * (r2**2-r1**2) + r1**2)
    theta = torch.rand([npt]) * 2 * np.pi
    return pt*torch.cos(theta), pt*torch.sin(theta)

def sample_donut(list_osc, r1, r2):
    ''' given a list of oscs, introduce noise to the gc except first one. gc = gc0 *cos(x) , x sampled from gaussian'''
    npt = len(list_osc)
    X, Y = donut(r1, r2, npt)
    ####################

    thetas = np.random.normal(mean, std_dev, N-1).tolist()
    cos_thetas = np.cos(thetas)
    list_osc_samp = IL([list_osc[0]])
    for i in range(1,N):
        osc = list_osc[i]
        _cos_thetas = cos_thetas[i-1]
        list_osc_samp.append(oscillator(osc.wv, osc.gmv, osc.gv, osc.gc * _cos_thetas))
    gc_av, gc_sdv = find_gc_info(list_osc_samp)
    list_osc_samp.info = {'gc_av_sample':gc_av, 'gc_sdv_sample':gc_sdv}
    return list_osc_samp



if (__name__=='__main__'):
    #get_all_matrix('metal')
    Ey_amp = get_all_matrix('Ey_amplitude')
    print(torch.tensor(Ey_amp).shape)
    Ey_ph = torch.tensor(get_all_matrix('Ey_phase'))
    Ex_amp = torch.tensor(get_all_matrix('Ex_amplitude'),dtype=torch.cfloat)
    Ex_ph = torch.tensor(get_all_matrix('Ex_phase'),dtype=torch.cfloat)

    
    #plot(Ex_amp*np.cos(Ex_ph))   
    
    y = [float(i) * config['dy'] for i in range(config['Ny'])]
    phase = [get_efield_phase((0.0,yi)) for yi in y]
    phases = torch.tensor(np.array([[p]* config['Nx'] for p in phase]), dtype=torch.cfloat)
    
    E0 = config['E0'] * torch.ones(Ex_amp.shape, dtype=torch.cfloat)
    
    E_loc_x = Ex_amp*torch.exp(1j*Ex_ph)- E0 * torch.exp(1j*phases)
    print('E_local')
    print(E_loc_x, Ey_amp, Ey_ph)
   #plot(torch.abs(E_loc_x), cmap='hot')
   #plot(torch.angle(E_loc_x), cmap='bwr')
   #plot(Ey_ph-np.pi*torch.ones(Ex_amp.shape), cmap='bwr')

    #plot(phases)
    #get_all_matrix('Hz')
    mol = oscillator(20, 2, 6)
    mol2 = oscillator(20, 2, .3, .54)
    list_mol = [mol2 for i in range(5)]
    list_mol[0] = mol
    distrib = (0, 0.1*np.pi) 
    list_mol_samp = sample_atoms(list_mol, distrib)
    print(list_mol_samp)

