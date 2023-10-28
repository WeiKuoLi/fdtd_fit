from func.get import get_all_matrix, config, get_efield_phase, get_index
import numpy as np
import torch
import sys
sys.path.append("../plasmon-polariton")
from model import oscillator, sample_atoms 
from lab.equipment import InfoList as IL
from function.func import find_gc_info


def donut(r1, r2, npt):
    pt = torch.sqrt(torch.rand([npt]) * (r2**2-r1**2) + r1**2)
    theta = torch.rand([npt]) * 2 * np.pi
    return pt*torch.cos(theta), pt*torch.sin(theta)

def random_unit(dim=3, num=100):
    vecs = torch.randn(num, dim)
    vecs_len = torch.unsqueeze(torch.sqrt(torch.sum(vecs * vecs,1)),1)
    return torch.div(vecs, vecs_len)

def sample_donut(list_osc, dp,  r1, r2):
    ''' given a list of oscs, introduce noise to the gc except first one. gc = gc0 *cos(x) , x sampled from gaussian'''
    npt = len(list_osc) - 1
    X, Y = donut(r1, r2, npt)
    vecs = random_unit(3, npt)
    ev_to_J = config['ev_to_J']
    ap = 0.1 #???
    ####################
    Ey_amp = torch.tensor(get_all_matrix('Ey_amplitude'))
    Ey_ph = torch.tensor(get_all_matrix('Ey_phase'))
    Ex_amp = torch.tensor(get_all_matrix('Ex_amplitude'),dtype=torch.cfloat)
    Ex_ph = torch.tensor(get_all_matrix('Ex_phase'),dtype=torch.cfloat)

    
    #plot(Ex_amp*np.cos(Ex_ph))   
    
    y = [float(i) * config['dy'] for i in range(config['Ny'])]
    phase = [get_efield_phase((0.0,yi)) for yi in y]
    phases = torch.tensor(np.array([[p]* config['Nx'] for p in phase]), dtype=torch.cfloat)
    
    E0 = config['E0'] * torch.ones(Ex_amp.shape, dtype=torch.cfloat)
    
    E_loc_x = Ex_amp*torch.exp(1j*Ex_ph)- E0 * torch.exp(1j*phases)






    list_osc_samp = IL([list_osc[0]])
    for i in range(npt):
        _x, _y = X[i], Y[i]
        _i, _j = get_index((_x, _y))
        
        _E_0_x = (E0 * torch.exp(1j*phases))[_i, _j]        
        _E_loc_x  = E_loc_x[_i, _j]
        _E_loc_y = (Ey_amp*torch.exp(1j*Ey_ph))[_i, _j]
        gv = dp * vecs[i][0] * _E_0_x / config['E0'] / ev_to_J
        gc = dp *(vecs[i][0] * _E_loc_x + vecs[i][1] * _E_loc_y) / ap / ev_to_J
        osc = list_osc[i]
        list_osc_samp.append(oscillator(osc.wv, osc.gmv, gv, gc))
    gc_av, gc_sdv = find_gc_info(list_osc_samp)
    list_osc_samp.info = {'gc_av_sample':gc_av, 'gc_sdv_sample':gc_sdv}
    return list_osc_samp



if (__name__=='__main__'):
    '''
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
    '''
    mol = oscillator(2.7, .3, 6)
    mol2 = oscillator(2.7, .2, .3, .54)
    list_mol = [mol2 for i in range(1000)]
    list_mol[0] = mol
    dp = config['dp']
    r1 = config['radius']
    r2 = r1*1.1
    list_mol_samp = sample_donut(list_mol, dp,r1,r2)
    print(list_mol_samp)

