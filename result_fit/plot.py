import torch
import matplotlib.pyplot as plt
import numpy as np

import yaml
with open('../config.yaml', 'r') as config_file:
    config = yaml.safe_load(config_file)

def plot(data_array, cmap='hot'):
    Ny = data_array.shape[0]
    Nx = data_array.shape[1]
    # Create a colormap plot
    plt.imshow(np.reshape(data_array,(Ny,Nx)), cmap=cmap, aspect='auto')  # 'viridis' is just one of many available colormaps
    plt.colorbar()  # Add a colorbar to the plot

    # You can add labels, titles, or customize the plot further as neede:
    plt.xlabel('X Axis Label')
    plt.ylabel('Y Axis Label')
    plt.title('Colormap Plot')

    # Show the plot
    plt.show()
def parse(x):
    return float(x)

def get_coords(Nx_loc, Ny_loc, coords):
    file_path ="coords.dat"
    x=0
    y=0
    # Open the file for reading
    with open(file_path, 'r') as file:
        # Read each line from the file
        for line in file:
            # Split the line into individual values (assuming values are space-separated)
            values = line.strip().split()
            if(len(values)!=3):
                continue
            coords[int(values[0])] = (int(values[1]), int(values[2]))
            x = max(x, int(values[1]))
            y = max(y, int(values[2]))
    return (x+1)*Nx_loc, (y+1)*Ny_loc, coords
def read_matrix( file_path):
    # Define the file path

    # Initialize an empty nested list to store the data
    data = []
    Nx_loc, Ny_loc = 0, 0
    # Open the file for reading
    with open(file_path, 'r') as file:
        # Read each line from the file
        for line in file:
            # Split the line into individual values (assuming values are space-separated)
            values = line.strip().split()
            Nx_loc = len(values) 
        # Convert values to float and append to the nested list
            for value in values:
                data.append(parse(value))
            Ny_loc += 1
    return Nx_loc, Ny_loc, data

def get_all_matrix(name):
    Nx_loc, Ny_loc = 0,0
    arrs= []
    for X in range(48):
        if( X<10):
            _X = "0" + str(X)
        else:
            _X = str(X)
        Nx_loc, Ny_loc, arr=read_matrix(f'{name}_{_X}.dat')
        #print('loc',Nx_loc, Ny_loc)
        #plt.imshow(np.reshape(arr,(Ny_loc,Nx_loc)), cmap='viridis', aspect='auto')  # 'viridis' is just one of many available colormaps
        #plt.show()
        arrs.append(arr)
    #plot(np.array(arrs[10]).reshape(Ny_loc, Nx_loc))

    Nx, Ny, coords =get_coords(Nx_loc, Ny_loc, {})
    data_array = np.array([0.0] * Ny * Nx).reshape((Ny,Nx))
    for i, arr in enumerate(arrs):
        bx, by = coords[i]
        data_array[by* Ny_loc: (by+1)*Ny_loc ,bx* Nx_loc: (bx+1)*Nx_loc] = np.array(arr).reshape((Ny_loc, Nx_loc))

    plot(data_array)
    return data_array

def get_pos(ri):
    xi, yi = ri
    dx, dy, Nx, Ny = config['dx'], config['dy'], config['Nx'], config['Ny']
    x = -dx * Nx/2.0 + dx * xi
    y = -dy * Ny/2.0 + dy * yi
    return (x, y)

def get_index(r):
    x, y = r
    dx, dy, Nx, Ny = config['dx'], config['dy'], config['Nx'], config['Ny']
    xi = int((x + dx * Nx/2.0 ) / dx)
    yi = int((y + dy * Ny/2.0 ) / dy)
    return (xi, yi)

def get_efield_inc(r, t):
    x, y = r
    E0 = config['E0']
    omega = config['omega']
    ms, js = config['ms'], config['js']
    dt, c = config['dt'], config['c']
    omega = config['ev_to_radsec'] * config['omega_ev']
    t_measure = config['Nt_print'] * dt
    k = omega / c
    y_inc = (ms* config['Ny_loc'] + js)* config['dy']
    E = E0* np.cos(omega * (t + t_measure) + k* (y-y_inc) -np.pi/2) 
    return E

def get_efield_phase(r):
    x, y = r
    E0 = config['E0']
    ms, js = config['ms'], config['js']
    dt, c = config['dt'], config['c']
    omega = config['ev_to_radsec'] * config['omega_ev']
    t_measure = config['Nt_print'] * dt
    k = omega / c
    y_inc = (ms* config['Ny_loc'] + js)* config['dy']
#   E = E0* cos(omega * (t + t_measure) + k* (y-y_inc) -np.pi/2) 
    phase = omega* t_measure +k* (y-y_inc) -np.pi/2
    phase =  np.remainder(phase, np.pi * 2)
    return phase

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
    plot(torch.abs(E_loc_x), cmap='hot')
    plot(torch.angle(E_loc_x), cmap='bwr')
    plot(Ey_ph-np.pi*torch.ones(Ex_amp.shape), cmap='bwr')

    #plot(phases)
    #get_all_matrix('Hz')

    #Et = Ex.copy()
    #for j in range((Et.shape[0])):
    #    for i in range((Et.shape[1])):
    #        Et[j,i] = (Et[j, i] **2 + Ey[j, i] **2 )**0.5
    #plot(Et)
