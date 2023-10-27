import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from result import plot

from fit.fit import fit_movie               
from result.plot import plot

import os
import yaml


def get_coords(Nx_loc, Ny_loc, coords):
    file_path ="./result/coords.dat"
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

def read_matrix_torch(file_path):
    # Read the data using pandas
    df = pd.read_csv(file_path, delim_whitespace=True, header=None)

    # Convert the DataFrame to a torch tensor
    data = torch.tensor(df.values, dtype=torch.float32)

    # Get the dimensions of the tensor
    Nx_loc, Ny_loc = data.size(1), data.size(0)

    return Nx_loc, Ny_loc, data

def write_matrix_torch(file_path, data):
    # Convert the tensor to a NumPy array
    data_np = data.numpy()

    # Write the NumPy array to a text file
    with open(file_path, 'w') as file:
        for row in data_np:
            # Convert each row to a space-separated string
            row_str = ' '.join(map(str, row))
            file.write(row_str + '\n')


def get_file_path(dir_name, header, time_index, mpi_index):

    '''example Hz_t=00870_06.dat  Ex_t=00900_42.dat '''
    file_path = f"{dir_name}/{header}_t={time_index:05d}_{mpi_index:02d}.dat"
    return file_path


def get_movie(time_index_list, mpi_index, coords_limit, header):
    nx_min, nx_max, ny_min, ny_max = coords_limit
    
    pics = []
    for time_index in time_index_list:
        _fpath = get_file_path('result/', header, time_index, mpi_index)
        Nx_loc, Ny_loc, data_torch = read_matrix_torch(_fpath)
        pics.append(data_torch[ny_min:ny_max,  nx_min:nx_max])   
    movie = torch.stack(pics, dim=-1)
    return movie 

if __name__=="__main__":
    # Load the YAML configuration file
    with open('config.yaml', 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    dx, c = config['dx'], config['c']

    dt = dx / (2.0 * c)
    omega = config['ev_to_radsec'] * config['omega_ev']
    
    Nx_loc, Ny_loc = config['Nx_loc'],config['Ny_loc']
    # get coords and nprocs
    Nx, Ny, coords =get_coords(Nx_loc, Ny_loc, {})

    # get time_index_list time_list
    t_sk, t_nm = config['time_skip'], config['time_num']
    time_index_list = [t_sk * (1 + i) for i in range(t_nm)]
    time_list = [index * dt for index in time_index_list ]
    
    #coords_limit = [0, 10, 0, 12]
    coords_limit = [0, Nx_loc, 0, Ny_loc]
    dir_name='result_fit'
    header='Ey'
    try:
        os.mkdir(dir_name)
    except:
        pass
    for header in ['Ey', 'Ex']:
        for mpi_index in range(48):
            print(f"mpi_index = {mpi_index}")
            movie =  get_movie(time_index_list, mpi_index, coords_limit, header)
            
            #fitter    
            A_fit, B_fit, C_fit =  fit_movie(time_list, movie, omega)
       
            file_path = f"{dir_name}/{header}_amplitude_{mpi_index:02d}.dat"
            write_matrix_torch(file_path, A_fit)

            file_path = f"{dir_name}/{header}_phase_{mpi_index:02d}.dat"
            write_matrix_torch(file_path, B_fit)
            '''
            # Plot the true and fitted data
            plt.figure(figsize=(10, 5))
            plt.scatter(time_list, movie.numpy()[0][0], label="True Data", alpha=0.5)
            plt.plot(time_list, C_fit.detach().numpy()[0][0], label="Fitted Data", color='red')
            plt.legend()
            plt.xlabel("Time")
            plt.ylabel("Value")
            plt.title("True vs. Fitted Data")
            plt.show()
            plt.scatter( torch.flatten(movie).numpy(),torch.flatten(C_fit).detach().numpy() ,label="prediction and truth", alpha=0.5)
            plt.show()
            #plot(A_fit.detach().numpy())
            #plot(B_fit.detach().numpy())
            '''

