import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
import math
from fit.test import create_test_movie


# Step 2: Define the model
class CosineModel(nn.Module):
    def __init__(self, n_rows, n_cols, omega):
        super(CosineModel, self).__init__()
        self.Amn = nn.Parameter(torch.abs(torch.randn(n_rows, n_cols)* 0.1))  # Learnable parameter for Amn
        self.Bmn = nn.Parameter(torch.rand(n_rows, n_cols) * 2 * math.pi)  # Learnable parameter for Bmn
        self.omega = omega

    def forward(self, t):
        # Calculate Cmnj using learned parameters Amn and Bmn
        Cmnj = self.Amn.unsqueeze(2) * torch.cos(self.omega * t.unsqueeze(0).unsqueeze(0) + self.Bmn.unsqueeze(2))
        return Cmnj


def fit_movie(time_list, movie, omega):
    n_rows, n_cols, _ = movie.shape
    model = CosineModel(n_rows=n_rows, n_cols=n_cols, omega=omega)

    # Step 3: Train the model
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-1) #5

    num_epochs = 500
    t = torch.tensor((time_list), dtype=torch.float32)

    for epoch in range(num_epochs):
        optimizer.zero_grad()
        movie_fit = model(t)
        loss = criterion(movie_fit, movie)
        loss.backward()
        optimizer.step()

    # Step 4: Verify the results
    A_fit = model.Amn.detach().to(torch.cfloat)
    B_fit = model.Bmn.detach().to(torch.cfloat)
    
    
    # A_fit_abs = torch.abs(A_fit)
    Z = A_fit * torch.exp(1j*B_fit)

    return Z.abs(), torch.remainder(Z.angle(), 2 * np.pi), movie_fit

if (__name__=='__main__'):
    # Step 1: Create random tensors A and B
    n_rows = 2  # Replace with your desired number of rows
    n_cols = 3  # Replace with your desired number of columns
    n_times = 53 # Replace with your desired number of time steps
    omega=2* np.pi
    noise_std = 0.055 # Adjust the noise level as needed

    A, B, C = create_test_movie(n_rows, n_cols, n_times, omega, noise_std)

    #t = torch.linspace(0, 1, n_times)  # Time values
    t = [i/(n_times+0.0) for i in range(n_times)]    

    A_fit, B_fit, C_fit = fit_movie(t, C, omega)
    print(f"True A: {A},\n Fitted A: {A_fit}\n\n")
    print(f"True B: {B},\n Fitted B: {B_fit}\n\n")
    print(f"True C: {C},\n Fitted C: {C_fit}\n\n")

    # Plot the true and fitted data
    plt.figure(figsize=(10, 5))
    plt.scatter(t, C.numpy()[0][0], label="True Data", alpha=0.5)
    plt.plot(t, C_fit.detach().numpy()[0][0], label="Fitted Data", color='red')
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Value")
    plt.title("True vs. Fitted Data")
    plt.show()

    plt.scatter( torch.flatten(C).numpy(),torch.flatten(C_fit).detach().numpy() ,label="prediction and truth", alpha=0.5)
    plt.show()
