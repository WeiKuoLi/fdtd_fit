import torch
import random
import numpy as np

q=3
def create_test_movie(n_rows=2, n_cols=3, n_times=10, omega=2, noise_std=0.1):
    
    A = torch.randn(n_rows, n_cols)  # Gaussian sampled tensor centered at 0
    B = torch.rand(n_rows, n_cols) * 2 * np.pi  # Random values between 0 and 2*pi

    # Step 2: Create a time_list
    time_list = torch.linspace(0, 1, n_times)  # Adjust the time range as needed

    # Step 3: Calculate 3D tensor C
    C = torch.zeros(n_rows, n_cols, n_times)
    '''
    C_brute = torch.zeros(n_rows, n_cols, n_times)

    #for i in range(n_rows):
        for j in range(n_cols):
            for k in range(n_times):
                C_brute[i, j, k] = A[i, j] * np.cos(time_list[k] *omega + B[i, j])
    '''
    # Step 3: Calculate 3D tensor C using broadcasting
    A_expanded = A.unsqueeze(2)  # Add a third dimension to A
    B_expanded = B.unsqueeze(2)  # Add a third dimension to B
    time_list_expanded = time_list.unsqueeze(0).unsqueeze(0)  # Add dimensions for broadcasting

    C = A_expanded * torch.cos(time_list_expanded *omega+ B_expanded)

    # Step 4: Add noise to C
    noise = torch.randn(n_rows, n_cols, n_times) * noise_std
    C_noise = C + noise
    
    return A, B, C_noise



if __name__=="__main__":
    # Set the seed for reproducibility
    torch.manual_seed(42)

    # Step 1: Create random tensors A and B
    n_rows = 2  # Replace with your desired number of rows
    n_cols = 3  # Replace with your desired number of columns
    n_times = 2  # Replace with your desired number of time steps
    omega=2
    noise_std = 0.0  # Adjust the noise level as needed

    A, B, C = create_test_movie(n_rows, n_cols, n_times, omega, noise_std)
    # Print the tensors A, B, C, and C_noise
    print("A:")
    print(A)

    print("B:")
    print(B)

    print("C:")
    print(C)
