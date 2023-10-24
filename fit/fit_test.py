import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Generate synthetic data
# Let's assume A = 2, w = 2*pi, phi = pi/4, and we have time values from 0 to 10.

A_true = .20
w_true = .4 * np.pi
phi_true = np.pi / 4
t = torch.linspace(0, 10, 1000)  # Time values
y_true = A_true * torch.cos(w_true * t + phi_true) + torch.randn_like(t) * 0.03  # Add some noise

# Step 2: Define the model
class CosineModel(nn.Module):
    def __init__(self):
        super(CosineModel, self).__init__()
        self.A = nn.Parameter(torch.rand(1, requires_grad=True)) 
        self.w = nn.Parameter((torch.rand(1, requires_grad=True)*0.2 + 0.3) * np.pi)
        self.phi = nn.Parameter(torch.rand(1, requires_grad=True) * np.pi)

    def forward(self, x):
        return self.A * torch.cos(self.w * x + self.phi)

model = CosineModel()

# Step 3: Train the model
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.01)

num_epochs = 1000
for epoch in range(num_epochs):
    optimizer.zero_grad()
    y_pred = model(t)
    loss = criterion(y_pred, y_true)
    loss.backward()
    optimizer.step()

# Step 4: Verify the results
A_fit = model.A.item()
w_fit = model.w.item()
phi_fit = model.phi.item()

print(f"True A: {A_true}, Fitted A: {A_fit}")
while(phi_fit>0):
    phi_fit -= np.pi
phi_fit +=np.pi
print(f"True w: {w_true}, Fitted w: {w_fit}")
print(f"True phi: {phi_true}, Fitted phi: {phi_fit}")

# Plot the true and fitted data
plt.figure(figsize=(10, 5))
plt.scatter(t.numpy(), y_true.numpy(), label="True Data", alpha=0.5)
plt.plot(t.numpy(), y_pred.detach().numpy(), label="Fitted Data", color='red')
plt.legend()
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("True vs. Fitted Data")
plt.show()

