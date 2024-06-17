import numpy as np
import matplotlib.pyplot as plt

rho = np.load('rho.npy')
u = np.load('u.npy')
v = np.load('v.npy')
p = np.load('p.npy')

x = np.load('Mesh/x.npy')
y = np.load('Mesh/y.npy')

x = (x[:-1, :-1] + x[1:, 1:]).transpose() / 2
y = (y[:-1, :-1] + y[1:, 1:]).transpose() / 2

rho = rho[2:-2, 2:-2]
u = u[2:-2, 2:-2]
v = v[2:-2, 2:-2]
p = p[2:-2, 2:-2]

gamma = 1.4
e = p / (gamma - 1) + 0.5 * rho * (u ** 2 + v ** 2)
ma = np.sqrt(u ** 2 + v ** 2) / np.sqrt(gamma * p / rho)


plt.contourf(x, y, p, levels=100)
plt.imshow(p)

plt.colorbar()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Pressure')

# Wave angle should be 45.3436167
x_end = 1 + 2.4 / np.tan(np.deg2rad(45.3436167))
plt.plot([1, x_end], [0, 2.4], 'cyan', linewidth=1, linestyle='--')
plt.legend(['Analytical Shockwave position'])

nx = p.shape[0]
ny = p.shape[1]

print("Pressure Ratio before and after shockwave: ", p[:,int(nx*11/16)].max() / p[:,int(nx*1/8)].max())

plt.savefig('Pressure.png')
plt.show()

plt.close()

plt.contourf(x, y, ma, levels=100)
# plt.imshow(p)

plt.colorbar()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Mach Number')

# Wave angle should be 45.3436167
x_end = 1 + 2.4 / np.tan(np.deg2rad(45.3436167))
plt.plot([1, x_end], [0, 2.4], 'cyan', linewidth=1, linestyle='--')
plt.legend(['Analytical Shockwave position'])

print("Mach Number after shockwave: ", ma[8,:].min())

plt.savefig('Mach.png')

plt.close()

