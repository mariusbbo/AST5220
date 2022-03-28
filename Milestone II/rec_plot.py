import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

# Read x and z values for recombination and decoupling
rec_times_file = open("rec_times.txt", "r")
x_vals = np.array(rec_times_file.readline().split(), dtype=float)
z_vals = np.array(rec_times_file.readline().split(), dtype=float)
Xe_freeze_out = float(rec_times_file.readline()) # Freeze-out abundance Xe

# Create arrays with results from simulation
rec_file = np.loadtxt("recombination.txt")
x = rec_file[:,0]
Xe = rec_file[:,1]
Xe_saha = rec_file[:,2]
ne = rec_file[:,3]
tau = rec_file[:,4]
dtaudx = rec_file[:,5]
ddtauddx = rec_file[:,6]
g = rec_file[:,7]
dgdx = rec_file[:,8]
ddgddx = rec_file[:,9]
T = rec_file[:,10]

a = np.exp(x)
z = 1/a - 1

# Find index of x-array when x ~ -7.3
idx_x_drop = np.argmin(abs(x+7.3))

k_B = 1.38064852e-23 # Boltzmann constant
J_per_ev = 1.60217653e-19 # Joule per electronvolt

# print("Temperature when Xe starts to drop: {} K, {:.2f} eV".format(T[idx_x_drop], T[idx_x_drop]*k_B/J_per_ev))

# Plot the fractional electron density Xe as function of x
plt.plot(x, Xe)
plt.hlines(Xe_freeze_out, xmin=x[0], xmax=x[-1], color="k", ls="--", label=r"$X_e = 2.026\cdot$10$^{-4}$")
plt.title("Free electron fraction")
plt.xlabel("x")
plt.ylabel(r"$X_e$")
plt.yscale("log")
plt.xlim(-12,0)
plt.ylim(1e-4,1.5)
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.grid()
plt.legend(loc="upper right", framealpha=1)
# plt.savefig("Images/Xe_of_x.pdf")
plt.show()


# Plot fractional electron density Xe computed only with Saha equation
# and Xe computed with Saha and Peebles
plt.plot(x, Xe_saha, label="Saha")
plt.plot(x, Xe, label="Saha + Peebles")
plt.hlines(0.5, xmin=x[0], xmax=x[-1], color="k", ls="--", label=r"$X_e$ = 0.5")
plt.vlines(x_vals[0], ymin=1e-4, ymax=0.3, color="k", ls="dashdot", label=r"$\tau$ = 1")
plt.title("Free electron fraction")
plt.xlabel("x")
plt.ylabel(r"$X_e$")
plt.yscale("log")
plt.xlim(-8,-5.5)
plt.ylim(1e-4,1.5)
plt.grid()
plt.legend(loc="upper right", framealpha=1)
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/Xe_comp.pdf")
plt.show()

# Index of x at matter-radiation equality
idx_mr = np.argmin(abs(x+8.13))

dtaudx_rad = np.exp(-x) / 5
dtaudx_mat = np.exp(-3*x/2) / 1e6

# Plot optical depth tau and derivatives dtau/dx and d^2tau/dx^2
# Also plot the approximations to dtau/dx
plt.plot(x, tau, label=r"$\tau(x)$")
plt.plot(x, -dtaudx, label=r"$-\tau'(x)$")
plt.plot(x, ddtauddx, label=r"$\tau''(x)$")
plt.plot(x[:idx_mr], dtaudx_rad[:idx_mr], color="k", ls="--", label=r"$e^{-x}/5$")
plt.plot(x[idx_mr:], dtaudx_mat[idx_mr:], color="k", ls="dashdot", label=r"$e^{-3x/2}\cdot10^{-6}$")
plt.title("Optical depth")
plt.xlabel("x")
plt.ylabel(r"$\tau$, $\tau'$, $\tau''$")
plt.yscale("log")
plt.xlim(-12, 0)
plt.ylim(1e-8, 1e6)
plt.legend(framealpha=1)
figure = plt.gcf()
figure.set_size_inches(8, 6)
plt.grid()
# plt.savefig("Images/tau.pdf")
plt.show()


# Find x and z for the peak of the visibility function
idx_max_g = np.argmax(g)
x_g = x[idx_max_g]
z_g = np.exp(-x_g) - 1
# print("x at peak of visibility function g: {}".format(x_g))
# print("z at peak of visibility function g: {}".format(z_g))

# Plot visibility function g and derivatives dg/dx and d^2g/dx^2
fig, ax = plt.subplots()
ax.plot(x, g, label=r"$g(x)$")
ax.plot(x, dgdx/10, ls="--", label=r"$g'(x)$")
ax.plot(x, ddgddx/150, ls="--", label=r"$g''(x)$")
ax.set_title("Visibility function")
ax.set_xlabel("x")
ax.set_ylabel(r"$\tilde{g}$, $\tilde{g}'$, $\tilde{g}''$")
ax.set_xlim(-7.5,-6)
ax.legend(framealpha=1)
ax.grid()
ax.yaxis.set_label_coords(-0.1, 0.55)
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/visibility_function.pdf")
plt.show()
