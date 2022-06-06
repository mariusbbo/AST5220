import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({'font.size': 16})

# Open and read file with values for matter power spectrum and theta_l
file_mpk = open("mpk_keq.txt", "r")
keq = float(file_mpk.readline()) # Value of k corresponding to the particle horizon at matter-radiation equality
leta0 = np.array(file_mpk.readline().split(), dtype=float) # k = l*eta0

lines = np.array(file_mpk.read().split(), dtype=float)
k = lines[::8]
mpk = lines[1::8]
ck_over_H0 = lines[2::8]

# Data from observations
SDSS_data = np.loadtxt("reid_DR7.txt")
WMAP_data = np.loadtxt("wmap_act.txt")

idx_peak = np.argmin(max(mpk)-mpk)

# Plot matter power spectrum
plt.plot(k, mpk)
plt.plot(k[:idx_peak], k[:idx_peak]**0.965 * 3e6, ls="--", label="$k^{n_s}$")
plt.vlines(keq, ymin=min(mpk), ymax=max(mpk)+0.2*max(mpk), color="k", ls="--", label=(r"$k_{\rm eq}$ = " + "{:.2f}".format(keq)))
# plt.errorbar(SDSS_data[:,0], SDSS_data[:,1], SDSS_data[:,2], ecolor="k", fmt="o", color="k",
#              label="SDSS Galaxies (DR7 LRG)", markersize=3)
# plt.errorbar(WMAP_data[:,0], WMAP_data[:,1], np.array([np.zeros(len(WMAP_data)), WMAP_data[:,2]]), ecolor="y", fmt="o", color="y",
#              label="CMB (WMAP+ACT)", markersize=3)
plt.title("Total matter power spectrum")
plt.xlabel(r"Wavenumber $k$ [$h$/Mpc]")
plt.ylabel(r"$P(k)$ [(Mpc/$h)^3$]")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper left")
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/mpk.pdf")
plt.show()

# List with values of ell for which theta_l will be plotted
l_list = [6, 100, 200, 500, 1000]

# Plot theta
fig, ax = plt.subplots()

for i in range(3, 8):
    idx_peak = np.argmin(max(lines[i::8])-lines[i::8])
    ax.plot(k, lines[i::8]*100, label="$\ell$ = {}".format(l_list[i-3]))
    plt.axvline(leta0[i-3], ymin=-1, ymax=2, color="k", ls="--", linewidth=1.5)

ax.set_title("Temperature multipoles")
ax.set_xlabel(r"$k$ [$h$/Mpc]")
ax.set_ylabel(r"$\theta_\ell(k) \times 10^{-2}$")
ax.set_xlim(-5e-3,0.2)
ax.set_ylim(-1,2)
ax.legend()
ax.yaxis.set_label_coords(-0.1, 0.55)
ax.xaxis.set_major_locator(MaxNLocator(5))
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/theta.pdf")
plt.show()


# Plot theta^2/k
fig, ax = plt.subplots()
for i in range(3, 7):
    ax.plot(k, lines[i::8]**2 / k * 1000, label=r"$\ell$={}".format(l_list[i-3]))
    plt.axvline(leta0[i-3], ymin=-1, ymax=2, color="k", ls="--", linewidth=1.5)

ax.set_title("Spectrum integrand")
ax.set_xlabel(r"$k$ [$h$/Mpc]")
ax.set_ylabel(r"$\frac{\theta_\ell^2(k)}{k} \times 10^{-3}$")
ax.set_xlim(-5e-3,0.1)
ax.set_ylim(-0.1, max(lines[4::8]**2 / k * 1000)+0.1*max(lines[4::8]**2 / k * 1000))
ax.legend()
ax.yaxis.set_label_coords(-0.08, 0.55)
ax.xaxis.set_major_locator(MaxNLocator(5))
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/theta2_l100.pdf")
plt.show()


# Load file with values of the CMB power spectrum
file_CMB = np.loadtxt("cells.txt")
ell = file_CMB[:,0]
Cell = file_CMB[:,1]
Cell_noscaling = file_CMB[:,2]

# # Plot CMB power spectrum without scaling with l(l+1)/2pi
# plt.plot(ell, Cell_noscaling)
# plt.xlabel("l")
# plt.ylabel("Cl")
# plt.xscale("log")
# plt.yscale("log")
# plt.show()

# Load CMB data from observations
file_planck = np.loadtxt("planck_cell_low.txt")
ellp = file_planck[:,0]
Cellp = file_planck[:,1]
err_up = file_planck[:,2]
err_down = file_planck[:,3]

file_planck_high_l = np.loadtxt("COM_PowerSpect_CMB-TT-binned_R3.01.txt")
ellp_h = file_planck_high_l[:,0]
err_up_h = file_planck_high_l[:,2]
Cellp_h = file_planck_high_l[:,1]
err_down_h = file_planck_high_l[:,3]

# Plot the CMB power spectrum
plt.plot(ell, Cell)
plt.errorbar(ellp, Cellp, np.array([err_down, err_up]), ecolor="k", fmt="o", color="k", label="Planck data", markersize=3)
plt.errorbar(ellp_h, Cellp_h, np.array([err_down_h, err_up_h]), ecolor="k", fmt="o", color="k", markersize=3)
plt.title("CMB angular power spectrum")
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell(\ell+1)C_\ell$ [$\mu$K$^2$]")
plt.xscale("log")
plt.legend(loc="upper left")
figure = plt.gcf()
figure.set_size_inches(8.5, 6)
# plt.savefig("Images/CMB.pdf")
plt.show()


# # Plot function
# def plot_CMB(x, y, label):
#     plt.plot(x, y, label="{}".format(label))
#     plt.title("CMB angular power spectrum")
#     plt.xlabel(r"Multipole $\ell$")
#     plt.ylabel(r"$\ell(\ell+1)C_\ell$ [$\mu$K$^2$]")
#     plt.xscale("log")
#     plt.legend(loc="upper left")
#     figure = plt.gcf()
#     figure.set_size_inches(8.5, 6)
#
# # Plot full CMB power spectrum
# plot_CMB(ell, Cell, "Total")
#
# # Plot CMB power spectrum with only the SW term in the source function
# file_CMB_SW = np.loadtxt("cells_sf_SW.txt")
# ell_SW = file_CMB_SW[:,0]
# Cell_SW = file_CMB_SW[:,1]
# plot_CMB(ell_SW, Cell_SW, "SW")
#
# # Plot CMB power spectrum with only the SW term in the source function
# file_CMB_ISW = np.loadtxt("cells_sf_ISW.txt")
# ell_ISW = file_CMB_ISW[:,0]
# Cell_ISW = file_CMB_ISW[:,1]
# plot_CMB(ell_ISW, Cell_ISW, "ISW")
#
# # Plot CMB power spectrum with only the SW term in the source function
# file_CMB_Doppler = np.loadtxt("cells_sf_Doppler.txt")
# ell_Doppler = file_CMB_Doppler[:,0]
# Cell_Doppler = file_CMB_Doppler[:,1]
# plot_CMB(ell_Doppler, Cell_Doppler, "Doppler")
#
# # Plot CMB power spectrum with only the SW term in the source function
# file_CMB_dd = np.loadtxt("cells_sf_dd.txt")
# ell_dd = file_CMB_dd[:,0]
# Cell_dd = file_CMB_dd[:,1]
# plot_CMB(ell_dd, Cell_dd, "dd")
# plt.show()

# # Plot spherical Bessel function
# file_bessel = np.loadtxt("bessel.txt")
# z = file_bessel[:,0]
#
# l_list2 = [2, 15, 60, 900, 1900]
#
# plt.plot(z, file_bessel[:,1], label="l=2")
# plt.plot(z, file_bessel[:,2])
# plt.plot(z, file_bessel[:,3], label="l=15")
# plt.plot(z, file_bessel[:,4])
# plt.plot(z, file_bessel[:,5], label="l=60")
# plt.plot(z, file_bessel[:,6])
# plt.plot(z, file_bessel[:,7], label="l=900")
# plt.plot(z, file_bessel[:,8])
# plt.plot(z, file_bessel[:,9], label="l=1900")
# plt.plot(z, file_bessel[:,10])
# plt.xlabel("z")
# plt.ylabel("j_ell")
# plt.xscale("log")
# plt.legend()
# plt.show()
#
#
# # Plot source function
# file_sf = np.loadtxt("perturbations_k0.1_sf.txt")
# x = file_sf[:,0]
# sf = file_sf[:,1]
#
# plt.plot(x, sf)
# plt.xlabel("x")
# plt.xlabel("S")
# plt.show()
#
# def plot(ax, x, y, label, title, xlabel, ylabel):
#     ax.plot(x, y, label=label)
#     ax.set_title(title)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#
# file_sf1 = np.loadtxt("perturbations_k0.1_sf_test.txt")
# file_sf2 = np.loadtxt("perturbations_k0.01_sf_test.txt")
# file_sf3 = np.loadtxt("perturbations_k0.001_sf_test.txt")
# file_sf4 = np.loadtxt("perturbations_k0.0001_sf_test.txt")
# x = file_sf1[:,0]
# sf1 = file_sf1[:,1]
# sf2 = file_sf2[:,1]
# sf3 = file_sf3[:,1]
# sf4 = file_sf4[:,1]
#
# plt.plot(x, sf1, ls="--", label="k=0.1")
# plt.plot(x, sf2, ls="--", label="k=0.01")
# plt.plot(x, sf3, ls="--", label="k=0.001")
# plt.plot(x, sf4, ls="--", label="k=0.0001")
# plt.xlabel("x")
# plt.ylabel("S")
# plt.legend()
# plt.show()
#
#
# fig, ax = plt.subplots(2,2,figsize=(8,8))
#
# l_list = [6, 100, 200, 500, 1000]
#
# for j in range(2, 7):
#     # plot(ax[0][0], x, file_sf1[:,j], "l = {}".format(l_list[j-2]), "k = 0.1", "x", "j_l")
#     plot(ax[0][0], x, file_sf1[:,j], None, "k = 0.1", "x", "j_l")
#     plot(ax[0][1], x, file_sf2[:,j], None, "k = 0.01", "x", "j_l")
#     plot(ax[1][0], x, file_sf3[:,j], None, "k = 0.001", "x", "j_l")
#     plot(ax[1][1], x, file_sf4[:,j], None, "k = 0.0001", "x", "j_l")
#
# ax[0][0].legend(loc="upper left")
# plt.tight_layout()
# plt.show()
#
# fig, ax = plt.subplots(2,2,figsize=(8,8))
#
# for j in range(7, 12):
#     # plot(ax[0][0], x, file_sf1[:,j], "l = {}".format(l_list[j-7]), "k = 0.1", "x", "S*j_l")
#     plot(ax[0][0], x, file_sf1[:,j], None, "k = 0.1", "x", "S*j_l")
#     plot(ax[0][1], x, file_sf2[:,j], None, "k = 0.01", "x", "S*j_l")
#     plot(ax[1][0], x, file_sf3[:,j], None, "k = 0.001", "x", "S*j_l")
#     plot(ax[1][1], x, file_sf4[:,j], None, "k = 0.0001", "x", "S*j_l")
#
# ax[0][0].legend(loc="upper left")
# plt.tight_layout()
# plt.show()
