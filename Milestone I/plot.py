import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

# Read parameter values from file
file = open("cosmology.txt", "r")
lines = np.array(file.read().split(), dtype=float)
x = lines[::22]
eta_in_Mpc = lines[1::22]
H_over_H0 = lines[2::22]
Hp = lines[3::22]
dHpdx_over_Hp = lines[4::22]
ddHpddx_over_Hp = lines[5::22]
OmegaB = lines[6::22]
OmegaCDM = lines[7::22]
OmegaLambda = lines[8::22]
OmegaGamma = lines[9::22]
OmegaNu = lines[10::22]
OmegaK = lines[11::22]
dHpdx_analytical = lines[12::22]
ddHpddx_analytical = lines[13::22]
t = lines[14::22]
comoving_dist = lines[15::22]
angular_dist = lines[16::22]
luminosity_dist = lines[17::22]
H = lines[18::22]
eta = lines[19::22]
x_of_t = lines[20::22]
a_of_t = lines[21::22]

a = np.exp(x) # scale factor
z = 1/a - 1 # redshift
c = 299792458 # speed of light in m/s
H0 = 67 * 10**3 # Hubble constant in unit m/s/Mpc
s_in_Gyr = 10**9 * 365 * 24 * 60 * 60 # seconds in a gigayear

index_t0 = np.argmin(abs(x)) # find index of x today
H0_s = H[index_t0] # Hubble constant in unit 1/s

print("Age of Universe today: {:.2f} Gyrs".format(t[index_t0]/s_in_Gyr))
print("Conformal time today eta0/c: {} Gyrs".format(eta[index_t0]/c/s_in_Gyr))

# Read x, a, z and t for matter-radiation equality, matter-dark energy equality
# and when the Universe started accelerating
time_file = open("times.txt", "r")
time_lines = np.array(time_file.read().split(), dtype=float)
x_list = time_lines[::4]
z_list = time_lines[1::4]
a_list = time_lines[2::4]
t_list = time_lines[3::4]

# Density parameters
OmegaMatter = OmegaB + OmegaCDM
OmegaRad = OmegaGamma + OmegaNu
OmegaTot = OmegaB + OmegaCDM + OmegaGamma + OmegaNu + OmegaLambda + OmegaK

# Plot evolution of Omegas
plt.plot(x, OmegaMatter/OmegaTot, label=r"$\Omega_{\rm m}$")
plt.plot(x, OmegaRad/OmegaTot, label=r"$\Omega_{\rm r}$")
plt.plot(x, OmegaLambda/OmegaTot, label=r"$\Omega_{\Lambda}$")
plt.plot(x, OmegaK/OmegaTot, label=r"$\Omega_{k}$")
plt.vlines(x_list[0], ymin=0, ymax=1, color="k", ls="--",
           label="Matter-radiation equality")
plt.vlines(x_list[2], ymin=0, ymax=1, color="k", ls="dotted",
           label="Universe starts accelerating")
plt.vlines(x_list[1], ymin=0, ymax=1, color="k", ls="dashdot",
           label="Matter-dark energy equality")
plt.title("Evolution of density parameters")
plt.xlabel("x")
plt.ylabel(r"$\Omega$")
plt.grid(alpha=0.3)
plt.legend()
figure = plt.gcf()
figure.set_size_inches(20, 7)
# plt.savefig("Images/omegas_of_x.pdf")
plt.show()

km_per_Mpc = 3.08567758e17 # km in a Mpc
H_new = H * km_per_Mpc # Hubble parameter in units of km/s/Mpc

# Plot Hubble parameter H(x)
plt.plot(x, H_new)
plt.vlines(x_list[0], ymin=min(H_new), ymax=1e11, color="k", ls="--",
           label="Matter-radiation equality")
plt.vlines(x_list[2], ymin=min(H_new), ymax=1e11, color="r", ls="--",
           label="Universe starts accelerating")
plt.vlines(x_list[1], ymin=min(H_new), ymax=1e11, color="g", ls="--",
           label="Matter-dark energy equality")
# plt.title("Evolution of Hubble parameter")
plt.xlabel("x")
plt.ylabel(r"$H(x)$ [100 km/s/Mpc]")
plt.yscale("log")
plt.grid(alpha=0.3)
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/H_new_unit.pdf")
plt.show()


Hp_new = Hp * km_per_Mpc # Hp in units km/s/Mpc

# Plot scaled Hubble parameter Hp(x)
plt.plot(x, Hp_new)
plt.vlines(x_list[0], ymin=0.4, ymax=8e4, color="k", ls="--",
           label="Matter-radiation equality")
plt.vlines(x_list[2], ymin=0.4, ymax=8e4, color="r", ls="--",
           label="Universe starts accelerating")
plt.vlines(x_list[1], ymin=0.4, ymax=8e4, color="g", ls="--",
           label="Matter-dark energy equality")
plt.title("Evolution of $\mathcal{H}$")
plt.xlabel("x")
plt.ylabel(r"$\mathcal{H}(x)$ [100 km/s/Mpc]")
plt.yscale("log")
plt.grid(alpha=0.3)
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/Hp_new_unit.pdf")
plt.show()


# Plot first and second derivative of Hp(x), numerical
# plt.plot(x, dHpdx_over_Hp, label=r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$")
# plt.plot(x, ddHpddx_over_Hp, label=r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
# plt.vlines(x_list[0], ymin=min(ddHpddx_over_Hp), ymax=max(ddHpddx_over_Hp), color="k", ls="--",
#            label="Matter-radiation equality")
# plt.vlines(x_list[2], ymin=min(ddHpddx_over_Hp), ymax=max(ddHpddx_over_Hp), color="r", ls="--",
#            label="Universe starts accelerating")
# plt.vlines(x_list[1], ymin=min(ddHpddx_over_Hp), ymax=max(ddHpddx_over_Hp), color="g", ls="--",
#            label="Matter-dark energy equality")
# plt.title("Derivatives of $\mathcal{H}$, numerical")
# plt.xlabel("x")
# plt.ylabel(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$ , $\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
# plt.grid(alpha=0.3)
# plt.legend()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# # plt.savefig("Images/dHpdx.pdf")
# plt.show()


# Plot first and second derivative of Hp(x), analytical
fig, ax = plt.subplots()

ax.plot(x, dHpdx_analytical/Hp, label=r"$\frac{\mathcal{H}^{\prime}}{\mathcal{H}}$")
ax.plot(x, ddHpddx_analytical/Hp, label=r"$\frac{\mathcal{H}^{\prime\prime}}{\mathcal{H}}$")
ax.vlines(x_list[0], ymin=min(dHpdx_over_Hp), ymax=max(ddHpddx_over_Hp), color="k", ls="--")
ax.vlines(x_list[2], ymin=min(dHpdx_over_Hp), ymax=max(ddHpddx_over_Hp), color="r", ls="--")
ax.vlines(x_list[1], ymin=min(dHpdx_over_Hp), ymax=max(ddHpddx_over_Hp), color="g", ls="--")
ax.set_title("Derivatives of $\mathcal{H}$")
ax.hlines(0, x[0], x[-1], color="k", ls="--")
ax.set_xlabel("x")
ax.set_ylabel(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$ , $\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
ax.set_xlim(-20,3)
ax.set_ylim(-1.8,1.4)
plt.grid(alpha=0.3)
plt.legend(loc="lower right")
ax.yaxis.set_label_coords(-0.08, 0.55)
figure = plt.gcf()
figure.set_size_inches(8, 6.2)
# plt.savefig("Images/derivatives_Hp.pdf")
plt.show()


# Plot conformal time eta(x)
plt.plot(x, eta_in_Mpc)
plt.vlines(x_list[0], ymin=min(eta_in_Mpc), ymax=max(eta_in_Mpc), color="k", ls="--")
plt.vlines(x_list[2], ymin=min(eta_in_Mpc), ymax=max(eta_in_Mpc), color="r", ls="--")
plt.vlines(x_list[1], ymin=min(eta_in_Mpc), ymax=max(eta_in_Mpc), color="g", ls="--")
plt.title("Evolution of conformal time")
plt.xlabel("x")
plt.ylabel(r"$\eta$ [Mpc]")
plt.yscale("log")
plt.grid(alpha=0.3)
plt.tight_layout()
figure = plt.gcf()
figure.set_size_inches(8, 6.2)
# plt.savefig("Images/eta.pdf")
plt.show()


# Plot eta(x) * Hp(x) / c
plt.plot(x, eta*Hp/c)
plt.vlines(x_list[0], ymin=min(eta*Hp/c), ymax=max(eta*Hp/c), color="k", ls="--")
plt.vlines(x_list[2], ymin=min(eta*Hp/c), ymax=max(eta*Hp/c), color="r", ls="--")
plt.vlines(x_list[1], ymin=min(eta*Hp/c), ymax=max(eta*Hp/c), color="g", ls="--")
plt.xlabel("x")
plt.ylabel(r"$\frac{\eta(x)\mathcal{H}(x)}{c}$")
plt.xlim(-15, 1)
plt.ylim(0.5, 6)
plt.grid(alpha=0.3)
plt.tight_layout()
figure = plt.gcf()
figure.set_size_inches(8, 6.2)
# plt.savefig("Images/etaHp_over_c.pdf")
plt.show()


# Plot physical time t(x)
fig, ax = plt.subplots()

ax.plot(x, t/s_in_Gyr)
ax.vlines(x_list[0], ymin=min(t/s_in_Gyr), ymax=max(t/s_in_Gyr), color="k", ls="--")
ax.vlines(x_list[2], ymin=min(t/s_in_Gyr), ymax=max(t/s_in_Gyr), color="r", ls="--")
ax.vlines(x_list[1], ymin=min(t/s_in_Gyr), ymax=max(t/s_in_Gyr), color="g", ls="--")
ax.set_title("Age of Universe")
ax.set_xlabel("x")
ax.set_ylabel(r"$t(x)$ [Gyr]")
ax.set_yscale("log")
ax.grid(alpha=0.3)
ax.yaxis.set_label_coords(-0.1, 0.5)
figure = plt.gcf()
figure.set_size_inches(8, 6.2)
# plt.savefig("Images/time_no_aoUtoday.pdf")
plt.show()


# Read supernova data
file2 = open("Sn_data.txt", "r")
lines2 = np.array(file2.read().split(), dtype=float)
z_Sn = lines2[::3]
dL_Sn = lines2[1::3]
error = lines2[2::3]

# Compare computed luminosity distance with supernova data
plt.plot(z, luminosity_dist/1e3, color="r", label="Theoretical")
plt.errorbar(z_Sn, dL_Sn, yerr=error, ecolor="k", fmt="o", color="k", label="Supernova data", markersize=3)
plt.title("Compare supernova data to luminosity distance")
plt.xlabel("Redshift $z$")
plt.ylabel(r"Distance [Gpc]")
plt.xscale("log")
plt.yscale("log")
plt.xlim(5e-3, 2)
plt.ylim(1e-2, 2e1)
plt.grid(alpha=0.3)
plt.legend()
figure = plt.gcf()
figure.set_size_inches(8, 6.2)
plt.savefig("Images/lum_dist_sn_data.pdf")
plt.show()


# # Plot distance measures
plt.plot(z, c/H0*z, label=r"Hubble distance ($d=\frac{c}{H_0}z$)", color="k", ls="--")
plt.plot(z, luminosity_dist, label="Luminosity distance $d_L$")
plt.plot(z, comoving_dist, label="Comoving distance $r(\chi)$")
plt.plot(z, angular_dist, label="Angular diameter distance $d_A$")
plt.title("Evolution of cosmological distance measures")
plt.xlabel("Redshift $z$")
plt.ylabel(r"Distance [Mpc]")
plt.xscale("log")
plt.yscale("log")
plt.xlim(1e-2, 1e3)
plt.ylim(10, 1e8)
plt.grid(alpha=0.3)
plt.legend()
plt.show()


# Test a(t)
a_rad = (2*H[index_t0]*np.sqrt(OmegaRad[index_t0]))**0.5 * t**0.5
a_mat = (3/2 * H[index_t0] * np.sqrt(OmegaMatter[index_t0]))**(2/3) * t**(2/3)
a_de = np.exp(H[index_t0] * OmegaLambda[index_t0]**0.5 *t)

index_mr_eq = np.argmin(abs(a_of_t-a_list[0]))
index_mde_eq = np.argmin(abs(a_of_t-a_list[1]))

# plt.plot(t, a_of_t, label="numerical")
# plt.plot(t[:index_mr_eq], a_rad[:index_mr_eq], label="rad")
# plt.plot(t[index_mr_eq:index_mde_eq], a_mat[index_mr_eq:index_mde_eq], label="mat")
# plt.plot(t[index_mde_eq:], a_de[index_mde_eq:], label="de")
# plt.hlines(a_list[0], min(t), max(t), color="k", ls="--",
# label="Matter-radiation equality")
# plt.hlines(a_list[1], min(t), max(t), color="k", ls="dashdot",
# label="Matter-dark energy equality")
# plt.hlines(a_list[2], min(t), max(t), color="k", ls="dotted",
# label="Universe starts accelerating")
# plt.xlabel("t")
# plt.ylabel("a")
# plt.xscale("log")
# plt.yscale("log")
# plt.legend()
# plt.show()
