import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

# Load files with perturbation variables for different values of k
file_k0_1 = np.loadtxt("perturbations_k0.1.txt")
file_k0_01 = np.loadtxt("perturbations_k0.01.txt")
file_k0_001 = np.loadtxt("perturbations_k0.001.txt")

pert_array = np.array([file_k0_1, file_k0_01, file_k0_001])

x = pert_array[0][:,0]
delta_cdm = pert_array[:,:,1]
delta_b = pert_array[:,:,2]
v_cdm = pert_array[:,:,3]
v_b = pert_array[:,:,4]
theta0 = pert_array[:,:,5]
theta1 = pert_array[:,:,6]
theta2 = pert_array[:,:,7]
phi = pert_array[:,:,8]
psi = pert_array[:,:,9]
pi = pert_array[:,:,10]
source_f = pert_array[:,:,11]


# Gravitational potential perturbations Phi and -Psi
plt.plot(x, phi[0], "b", label=r"k=0.1/Mpc")
plt.plot(x, phi[1], "r", label=r"k=0.01/Mpc")
plt.plot(x, phi[2], "g", label=r"k=0.001/Mpc")
plt.plot(x, -psi[0], "b--")
plt.plot(x, -psi[1], "r--")
plt.plot(x, -psi[2], "g--")
plt.vlines(-8.13, ymin=phi[0].min(), ymax=phi[0].max(), color="k", label=r"Mat-rad eq")
plt.vlines(-6.988, ymin=phi[0].min(), ymax=phi[0].max(), color="k", ls="--", label=r"Decoupling")
plt.title(r"Newtonian potentials")
plt.xlabel("x")
plt.ylabel(r"$\Phi$ (solid),  $-\Psi$ (dashed)")
plt.xlim(-19,1)
plt.legend()
plt.grid()
figure = plt.gcf()
figure.set_size_inches(8, 6)
# plt.savefig("Images/phi_and_psi_new.pdf")
plt.show()



# Density perturbations, delta_cdm, delta_b
plt.plot(x, delta_cdm[0], "b", label=r"k=0.1/Mpc")
plt.plot(x, abs(delta_b[0]), "b--")
plt.plot(x, abs(4*theta0[0]), color="b", ls="dotted")
plt.plot(x, delta_cdm[1], "r", label=r"k=0.01/Mpc")
plt.plot(x, abs(delta_b[1]), "r--")
plt.plot(x, abs(4*theta0[1]), color="r", ls="dotted")
plt.plot(x, delta_cdm[2], "g", label=r"k=0.001/Mpc")
plt.plot(x, abs(delta_b[2]), "g--")
plt.plot(x, abs(4*theta0[2]), color="g", ls="dotted")
plt.vlines(-8.13, ymin=2e-3, ymax=delta_cdm[0].max(), color="k", label=r"Mat-rad eq")
plt.vlines(-6.988, ymin=2e-3, ymax=delta_cdm[0].max(), color="k", ls="--", label=r"Decoupling")
plt.title(r"Density perturbations")
plt.xlabel("x")
plt.ylabel(r"$\delta_{\rm CDM}$ (solid),  |$\delta_{\rm b}$| (dashed),  |$\delta_\gamma$| (dotted)")
plt.xlim(-19,1)
plt.yscale("log")
plt.legend()
plt.grid()
figure = plt.gcf()
figure.set_size_inches(8.5, 6.5)
# plt.savefig("Images/delta_m_and_gamma_new.pdf")
plt.show()


file_bg = np.loadtxt("cosmology.txt")
x_bg = file_bg[:,0]
Hp = file_bg[:,3]

# Velocities, v_cdm, v_b
plt.plot(x_bg[100:-90], 1/Hp[100:-90]/5e17, color="orange", label=r"$\mathcal{H}^{-1}$")
plt.plot(x, v_cdm[0], "b", label=r"k=0.1/Mpc")
plt.plot(x, abs(v_b[0]), "b--")
# plt.plot(x, abs(-3*theta1[0]), color="b", ls="dotted")
plt.plot(x, v_cdm[1], "r", label=r"k=0.01/Mpc")
plt.plot(x, abs(v_b[1]), "r--")
# plt.plot(x, abs(-3*theta1[1]), color="r", ls="dotted")
plt.plot(x, v_cdm[2], "g", label=r"k=0.001/Mpc")
plt.plot(x, abs(v_b[2]), "g--")
# plt.plot(x, abs(-3*theta1[2]), color="g", ls="dotted")
plt.vlines(-8.13, ymin=v_cdm[2].min(), ymax=v_cdm[0].max(), color="k", label=r"Mat-rad eq")
plt.vlines(-6.988, ymin=v_cdm[2].min(), ymax=v_cdm[0].max(), color="k", ls="--", label=r"Decoupling")
plt.title(r"Velocity perturbations, matter")
plt.xlabel("x")
plt.ylabel(r"$v_{\rm CDM}$ (solid),  |$v_{\rm b}$| (dashed)")
plt.xlim(-19,1)
plt.ylim(1e-6,2000)
plt.yscale("log")
plt.legend(loc="upper left")
plt.grid()
figure = plt.gcf()
figure.set_size_inches(8.5, 6.5)
# plt.savefig("Images/v_m_new.pdf")
plt.show()

# Velocities, v_gamma = -3*theta_1
plt.plot(x_bg[100:-90], 1/Hp[100:-90]/5e17, color="orange", label=r"$\mathcal{H}^{-1}$")
plt.plot(x, abs(-3*theta1[0]), color="b", label=r"k=0.1/Mpc")
plt.plot(x, abs(-3*theta1[1]), color="r", label=r"k=0.01/Mpc")
plt.plot(x, abs(-3*theta1[2]), color="g", label=r"k=0.001/Mpc")
plt.vlines(-8.13, ymin=v_cdm[2].min(), ymax=v_cdm[0].max(), color="k", label=r"Mat-rad eq")
plt.vlines(-6.988, ymin=v_cdm[2].min(), ymax=v_cdm[0].max(), color="k", ls="--", label=r"Decoupling")
plt.title(r"Velocity perturbations, photons")
plt.xlabel("x")
plt.ylabel(r"|$v_\gamma$|")
plt.xlim(-19,1)
plt.ylim(1e-6,2000)
plt.yscale("log")
plt.legend(loc="upper left")
plt.grid()
figure = plt.gcf()
figure.set_size_inches(8.5, 6.5)
# plt.savefig("Images/v_gamma_new.pdf")
plt.show()


# # Gravitational potential perturbation Phi
# plt.plot(x, phi[0], label=r"k=0.1/Mpc")
# plt.plot(x, phi[1], label=r"k=0.01/Mpc")
# plt.plot(x, phi[2], label=r"k=0.001/Mpc")
# plt.vlines(-8.13, ymin=phi[0].min(), ymax=phi[0].max(), color="k", linewidth=2, label=r"Mat-rad eq")
# plt.title(r"Newtonian potential")
# plt.xlabel("x")
# plt.ylabel(r"$\Phi(x,k)$")
# plt.xlim(-19,1)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# # plt.savefig("Images/phi.pdf")
# plt.show()


# # Gravitational potential perturbation Psi
# plt.plot(x, psi[0], label=r"k=0.1/Mpc")
# plt.plot(x, psi[1], label=r"k=0.01/Mpc")
# plt.plot(x, psi[2], label=r"k=0.001/Mpc")
# plt.vlines(-8.13, ymin=psi[0].min(), ymax=psi[0].max(), color="k", linewidth=2, label=r"Mat-rad eq")
# plt.title(r"Newtonian potential")
# plt.xlabel("x")
# plt.ylabel(r"$\Psi(x,k)$")
# plt.xlim(-19,1)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8.5, 6)
# # plt.savefig("Images/psi.pdf")
# plt.show()


# # Photons density delta_gamma = 4*theta0
# plt.plot(x, 4*theta0[0], label=r"k=0.1/Mpc")
# plt.plot(x, 4*theta0[1], label=r"k=0.01/Mpc")
# plt.plot(x, 4*theta0[2], label=r"k=0.001/Mpc")
# plt.vlines(-8.13, ymin=min(4*theta0[0])-0.5, ymax=max(4*theta0[0])+0.5, color="k", linewidth=2, label=r"Mat-rad eq")
# plt.title("Photon density perturbation")
# plt.xlabel("x")
# plt.ylabel(r"$\delta_\gamma = 4\theta_0$")
# plt.xlim(-19,1)
# plt.legend(loc="lower left")
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# # plt.savefig("Images/delta_gamma.pdf")
# plt.show()


# # Photons velocity v_gamma = -3*theta1
# plt.plot(x, -3*theta1[0], label=r"k=0.1/Mpc")
# plt.plot(x, -3*theta1[1], label=r"k=0.01/Mpc")
# plt.plot(x, -3*theta1[2], label=r"k=0.001/Mpc")
# plt.vlines(-8.13, ymin=min(-3*theta1[0])-0.1, ymax=max(-3*theta1[0])+0.2, color="k", linewidth=2, label=r"Mat-rad eq")
# plt.title("Photon velocity perturbation")
# plt.xlabel("x")
# plt.ylabel(r"$v_\gamma = -3\theta_1$")
# plt.xlim(-19,1)
# plt.legend(loc="upper left")
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# # plt.savefig("Images/v_gamma.pdf")
# plt.show()


# # Phi+Psi
# plt.plot(x, phi[0]+psi[0], label=r"k=0.1/Mpc")
# plt.plot(x, phi[1]+psi[1], label=r"k=0.01/Mpc")
# plt.plot(x, phi[2]+psi[2], label=r"k=0.001/Mpc")
# # plt.title(r"$\Phi+\Psi$")
# plt.xlabel("x")
# plt.ylabel(r"$\Phi+\Psi$")
# plt.xlim(-19,1)
# # plt.xscale("log")
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# plt.savefig("Images/phi_pluss_psi.pdf")
# plt.show()


# # Temperature perturbation theta0
# plt.plot(x, theta0[0], label=r"k=0.1/Mpc")
# plt.plot(x, theta0[1], label=r"k=0.01/Mpc")
# plt.plot(x, theta0[2], label=r"k=0.001/Mpc")
# plt.xlabel("x")
# plt.ylabel(r"$\theta_0$")
# plt.xlim(-19,1)
# # plt.ylim(-0.5,0.9)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# plt.savefig("Images/theta0.pdf")
# plt.show()


# # Effective temperature perturbation theta0 + Psi
# plt.plot(x, theta0[0] + psi[0], label=r"k=0.1/Mpc")
# plt.plot(x, theta0[1] + psi[1], label=r"k=0.01/Mpc")
# plt.plot(x, theta0[2] + psi[2], label=r"k=0.001/Mpc")
# plt.title(r"$\theta_0+\Psi$")
# plt.xlabel("x")
# plt.ylabel(r"$\theta_0+\Psi$")
# plt.xlim(-19,1)
# # plt.ylim(-0.5,0.9)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# plt.savefig("Images/theta0_pluss_psi.pdf")
# plt.show()


# # Temperature perturbation theta1
# plt.plot(x, theta1[0], label=r"k=0.1/Mpc")
# plt.plot(x, theta1[1], label=r"k=0.01/Mpc")
# plt.plot(x, theta1[2], label=r"k=0.001/Mpc")
# plt.xlabel("x")
# plt.ylabel(r"$\theta_1$")
# plt.xlim(-19,1)
# # plt.ylim(-0.3,0.4)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# plt.savefig("Images/theta1.pdf")
# plt.show()


# # Temperature perturbation theta2
# plt.plot(x, theta2[0], label=r"k=0.1/Mpc")
# plt.plot(x, theta2[1], label=r"k=0.01/Mpc")
# plt.plot(x, theta2[2], label=r"k=0.001/Mpc")
# plt.xlabel("x")
# plt.ylabel(r"$\theta_2$")
# plt.xlim(-19,1)
# # plt.ylim(-0.3,0.4)
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# # plt.savefig("Images/theta2.pdf")
# plt.show()


# # Source function
# plt.plot(x, source_f[0], label=r"k=0.1/Mpc")
# plt.plot(x, source_f[1], label=r"k=0.01/Mpc")
# plt.plot(x, source_f[2], label=r"k=0.001/Mpc")
# plt.title(r"Source function")
# plt.xlabel("x")
# plt.ylabel(r"$S$")
# plt.xlim(-19,1)
# # plt.xscale("log")
# plt.legend()
# plt.grid()
# figure = plt.gcf()
# figure.set_size_inches(8, 6)
# plt.savefig("Images/source_function.pdf")
# plt.show()
