import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt
import sys

file_format = "csv"
file_name = "tracers_0." + file_format

print("Reading file:", file_name)

particles = np.loadtxt(file_name, delimiter=",", skiprows=1)


wall_bounce = True
n_particles = 3
t = []
tstep = []
x = [None] * n_particles
v = [None] * n_particles
d = [None] * n_particles
rho = [None] * n_particles
    
def exact_solution(x0, v0, t, d_p, rho_p, wall_bounce=False):
    rho_f = 1.0
    mu_f = 1.0

    tau_p = rho_p * d_p**2 / (18.0 * mu_f)

    Re0 = rho_f * d_p * abs(v0) / mu_f

    m = 0.687
    alpha = 0.15

    t = np.asarray(t)

    vel = (
        v0
        * np.exp(-t / tau_p)
        * (
            1.0
            + alpha * Re0**m
            * (1.0 - np.exp(-m * t / tau_p))
        )**(-1.0 / m)
    )

    disp = cumulative_trapezoid(vel, t, initial=0.0)
    if wall_bounce:
        disp[disp > 0.5 - d_p/2] = (0.5 - d_p/2)*2 - disp[disp > 0.5 - d_p/2]
        return x0 + disp
    else:
        return x0 + disp

print("Processing", n_particles, "particles")
for i in range(particles.shape[0]):
    tstep_current = particles[i,0]
    if tstep_current not in tstep:
        t.append(particles[i,1])
        tstep.append(tstep_current)
    pt_id = int(particles[i,2]) - 1
    if x[pt_id] is None:
        x[pt_id] = []
    x[pt_id].append(particles[i,3])
    if v[pt_id] is None:
        v[pt_id] = []
    v[pt_id].append(particles[i,6])
    if d[pt_id] is None:
        d[pt_id] = []
    d[pt_id].append(particles[i,9])
    if rho[pt_id] is None:
        rho[pt_id] = []
    rho[pt_id].append(particles[i,10])
for i in range(n_particles):
    if x[i] is None:
        x[i] = []
    x[i] = np.array(x[i])
    if v[i] is None:
        v[i] = []
    v[i] = np.array(v[i])
    if d[i] is None:
        d[i] = []
    d[i] = np.array(d[i])
    if rho[i] is None:
        rho[i] = []
    rho[i] = np.array(rho[i])
    
t = np.array(t)

print("Plotting particle trajectories")
plt.figure(figsize=(8, 6))
for i in range(n_particles):
    t_plot = t[:len(x[i])]
    x_plot = x[i]
    d_p = d[i][0]
    rho_p = rho[i][0]
    x_exact = exact_solution(x_plot[0], v[i][0], \
                             t_plot, d_p, rho_p, wall_bounce=wall_bounce)
    plt.plot(t_plot, x_plot,color="blue", linewidth=2)
    plt.plot(t_plot, x_exact, color="red", linewidth=2, linestyle="--")

plt.plot([],[], color="blue", linestyle="-", label="Numerical Solution")
plt.plot([],[], color="red", linestyle="--", label="Exact Solution")
plt.xlabel("Time")
plt.ylabel("Position")
plt.title("Particle Trajectories")
plt.legend()
plt.grid()
plt.savefig("particle_trajectories.png")

plt.figure(figsize=(8, 6))
for i in range(n_particles):
    t_plot = t[:len(x[i])]
    v_plot = v[i]
    plt.plot(t_plot, v_plot, linewidth=2, label=f"Particle {i+1}")
plt.xlabel("Time")
plt.ylabel("u")
plt.xlim(0, 4.0)
plt.title("Velocity history")
plt.legend()
plt.grid()
plt.savefig("velocity.png")