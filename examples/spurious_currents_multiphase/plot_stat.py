# Postprocessing script for spurious currents test.
# Stationary circular drop, La=12000 (rho=300, mu=0.1, sigma=1.0, D=0.4).
# Reads ekin.csv (t, Ekin, enstrophy) and plots:
#   1. Kinetic energy vs time
#   2. Enstrophy vs time
#   3. Spurious capillary number Ca_rms = mu * u_rms / sigma

import numpy as np
import matplotlib.pyplot as plt

# Physics parameters
mu = 0.1
sigma = 1.0
La = 12000  # rho * sigma * D / mu^2

data = np.genfromtxt('ekin.csv', delimiter=',', comments='#')
t = data[:, 0]
Ekin = data[:, 1]
enst = data[:, 2]

# Spurious capillary number estimate from kinetic energy.
# Ekin = 0.5 * <|u|^2>, so u_rms = sqrt(2*Ekin)
u_rms = np.sqrt(2.0 * Ekin)
Ca_rms = mu * u_rms / sigma

label = f'La={La}, $\\sigma$={sigma}, $\\mu$={mu}'

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

axes[0].plot(t, Ekin, 'o-', label=label)
axes[0].set_xlabel('time $t$')
axes[0].set_ylabel('kinetic energy $E_{kin}$')
axes[0].set_title('Kinetic Energy vs Time')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(t, enst, 'o-', label=label)
axes[1].set_xlabel('time $t$')
axes[1].set_ylabel('enstrophy $En$')
axes[1].set_title('Enstrophy vs Time')
axes[1].legend()
axes[1].grid(True)

axes[2].plot(t, Ca_rms, 'o-', label=label)
axes[2].set_xlabel('time $t$')
axes[2].set_ylabel(r'$Ca^*_{rms} = \mu\, u_{rms} / \sigma$')
axes[2].set_title('Spurious Capillary Number (rms estimate)')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.savefig('plot_stat.png', dpi=150)
plt.show()
