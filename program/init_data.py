import numpy as np

# Masses of planets
m_1 = 10.0
m_2 = 20.0
m_3 = 30.0

# Starting coordinates for planets
p1_start = np.array([-10.0, 10.0, -11.0])  # Position of planet 1
v1_start = np.array([-3.0, 0.0, 0.0])      # Velocity of planet 1

p2_start = np.array([0.0, 0.0, 0.0])       # Position of planet 2
v2_start = np.array([0.0, 0.0, 0.0])       # Velocity of planet 2

p3_start = np.array([10.0, 10.0, 12.0])    # Position of planet 3
v3_start = np.array([3.0, 0.0, 0.0])       # Velocity of planet 3

# Create position and velocity arrays
r = np.array([p1_start, p2_start, p3_start])  # Positions of the planets
v = np.array([v1_start, v2_start, v3_start])  # Velocities of the planets

# Masses of the planets
masses = np.array([m_1, m_2, m_3])

# Save data to files (if necessary)
r.T.tofile('data/positions.dat')  # Save positions (transpose for column-major format)
v.T.tofile('data/velocities.dat')  # Save velocities (transpose for column-major format)
masses.T.tofile('data/masses.dat')  # Save masses (transpose for column-major format)
