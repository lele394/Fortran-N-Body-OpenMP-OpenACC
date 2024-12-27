import numpy as np


#
# Available config :
# "Lagrange", "2-plane 8"
configuration = "Lagrange"

if configuration == "Lagrange":

    # Parameters
    M_idle = 1.0  # Mass of each body
    R = 1.0  # Distance from each body to the center of mass
    G = 1.0  # Gravitational constant
    v_mag = np.sqrt(G * M_idle / R)  # Orbital velocity for stable configuration

    # Positions: Equilateral triangle in the XY plane
    pos_data = np.array([
        [R, 0.0, 0.0],  # Body 1
        [-R / 2, np.sqrt(3) * R / 2, 0.0],  # Body 2
        [-R / 2, -np.sqrt(3) * R / 2, 0.0]  # Body 3
    ])

    # Velocities: Tangential to circular orbit
    vel_data = np.array([
        [0.0, v_mag, 0.0],  # Body 1
        [-v_mag * np.sqrt(3) / 2, -v_mag / 2, 0.0],  # Body 2
        [v_mag * np.sqrt(3) / 2, -v_mag / 2, 0.0]  # Body 3
    ])

    # Write data to files
    np.array([M_idle, M_idle, M_idle]).T.tofile('data/masses.dat')  # Masses
    pos_data.T.tofile('data/positions.dat')  # Positions (column-major order)
    vel_data.T.tofile('data/velocities.dat')  # Velocities (column-major order)



if configuration == "2-plane 8":

    # Initial conditions based on the data provided
    r = np.zeros((3, 3))  # Positions array (3 bodies, 3 coordinates)
    v = np.zeros((3, 3))  # Velocities array (3 bodies, 3 velocity components)

    # Body 0
    r[0] = [0.9700436, -0.24308753, 0]
    v[0] = [0.466203685, 0.43236573, 0]

    # Body 1 (symmetrical to body 0)
    r[1] = [-r[0][0], -r[0][1], -r[0][2]]
    v[1] = [v[0][0], v[0][1], v[0][2]]

    # Body 2 (center of mass, velocity is -2 * velocity of body 0)
    r[2] = [0, 0, 0]
    v[2] = [-2 * v[0][0], -2 * v[0][1], -2 * v[0][2]]

    # Masses for each body (equal masses)
    masses = np.array([1.0, 1.0, 1.0])

    # Save data to files (if necessary)
    r.T.tofile('data/positions.dat')  # Save positions (transpose for column-major format)
    v.T.tofile('data/velocities.dat')  # Save velocities (transpose for column-major format)
    masses.T.tofile('data/masses.dat')  # Save masses (transpose for column-major format)

    # Output the initial conditions for verification
    print("Positions (r):")
    print(r)
    print("Velocities (v):")
    print(v)
    print("Masses:")
    print(masses)


