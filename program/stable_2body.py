import numpy as np

M_idle = 1
pos_data = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
vel_data = np.array([   [0.0,   np.sqrt(M_idle/4)  , 0.0]   , [0.0, -np.sqrt(M_idle/4), 0.0]])
np.array([M_idle, M_idle]).T.tofile('data/masses.dat')


pos_data.T.tofile('data/positions.dat')
vel_data.T.tofile('data/velocities.dat')