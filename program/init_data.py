import numpy as np

radius = 1

n_stars = 400
tallness = 1 #galaxy sim shennanigans

def generate_points_sphere(n_points):
    points = []
    while len(points) < n_points:
        x, y, z = np.random.uniform(-radius, radius, 3)
        if x**2 + y**2 + z**2 <= radius**2:
            points.append([x, y, z])
    return points

def generate_velocity(pos_data, omega_z=1):
    x, y, z = pos_data[:, 0], pos_data[:, 1], pos_data[:, 2]
    vx = -omega_z * y
    vy = omega_z * x
    vz = np.zeros_like(z)  
    velocities = np.column_stack((vx, vy, vz))
    return velocities

pos_data = np.array(generate_points_sphere(n_stars))
pos_data.T.tofile('data/positions.dat')
vel_data = generate_velocity(pos_data)
vel_data.T.tofile('data/velocities.dat')


# masses
np.array([1/n_stars for _ in range(n_stars)]).T.tofile('data/masses.dat')

quit()



# import numpy as np
# from skyfield.api import load

# # Load ephemeris data for planetary positions and velocities
# eph = load('de421.bsp')  # Or 'de430t.bsp' for higher precision
# ts = load.timescale()
# t = ts.now()

# # List of planets
# planet_names = ['mercury', 'venus', 'earth', 'mars', 'jupiter barycenter', 'saturn barycenter', 'uranus barycenter', 'neptune barycenter']

# # Fetch positions and velocities relative to the Sun
# positions = []
# velocities = []
# for name in planet_names:
#     planet = eph[name]
#     # Get position and velocity vectors relative to the Sun in the ICRF frame
#     pos, vel = planet.at(t).observe(eph['sun']).position.km, planet.at(t).observe(eph['sun']).velocity.km_per_s
#     positions.append(pos)
#     velocities.append(vel)

# # Convert to numpy arrays
# pos_data = np.array(positions)  # Shape: (n_planets, 3)
# vel_data = np.array(velocities)  # Shape: (n_planets, 3)

# # Planetary masses (approximate values in units of solar mass)
# # These are normalized to sum to 1 for the solar system context
# planetary_masses = np.array([
#     3.3011e23,  # Mercury
#     4.8675e24,  # Venus
#     5.97237e24, # Earth
#     6.4171e23,  # Mars
#     1.8982e27,  # Jupiter
#     5.6834e26,  # Saturn
#     8.6810e25,  # Uranus
#     1.02413e26  # Neptune
# ])

# normalization = 1e9
# pos_data = np.array([i/normalization for i in pos_data])
# vel_data = np.array([i/normalization for i in vel_data])
# planetary_masses = np.array([i/sum(planetary_masses) for i in planetary_masses])


# print(pos_data)
# print(vel_data)


# velocity_conversion_factor = 1/86400 # number of seconds in a day
velocity_conversion_factor = 1/86400 # number of seconds in a day
import numpy as np
# Saturn, Mimas, Tethys, Titan
pos_data = np.array([ [0,0,0], [-1.781644589918471E+05, 4.978604034734525E+04, 5.015471754576375E+03], [2.874974592673592E+05, 6.446905997181479E+04, 5.205343659123462E+03], [-7.931453936081687E+05, 9.251201308131195E+05, -7.317676996395952E+03] ])
# vel_data = np.array([ [0,0,0], [-3.101983306920243E+05, -1.200735430939711E+06, 4.624084160095980E+03], [-2.142954962270242E+05, 9.570034498442297E+05, -7.058677519081212E+03], [-3.573469027363011E+05, -3.245883957332915E+05, -1.661041185020754E+03] ]) * velocity_conversion_factor
vel_data = np.array([ [0,0,0], [0,0,0], [0,0,0], [0,0,0] ]) * velocity_conversion_factor
# planetary_masses = np.array([37931206.234, 2.503489, 41.21, 8978.14 ])
planetary_masses = np.array([379312060000.234, 2.503489, 41.21, 8978.14 ])



# pos_data = np.array([ [-1.0,0, 0], [1.0,0, 0], [0.0, 0.0, 0.0]])
# vel_data = np.array([ [1,0.5,0], [-1,-0.5,0], [0.0, 0.0, 0.0]])
# planetary_masses = np.array([1.0, 1.0, 2.0])

print(pos_data)
print(vel_data)



# Save to files
pos_data.T.tofile('data/positions.dat')  # Transposed to save in column-major order
vel_data.T.tofile('data/velocities.dat')
planetary_masses.T.tofile('data/masses.dat')

