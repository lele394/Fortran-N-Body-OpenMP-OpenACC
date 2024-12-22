import numpy as np

def generate_data(n_stars = 50000, radius = 1.0, omega_z = 1.0, pos_file='positions.dat', vel_file='velocities.dat', transpose=True):

    # Gravitational constant in kpc (km/s)^2 M_sun^-1
    G = 4.3e-6

    try:
        def generate_points_sphere(n_points, radius=1.0, offset=(0, 0, 0)):
            points = []
            while len(points) < n_points:
                x, y, z = np.random.uniform(-radius, radius, 3)
                if x**2 + y**2 + z**2 <= radius**2:
                    points.append([x + offset[0], y + offset[1], z + offset[2]])
            return points

        pos_data = np.array(
        generate_points_sphere(n_stars, radius)
        )

        def generate_velocity(pos_data, omega_z=1.0):
            x, y, z = pos_data[:, 0], pos_data[:, 1], pos_data[:, 2]
        
        # Circular velocity components due to rotation
            vx = -omega_z * y
            vy = omega_z * x
            vz = np.zeros_like(z)  
        
            velocities = np.column_stack((vx, vy, vz))
        
            return velocities

        vel_data = generate_velocity(pos_data, omega_z=omega_z)

        # Ensure gravitational binding (optional validation step)
        total_mass = n_stars  # Assume 1 M_sun per star
        velocity_dispersion = np.mean(np.sqrt(vel_data[:, 0]**2 + vel_data[:, 1]**2))
        binding_mass = (velocity_dispersion**2 * radius) / G
        print(f"Total Mass: {total_mass} M_sun, Binding Mass: {binding_mass:.2e} M_sun")
        if total_mass < binding_mass:
            print("Warning: System may not be gravitationally bound!")

        # Save data
        if transpose:
            pos_data.T.tofile(pos_file)
            vel_data.T.tofile(vel_file)
        else:
            pos_data.tofile(pos_file)
            vel_data.tofile(vel_file)
            
        return 0

    except Exception as e:
        raise e

