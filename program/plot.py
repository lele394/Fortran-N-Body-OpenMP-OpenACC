#!/usr/bin/env python3

# Author: Barnabé DEFORET

speed = int(20)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation

####################################
#             Functions            #
####################################

# def load_file(filename: str, nbline: int = 0, timestep: int = -1, param: list = []):
#     # This function load the data from the position file. It is designed to read data for a specific time step.
#     #
#     # Args:
#     #     filename (str): path to the file
#     #     nbline (int): number of particles in a line
#     #     timestep (int): time step to read
#     #     param, (list): parameters, i.e number of particles, number of time steps, float precision.
#     #
#     # Return:
#     #     (np.array): Array containing the x, y, z positions of the particles

#     if timestep == -1:
#         # Read the first 3 integers of the file (number of particles, number of time steps, float precision)
#         data = np.loadtxt(filename, max_rows=1, dtype=np.int32)
#     else:
#         # Read the data for the time step t
#         data = np.loadtxt(filename, skiprows=1+timestep, max_rows=1, dtype=np.float64)
#         data = data.reshape(3, nbline)
#     return data
















def load_file(filename: str, nbline: int = 0, timestep: int = -1, param: list = []):
    with open(filename, "rb") as file:
        if timestep == -1:
            # Read the first 3 integers of the file (metadata)
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            metadata = np.fromfile(file, dtype=np.int32, count=3)    # Read metadata
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            return metadata
        else:
            # Validate 'param'
            if param is None or len(param) < 3:
                raise ValueError("Parameter list 'param' must include [num_particles, num_timesteps, float_precision].")
            
            num_particles, num_timesteps, float_precision = param
            dtype = np.float32 if float_precision == 4 else np.float64
            record_bytes = 4 + (num_particles * 3 * dtype().nbytes) + 4
            
            # Skip records before the desired timestep
            file.seek(4 + (record_bytes * timestep), 0)  # 4 bytes for record size + N records
            
            # Read the data for the timestep
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            data = np.fromfile(file, dtype=dtype, count=num_particles * 3)  # Read positions
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            
            # Reshape data into (3, nbline) format (x, y, z)
            data = data.reshape(3, nbline)
            return data









size = 3e0

def anim_realTime(file: str, figs: tuple):
    # This function create an animation of the particles in real time.
    #
    # Args:
    #     file (str): path to the file
    #     figs (tuple): size (x, y) of the output GIF (in pixels)


    def update_graph(t: int):
        # Thus function update the position of the particles for time step t
        #
        # Args:
        #     t (int): time step to read

        data = load_file(file, param[0], t, param)
        # print(t)
        # print(data)
        graph._offsets3d = (data[0, :], data[1, :], data[2, :])
        # graph._offsets3d = (data[2, :], data[1, :], data[0, :])

    # Get parameters
    param = load_file(file)

    # Load first time step data
    data = load_file(file, param[0], 0, param)

    # Convert figure size to inches
    dpi_v = 100
    sfig = (figs[0]/dpi_v, figs[1]/dpi_v)

    # Create figure
    fig = plt.figure(figsize = sfig, facecolor = 'k', num = 1)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.axes.set_xlim3d(left=-size, right=size) 
    ax.axes.set_ylim3d(bottom=-size, top=size) 
    ax.axes.set_zlim3d(bottom=-size, top=size)
    ax.tick_params(axis='x', colors='w', which='both')
    ax.tick_params(axis='y', colors='w', which='both')
    ax.tick_params(axis='z', colors='w', which='both')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.yaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.zaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.set_facecolor("k")
    plt.tight_layout(pad=1)

    # Plot data
    # colors = ['y', 'r', 'g', 'c']
    colors = 'w'
    graph = ax.scatter(data[0, :], data[1, :], data[2, :], s=8, c=colors, linewidths=2)

    # Animate
    ani = matplotlib.animation.FuncAnimation(fig, update_graph, param[1],
                            interval=1000/60, blit=False)
                            
    plt.show()

if __name__ == '__main__':

    ####################################
    #            Parameters            #
    ####################################

    # Path to the file
    filepath_pos = 'out_data/position.dat'

    # Size of the output plot (in pixels)
    figs = (720, 720)

    ####################################
    #               Main               #
    ####################################

    anim_realTime(filepath_pos, figs)