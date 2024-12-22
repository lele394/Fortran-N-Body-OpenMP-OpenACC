import importlib.util
import subprocess
import os




# === Default flags
FFLAGS = ''
FFLAGS = '-O3'
OMP_THREADS = str(os.cpu_count()-2)
OMP_SCHEDULE= "static,50"
ACC_DEVICE_NUM="0"
CUDA_VISIBLE_DEVICES="0"



# === src files
src_cpu = 'src/MT_OpenMPsim.f90'
src_gpu = 'src/GPU_OpenACCsim.f90'

# === hidden var, set in .f90 directly
hidden_nstars = 3000



os.system("clear")

def is_installed(package_name):
    """Check if a package is installed."""
    return importlib.util.find_spec(package_name) is not None

# ANSI color codes
VYAN = "\033[96m"
MAGENTA = "\033[95m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"
CHECK = "\u2714"  # ✔
CROSS = "\u2718"  # ✘


print(MAGENTA, " === Python dependencies check ===", RESET)

# Python packages checks 
packages = ["numpy", "matplotlib"]

deps_met = True

for pkg in packages:
    if is_installed(pkg):
        print(f"{GREEN}{CHECK} {pkg} {RESET}")
    else:
        print(f"{RED}{CROSS} {pkg} {RESET}")
        deps_met = False



print(MAGENTA, " === Compilers check ===", RESET)

# compilers check
def check_compiler(compiler):
    try:
        subprocess.run(["which", compiler], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False



compiler = "gfortran"
if check_compiler(compiler):
    print(f"{GREEN}{CHECK} {compiler} \n     => CPU and OpenMP compilation is available{RESET}")
else:
    print(f"{RED}{CROSS} {compiler} \n     => CPU and OpenMP compilation is NOT available. THIS IS BLOCKING.{RESET}")
    deps_met = True


is_gpu_available = False
compiler = "pgfortran"
if check_compiler(compiler):
    print(f"{GREEN}{CHECK} {compiler} \n     => Legacy GPU compilation is available{RESET}")
    is_gpu_available = True
else:
    print(f"{RED}{CROSS} {compiler} \n     => Legacy GPU compilation is NOT available{RESET}")

# compiler = "nvfortran"
# if check_compiler(compiler):
#     print(f"{GREEN}{CHECK} {compiler} \n     => Untested GPU compilation is available{RESET}")
#     is_gpu_available = True
# else:
#     print(f"{RED}{CROSS} {compiler} \n     => Untested GPU compilation is NOT available{RESET}")
# print(YELLOW, """
#  => nvfortran is non blocking as it is not used.  <=
#  => You will need to change the compiler files if <=
#  =>        you do not have pgfortran              <=""")
print(YELLOW," ====> GPU dpendencies are non-blocking here.\n\n\n", RESET)




if deps_met:
    print(GREEN, " All dependencies seem to be installed!")
else:
    print(RED, "  You are missing essential dependencies  \n  Please install the missing deps or modify main.py to skip deps check.", RESET)
    quit()





ansi_art = """
\033[49m       \033[38;5;15;49m▄▄▄\033[48;5;15m  \033[38;5;15;49m▄▄▄▄\033[49m                                            \033[m
\033[49m    \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m        \033[49;38;5;15m▀▀\033[38;5;15;48;5;15m▄\033[38;5;15;49m▄▄\033[49m                               \033[38;5;33;49m▄\033[48;5;33m  \033[38;5;33;49m▄\033[49m     \033[m
\033[49m   \033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m              \033[49;38;5;15m▀▀\033[38;5;15;49m▄\033[49m                    \033[38;5;15;49m▄▄\033[48;5;15m \033[49;38;5;15m▀▀▀▀▀\033[48;5;33m      \033[49m    \033[m
\033[49m \033[38;5;15;49m▄\033[48;5;15m \033[49m                   \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m             \033[38;5;15;49m▄▄\033[49;38;5;15m▀▀\033[49m        \033[49;38;5;33m▀\033[48;5;33m    \033[38;5;15;48;5;33m▄\033[38;5;15;49m▄\033[49m   \033[m
\033[49m \033[48;5;15m \033[49m                       \033[49;38;5;15m▀\033[38;5;15;49m▄▄\033[49m       \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m                  \033[49;38;5;15m▀\033[48;5;15m \033[49m  \033[m
\033[49m \033[48;5;15m \033[49m                         \033[49;38;5;15m▀\033[38;5;9;48;5;202m▄\033[48;5;9m    \033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m                       \033[48;5;15m \033[49m \033[m
\033[49m \033[48;5;15m \033[49m                         \033[49;38;5;9m▀\033[48;5;9m     \033[49;38;5;9m▀\033[49m                        \033[48;5;15m \033[38;5;15;49m▄\033[m
\033[49m \033[48;5;15m \033[38;5;15;49m▄\033[49m                       \033[38;5;15;49m▄\033[48;5;15m \033[49;38;5;9m▀▀\033[38;5;9;48;5;9m▄\033[49;38;5;9m▀▀\033[48;5;15m \033[38;5;15;49m▄\033[49m                       \033[48;5;15m \033[49;38;5;15m▀\033[m
\033[49m  \033[49;38;5;15m▀\033[38;5;15;49m▄\033[49m                   \033[38;5;15;49m▄\033[48;5;15m \033[49;38;5;15m▀\033[49m         \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m                    \033[48;5;15m \033[49m \033[m
\033[49m    \033[38;5;47;48;5;15m▄\033[48;5;47m    \033[38;5;47;49m▄\033[49m         \033[38;5;15;49m▄▄\033[49;38;5;15m▀▀\033[49m               \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m               \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49m  \033[m
\033[49m    \033[48;5;47m      \033[38;5;15;49m▄▄▄▄▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀▀▀\033[49m                      \033[49;38;5;15m▀▀\033[38;5;15;49m▄▄▄\033[49m      \033[38;5;15;49m▄▄▄\033[49;38;5;15m▀▀\033[49m   \033[m
\033[49m     \033[49;38;5;47m▀▀▀▀\033[49m                                     \033[49;38;5;15m▀▀▀▀▀▀▀\033[49m       \033[m
\033[49m                                                            \033[m

\033[96m        -====-      N-Body Simulation      -====- \033[0m
  \033[0m                 Author - Léo BECHET
                  Licence - CC-BY-NC-SA
  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \033[0m
"""
print(ansi_art)





def SetDefaultSettings():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    FFLAGS = '-O3'
    OMP_THREADS = str(os.cpu_count()-2)
    OMP_SCHEDULE= "static,50"
    ACC_DEVICE_NUM="0"
    CUDA_VISIBLE_DEVICES="0"


def ShowEnvVar():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    print(f' FFLAGS _____________ : {FFLAGS}')
    print(f' OMP_THREADS          : {OMP_THREADS}')
    print(f' OMP_SCHEDULE _______ : {OMP_SCHEDULE}')
    print(f' ACC_DEVICE_NUM       : {ACC_DEVICE_NUM}')
    print(f' CUDA_VISIBLE_DEVICES : {CUDA_VISIBLE_DEVICES}')
    Menu()

def EditEnvVar():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    print(f' 0) Cancel')
    print(f' 1) FFLAGS _____________ : {FFLAGS}')
    print(f' 2) OMP_THREADS          : {OMP_THREADS}')
    print(f' 3) OMP_SCHEDULE _______ : {OMP_SCHEDULE}')
    print(f' 4) ACC_DEVICE_NUM       : {ACC_DEVICE_NUM}')
    print(f' 5) CUDA_VISIBLE_DEVICES : {CUDA_VISIBLE_DEVICES}')

    PASS = True
    while PASS:
        PASS = False
        inp = int(input(f'\n Which Env Var do you wish to edit? (number)\n > '))

        match inp:
            case 1:
                FFLAGS=str(input("value > "))
            case 2:
                OMP_THREADS=str(input("value > "))
            case 3:
                OMP_SCHEDULE=str(input("value > "))
            case 4:
                ACC_DEVICE_NUM=str(input("value > "))
            case 5:
                CUDA_VISIBLE_DEVICES=str(input("value > "))

            case _:
                PASS = True
                print(" ID not recognised")
    print("\n")
    Menu()



import numpy as np
def generate_points_sphere(n_points, radius=1.0, offset=(0,0,0)):
    n_points = int(n_points)
    points = []
    while len(points) < n_points:
        x, y, z = np.random.uniform(-radius, radius, 3)
        if x**2 + y**2 + z**2 <= radius**2:
            points.append([x+offset[0], y+offset[1], z+offset[2]])
    return points

def generate_velocity(pos_data, omega_z=1):
    x, y, z = pos_data[:, 0], pos_data[:, 1], pos_data[:, 2]
    
    vx = -omega_z * y
    vy = omega_z * x
    vz = np.zeros_like(z)  
    
    velocities = np.column_stack((vx, vy, vz))
    
    return velocities


def CompileAndExec_CPU():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES, N_STARS

    print(YELLOW, " Generating test data", RESET)
    pos_data = np.array(generate_points_sphere(hidden_nstars, 1.0, offset=(0, 0.25, 0.25)))
    pos_data.T.tofile('data/positions.dat')

    vel_data = generate_velocity(pos_data)
    vel_data.T.tofile('data/velocities.dat')

    print(" Compiling program")
    i = os.system(f'gfortran -fopenmp {FFLAGS} -o bin/CPU.o {src_cpu}')

    if i == 0 :
        os.environ['OMP_SCHEDULE'] = OMP_SCHEDULE
        os.environ['OMP_NUMT_THREADS'] = OMP_THREADS
        os.system("bin/CPU.o")
        plot_energy()
        anim_realTime('out_data/position.dat')
    else:
        print(RED, "Error during compilation, aborting.", RESET)

    Menu()









def CompileAndExec_GPU():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES, N_STARS
    print(YELLOW, "  This program does not check CUDA capabilities on your machine and assume it is installed and configured")
    print("  Please specify your GPU architecture.")
    print("  With XXX as your GPU architecture. For QUADRO RTX 5000, use cc75 (default, enter to use)", RESET)
    GPU_arch = input(MAGENTA+"default=cc75 > "+RESET)
    if GPU_arch == "" : GPU_arch="cc75"

    print(YELLOW, " Generating test data", RESET)
    pos_data = np.array(generate_points_sphere(hidden_nstars, 1.0, offset=(0, 0.25, 0.25)))
    pos_data.T.tofile('data/positions.dat')
    vel_data = generate_velocity(pos_data)
    vel_data.T.tofile('data/velocities.dat')
    print(YELLOW, " === Compiling program ===", RESET)
    i = os.system(f'pgfortran -acc {FFLAGS} -gpu={GPU_arch} -o bin/GPU.o {src_gpu}') # Modify nvfortran here if needed

    if i == 0 :
        os.environ['ACC_DEVICE_NUM'] = ACC_DEVICE_NUM
        os.environ['CUDA_VISIBLE_DEVICES'] = CUDA_VISIBLE_DEVICES
        os.system("bin/GPU.o")
        anim_realTime('out_data/position.dat')

    else:
        print(RED, "Error during compilation, aborting.", RESET)

    Menu()













import matplotlib.pyplot as plt
import matplotlib.animation

def plot_energy():
    data = np.loadtxt('out_data/energy.dat')
    potential_energy = data[:, 0]
    kinetic_energy = data[:, 1]
    total_energy = kinetic_energy + potential_energy
    plt.figure(figsize=(10, 6))
    plt.plot(kinetic_energy, label='Kinetic Energy', color='blue')
    plt.plot(potential_energy, label='Potential Energy', color='red')
    plt.plot(total_energy, label='Total Energy', color='green', linestyle='--')
    plt.xlabel('Time Step (T)')
    plt.ylabel('Energy')
    plt.title('Kinetic, Potential, and Total Energy vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()





















def load_file(filename: str, nbline: int = 0, timestep: int = -1, param: list = []):
    with open(filename, "rb") as file:
        if timestep == -1:
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            metadata = np.fromfile(file, dtype=np.int32, count=3)    # Read metadata
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            return metadata
        else:
            if param is None or len(param) < 3:
                raise ValueError("Parameter list 'param' must include [num_particles, num_timesteps, float_precision].")
            num_particles, num_timesteps, float_precision = param
            dtype = np.float32 if float_precision == 4 else np.float64
            record_bytes = 4 + (num_particles * 3 * dtype().nbytes) + 4
            file.seek(4 + (record_bytes * timestep), 0)  # 4 bytes for record size + N records
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            data = np.fromfile(file, dtype=dtype, count=num_particles * 3)  # Read positions
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            data = data.reshape(3, nbline)
            return data

def anim_realTime(file: str):
    figs = (720, 720)
    def update_graph(t: int):
        data = load_file(file, param[0], t, param)
        graph._offsets3d = (data[0, :], data[1, :], data[2, :])
    param = load_file(file)
    data = load_file(file, param[0], 0, param)
    dpi_v = 100
    sfig = (figs[0]/dpi_v, figs[1]/dpi_v)
    fig = plt.figure(figsize = sfig, facecolor = 'k', num = 1)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.axes.set_xlim3d(left=-3, right=3) 
    ax.axes.set_ylim3d(bottom=-3, top=3) 
    ax.axes.set_zlim3d(bottom=-3, top=3)
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
    graph = ax.scatter(data[0, :], data[1, :], data[2, :], s=2, c='w', linewidths=0)
    ani = matplotlib.animation.FuncAnimation(fig, update_graph, param[1],
                            interval=1, blit=False)
    plt.show()









def Menu():
    print(YELLOW, "\n ======   M E N U   ======", RESET)
    print(f'  0) Quit')
    print(f'  1) Show Env Var')
    print(f'  2) Edit Env Var')
    print(GREEN, f' 3) Compile and exec CPU version', RESET)

    if is_gpu_available:
        print(GREEN,f' 4) Compile and exec GPU version',RESET)
    else:
        print(RED, f' 4) Compile and exec GPU version [UNAVAILABLE]', RESET)

    print(f'  5) Reset Env Var')
    

    print("\n")
    
    PASS = True
    while PASS:
        inp = int(input(" > "))
        match inp:
            case 0:
                quit()
            case 1:
                ShowEnvVar()
            case 2:
                EditEnvVar()
                ShowEnvVar()
            case 3:
                CompileAndExec_CPU()
            case 4:
                CompileAndExec_GPU()
            case 5:
                SetDefaultSettings()
                ShowEnvVar()

            case _:
                PASS = True
                print(" ID not recognised")


ShowEnvVar()

Menu()















