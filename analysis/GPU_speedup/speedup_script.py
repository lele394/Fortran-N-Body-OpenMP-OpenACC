import os
import time
import json

# Use 4000 stars, 200 steps in sources


ST_source = 'src/MT_src.f90'
ST_flags = '-O3 -fopenmp'
ST_compile = f'gfortran {ST_flags} -o bin/ST_out.o {ST_source}'


MT_source = 'src/MT_src.f90'
MT_flags = '-O3'
MT_compile = f'gfortran {MT_flags} -fopenmp -o bin/MT_out.o {MT_source}'

GPU_source = 'src/GPU_src.f90'
GPU_flags = '-O3 -acc -gpu=cc75'
GPU_compile = f'pgfortran {GPU_flags} -o bin/GPU_out.o {GPU_source}'



# Result dictionnary
results = {}



# Generate Starting data
from init_files import generate_data
generate_data(10000, 1, 1, 'positions.dat', 'velocities.dat')





def Bench(exec, env_var=[], output_null=False):
    if output_null:
        print(" /!\ outputting to /dev/null /!\ ")
        exec +=  ' >> /dev/null'

    for i in env_var:
        os.environ[i[0]] = i[1]

    print(exec)
    t = time.time()
    os.system(exec)
    return  time.time() - t




# Single Threaded Bench using no o
print(ST_compile)
os.system(ST_compile)

env_var = [
    ["OMP_NUM_THREADS", "1"]
]

ST_time = Bench('bin/ST_out.o', output_null=False, env_var=env_var)
print(f'ST bench : {ST_time}')
results["ST"] = ST_time






# MT Bench  23 threads, dynamic,50
print(MT_compile)
os.system(MT_compile)

env_var = [
    ["OMP_NUM_THREADS", "23"],
    ["OMP_SCHEDULE", f'dynamic,50']
]
MT_time = Bench('bin/MT_out.o', output_null=False, env_var=env_var)
print(f'MT bench : {MT_time}')
results["MT_23"] = MT_time






# MT Bench  48 threads, dynamic,50
print(MT_compile)
os.system(MT_compile)

env_var = [
    ["OMP_NUM_THREADS", "48"],
    ["OMP_SCHEDULE", f'dynamic,50']
]
MT_time = Bench('bin/MT_out.o', output_null=False, env_var=env_var)
print(f'MT bench : {MT_time}')
results["MT_48"] = MT_time



# GPU Bench  
print(GPU_compile)
os.system(GPU_compile)

env_var = [
    ["ACC_DEVICE_NUM", "0"],
    ["CUDA_VISIBLE_DEVICES", '0']
]
GPU_time = Bench('bin/GPU_out.o', output_null=False, env_var=env_var)
print(f'GPU bench : {GPU_time}')
results["GPU"] = GPU_time




with open('results.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)





