import os
import time
import json

# Use 4000 stars, 200 steps in sources


ST_source = 'src/ST_src.f90'
ST_flags = '-O3'
ST_compile = f'gfortran {ST_flags} -o bin/ST_out.o {ST_source}'


MT_source = 'src/MT_src.f90'
MT_flags = '-O3'
MT_compile = f'gfortran {MT_flags} -fopenmp -o bin/MT_out.o {MT_source}'

# Moved to GPU folder
# GPU_source = 'src/GPU_src.f90'
# GPU_flags = '-O3 -acc '
# GPU_compile = f'pgfortran {GPU_flags} -o bin/GPU_out.o {GPU_source}'


Thread_sweep = [i+1 for i in range(48)]

foo = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Schedule_chunks = foo + [i*10 for i in foo] + [i*100 for i in foo] + [2000]

Schedule_mode = ["static", "dynamic", "guided"]

# Result dictionnary
results = {
    "Thread_sweep" : Thread_sweep,
    "Schedule_chunks": Schedule_chunks,
    "Schedule_mode": Schedule_mode,

    "ST_compile" : ST_compile,
    "MT_compile" : MT_compile,
    # "GPU_compile": GPU_compile
    }



print(f'Benching with the following :\nSchedules : {Schedule_mode}\nChunks    : {Schedule_chunks}\nThreads   : {Thread_sweep}')


# Generate Starting data
from init_files import generate_data
generate_data(4000, 1, 1, 'positions.dat', 'velocities.dat')





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




# Single Threaded Bench
print(ST_compile)
os.system(ST_compile)
ST_time = Bench('bin/ST_out.o', output_null=True)
print(f'ST bench : {ST_time}')
results["ST"] = ST_time


# MT Bench
print(MT_compile)
os.system(MT_compile)

bench_dict = {}
for threads in Thread_sweep:
    print(f' #######  THREADS  {threads}  ##############################')
    thread_dict = {}
    for mode in Schedule_mode:
        print(f' =======  MODE  {mode}  ======================')
        mode_dict = {}
        for chunk in Schedule_chunks:
            print(f' - - Chunk {chunk} - - - -')
            env_var = [
                ["OMP_NUM_THREADS", str(threads)],
                ["OMP_SCHEDULE", f'{mode},{chunk}']
            ]

            t = Bench('bin/MT_out.o',env_var=env_var, output_null=True)
            print(t)
            mode_dict[str(chunk)] = t
        
        thread_dict[mode] = mode_dict
    bench_dict[str(threads)]=thread_dict

results['MT'] = bench_dict






with open('results.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)





