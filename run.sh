cd ./program
python init_data.py


# gfortran -fopenmp -O3 -o bin/CPU.o src/MT_OpenMPsim.f90
gfortran -O3 -o bin/CPU.o src/MT_OpenMPsim.f90

bin/CPU.o

python plot.py
