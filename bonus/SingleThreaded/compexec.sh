echo "======= Compiling Main Program ======="
gfortran -Wall -O3 sim_ST.f90 -o out.o
gfortran -Wall -O3 sim_ST_optimized.f90 -o out_opt.o
echo done


echo "======= Running Program =============="
time ./out.o
time ./out_opt.o
# python3 plot.py