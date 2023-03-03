# Makefile for compiling the Serial Spectral Code Version 2.1 
# This uses the Intel Fortran Compiler, with FFTW 
# ---------------------------------------------------------------------

compile = mpif90 -g -fopenmp -O1 -fno-strict-aliasing -Wno-lto-type-mismatch
link = -L/home/shubham/Dependencies/fftw3_install/lib -lfftw3_mpi -lfftw3 -lfftw3_threads  
incl = -I/home/shubham/Dependencies/fftw3_install/include


#-----------------------------------------------------------------------

#OBJ_FILE=spectral_serial.o  global_array.o  \
#energy_serial.o rnkt_serial.o force_serial.o adbsh_serial.o nonlin_serial.o \
#energy_injection.o dermat_serial.o  omplib.o decimation.o

OBJ_FILE=spectral_serial.o  global_array.o energy_serial.o rnkt_serial.o force_serial.o nonlin_serial.o energy_injection.o dermat_serial.o  omplib.o decimation.o


MOD_FILE=mod_serial_fluid.o
MOD_FILE1=mod_initial.o
MOD_FILE2=mod_particle.o
#test : 
#	echo $(OBJ_FILE) 
# ---------------------------------------------------------------------
part.exe:$(OBJ_FILE) $(MOD_FILE) $(MOD_FILE1)  
	$(compile) $(OBJ_FILE) $(MOD_FILE) $(MOD_FILE1)  $(link) -o part.exe 
#                 ------------------------------------               
spectral_serial.o: spectral_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
	$(compile) $(link) -c -o spectral_serial.o \
	spectral_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
#         
global_array.o: global_array.f90 $(MOD_FILE) 
	$(compile) -c -o global_array.o \
	global_array.f90 $(MOD_FILE)
#         
#fourier.o: fourier.f90 $(MOD_FILE) 
#	$(compile) -c -o fourier.o \
#fourier.f90 -lm $(MOD_FILE) $(link) $(incl)
#         

energy_serial.o: energy_serial.f90 $(MOD_FILE) 
	$(compile) -c -o energy_serial.o \
	energy_serial.f90 $(MOD_FILE)
#         
adbsh_serial.o: adbsh_serial.f90 $(MOD_FILE) 
	$(compile) -c -o adbsh_serial.o \
	adbsh_serial.f90 $(MOD_FILE)
#         
nonlin_serial.o: nonlin_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
	$(compile) -c -o nonlin_serial.o  \
	nonlin_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
#         
force_serial.o: force_serial.f90 $(MOD_FILE) 
	$(compile) -c -o force_serial.o \
	force_serial.f90 $(MOD_FILE)
#
energy_injection.o: energy_injection.f90 
	$(compile) -c -o energy_injection.o \
	energy_injection.f90 $(MOD_FILE)

#         
rnkt_serial.o: rnkt_serial.f90 $(MOD_FILE) 
	$(compile) -c -o rnkt_serial.o \
	rnkt_serial.f90 $(MOD_FILE)
#         
dermat_serial.o: dermat_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
	$(compile) -c -o dermat_serial.o \
	dermat_serial.f90 $(MOD_FILE) $(MOD_FILE1) 
#        
mod_serial_fluid.o:  mod_serial_fluid.f90 
	$(compile) -c -o mod_serial_fluid.o mod_serial_fluid.f90 $(incl)
# ------------------------------------------------------------------------ 
mod_initial.o:  mod_initial.f90 $(MOD_FILE) omplib.o
	$(compile) -c -o mod_initial.o mod_initial.f90 -lm $(MOD_FILE) omplib.o $(link) $(incl)
# ------------------------------------------------------------------------ 
#mod_particle.o:  mod_particle.f90 $(MOD_FILE)
#	$(compile) -c mod_particle.f90 -o mod_particle.o -lm $(MOD_FILE) $(link) $(incl)  
# ------------------------------------------------------------------------ 
omplib.o:  omplib.f90 $(MOD_FILE) 
	$(compile) -c omplib.f90  $(MOD_FILE) -o omplib.o
# ------------------------------------------------------------------------ 
decimation.o: decimation.f90 $(MOD_FILE) 
	$(compile) -c -o decimation.o \
	decimation.f90 $(MOD_FILE)
#-------------------------------------------         
clean:
	rm -f *.o *.mod *.exe
total_clean:
	rm -f *.o *.mod core *.exe *.out
clean_out:
	rm -f *.out
## --------------------------------------------------------------------##  
