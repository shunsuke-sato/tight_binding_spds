FC = mpif90 -O2 ## gfotran
#FC = mpif90 -O0 -fbounds-check ## gfotran
#FC = mpif90 -O2 -pg ## gfotran
#FC = mpiifort -O3 -ipo -xHOST  -L$MKL_HOME/lib/intel64  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64 ## draco
#FC = mpiifort -O3 -ipo -xHOST -L/mpcdf/soft/SLES122/common/intel/ps2017.7/17.0/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -Wl,-rpath,/mpcdf/soft/SLES122/common/intel/ps2017.7/17.0/linux/mkl/lib/intel64 ## draco
#FC = mpiifort -O3 -ipo -xHOST ## draco

LN = -llapack -lblas
#LN = 

VPATH = src:object
#SRC = $(shell cd src ;ls *.f90 ;cd ..)
#OBJ = $(SRC:.f90=.o)
#OBJ_dir = $(addprefix object/,$(OBJ))


PROG = tb_model
OBJ = object/math.o object/parallel.o object/communication.o object/constants.o object/inputoutput.o object/electronic_system.o object/main.o

$(PROG):math.o \
        parallel.o \
        communication.o \
        constants.o \
        inputoutput.o \
        electronic_system.o \
        main.o
	$(FC) -o $(PROG) $(OBJ) $(LN)

math.o:math.f90
	$(FC) -c $< $(LN);mv $@  object 

parallel.o:parallel.f90
	$(FC) -c $< $(LN);mv $@  object 

constants.o:constants.f90
	$(FC) -c $< $(LN);mv $@  object 

communication.o:communication.f90 parallel.o
	$(FC) -c $< $(LN);mv $@  object 

inputoutput.o:inputoutput.f90 parallel.o communication.o
	$(FC) -c $< $(LN);mv $@  object 

electronic_system.o:electronic_system.f90 parallel.o communication.o math.o constants.o inputoutput.o
	$(FC) -c $< $(LN);mv $@  object 

main.o:main.f90 parallel.o inputoutput.o electronic_system.o
	$(FC) -c $< $(LN);mv $@  object 

clean:
	rm  -f  object/*.o  *.mod ${PROG}
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod ${PROG} */#* *.out *.log
