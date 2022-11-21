FC = mpif90 -O2 ## gfotran
#FC = mpif90 -O0 -fbounds-check ## gfotran
#FC = mpif90 -O2 -pg ## gfotran
#FC = mpiifort -O3 -ipo -xHOST  -L$MKL_HOME/lib/intel64  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64 ## draco
#FC = mpiifort -O3 -ipo -xHOST -L/mpcdf/soft/SLES122/common/intel/ps2017.7/17.0/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -Wl,-rpath,/mpcdf/soft/SLES122/common/intel/ps2017.7/17.0/linux/mkl/lib/intel64 ## draco
#FC = mpiifort -O3 -ipo -xHOST ## draco

#LN = -llapack -lblas
LN = 

VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = tbmodel

$(PROG):math.o \
        parallel.o \
        communication.o \
        constants.o \
        io_mod.o \
        global_variables.o \
        main.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod ${PROG}
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod ${PROG} */#* *.out *.log
