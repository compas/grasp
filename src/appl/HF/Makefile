.SUFFIXES: .f90 .mod

BIN = ${GRASP}/bin

UTIL = hf

install: EXE 

EXE :        $(BIN)/hf 

OBJ : hf .o

$(BIN)/hf: HF.o
	$(FC) -o $(BIN)/hf  HF.o

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core *.mod



