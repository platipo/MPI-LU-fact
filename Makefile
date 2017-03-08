CC=gcc
MPI=mpicc
PAR=LU_mpi.c
EX_PAR=par.out
SEQ=LU_seq.c
EX_SEQ=seq.out
ifdef N
  SAMPLE:=$(N)
else
  SAMPLE=14
endif
PRINT_ALU := $(shell if [[ ${SAMPLE} -le 15 ]]; then echo -DALU; fi)

default: seq mpi

.PHONY: seq mpi clean

$(EX_SEQ): $(SEQ)
	$(CC) $(SEQ) $(PRINT_ALU) -o $(EX_SEQ)

$(EX_PAR): $(PAR)
	$(MPI) $(PAR) $(PRINT_ALU) -o $(EX_PAR)

mpi: $(EX_PAR) 
	mpirun -np 4 ./$(EX_PAR) $(SAMPLE)

seq: $(EX_SEQ)
	./$(EX_SEQ) $(SAMPLE)

clean:
	$(RM) -f *.out
