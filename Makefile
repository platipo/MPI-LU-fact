CC=gcc
MPI=mpicc
PAR=LU_mpi.c
EX_PAR=par.out
SEQ=LU_seq.c
EX_SEQ=seq.out
SAMPLE=10

default: sequential mpi

sequential: $(SEQ)
	$(CC) $(SEQ) -o $(EX_SEQ)
	./$(EX_SEQ) $(SAMPLE)

mpi: $(PAR)
	$(MPI) $(PAR) -o $(EX_PAR)
	mpirun -np 4 ./$(EX_PAR) $(SAMPLE)

clean:
	$(RM) -f *.out
