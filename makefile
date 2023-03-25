NTL = -lntl -lgmp
OBJ = squfof.o IsSquare.o ZZlib.o morbri.o

main: main.o $(OBJ)
	g++ main.o $(OBJ) $(NTL)
fig4: fig4.o mpqs.o rho.o $(OBJ)
	g++ fig4.o mpqs.o rho.o $(OBJ) $(NTL)
