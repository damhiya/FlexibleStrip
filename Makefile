default: bin/FlexibleStrip


# obj, mod
obj/FlexibleStrip.o: src/FlexibleStrip.f90 mod/solver.mod mod/util.mod
	gfortran -c src/FlexibleStrip.f90 -O3 -I./mod -o obj/FlexibleStrip.o

obj/solver.o: src/solver.f90 mod/util.mod
	gfortran -c src/solver.f90 -O3 -J./mod -I./mod -o obj/solver.o

mod/solver.mod: obj/solver.o

obj/util.o: src/util.f90
	gfortran -c src/util.f90 -O3 -J./mod -o obj/util.o

mod/util.mod: obj/util.o


# bin
bin/FlexibleStrip: obj/FlexibleStrip.o obj/solver.o obj/util.o
	gfortran obj/FlexibleStrip.o obj/solver.o obj/util.o -O3 -o bin/FlexibleStrip


# else
clean:
	-rm obj/*.o
	-rm mod/*.mod
	-rm bin/FlexibleStrip
	-rm result/result*

plot0: bin/FlexibleStrip
	bin/FlexibleStrip < data/data0 > result/result0
	gnuplot plot0 -p

plot1: bin/FlexibleStrip
	bin/FlexibleStrip < data/data1 > result/result1
	gnuplot plot1 -p

plot2: bin/FlexibleStrip
	bin/FlexibleStrip < data/data2 > result/result2
	gnuplot plot2 -p

plot3: bin/FlexibleStrip
	bin/FlexibleStrip < data/data3 > result/result3
	gnuplot plot3 -p
