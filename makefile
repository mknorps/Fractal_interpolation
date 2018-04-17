OBJECTS=affine.o fractal.o simple_test.o

run: $(OBJECTS) 
	gfortran $(OBJECTS) -O3 -o fractal.exe

.f.o:
	gfortran -c $<
clean:
	rm *.o
