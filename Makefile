CPP = /usr/local/opt/llvm/bin/clang
CPPFLAGS = -I/usr/local/opt/llvm/include -fopenmp -I/Library/gurobi811/mac64/include/
LDFLAGS = -L/usr/local/opt/llvm/lib -L/Library/gurobi811/mac64/ -lgurobi81

lsq: ./src/lsq.c
	$(CPP) $(CPPFLAGS) $^ -o $@.o $(LDFLAGS)
dual_decomp_lsq: ./src/dual_decomp_lsq.c
	$(CPP) $(CPPFLAGS) $^ -o $@.o $(LDFLAGS)
dual_ascent: ./src/dual_ascent.c
	$(CPP) $(CPPFLAGS) $^ -o $@.o $(LDFLAGS)
clean:
	rm *.o *.lp *.log