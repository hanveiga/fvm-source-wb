CF=gfortran

FFLAGS=

OBJS=dg_commons.o legendre.o

exe: $(OBJS) dg_with_source.f90
	$(CF) $(FFLAGS) -o dg dg_with_source.f90 $(OBJS)

dg_commons.o: dg_commons.f90
	$(CF) $(FFLAGS) -c $<

legendre.o: legendre.f90
	$(CF) $(FFLAGS) -c $<

clean:	
	rm -f *.f90~ 
	rm -f *.mod
	rm dg
