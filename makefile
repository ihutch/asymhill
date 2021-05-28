
libraries = -L/usr/X11R6/lib/ -L/home/hutch/accis/ -laccisX -lX11 -llapack
FORTRAN=gfortran
COMPILE-SWITCHES = -Wall -O2

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : %.f makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

%.o : %.F makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.F

% : %.f
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f  $(libraries)

% : %.f90
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90  $(libraries)

% : %.F
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.F  $(libraries)


