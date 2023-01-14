#########################################################################
# This makefile should build on linux provided it has access to the X11 library
# which typically requires the development package, a fortran compiler,
# and git. It may work on Mac if XQuartz is present. It won't work on MSWindows.
# Decide the FORTRAN compiler and verify/make the accis graphics routines:
#ACCISPARENT:=$(realpath .)
#export ACCISPARENT
include ACCIS.mk
#########################################################################
LIBRARIES := $(LIBRARIES)
LIBDEPS := $(LIBDEPS)
COMPILE-SWITCHES:=$(COMPILE-SWITCHES) -Wno-unused-dummy-argument -Wno-integer-division
#########################################################################
OBJECTS=BGKint.o Zfun.o

.PRECIOUS : $(OBJECTS)
# Not needed if fortran default is defeated.
#########################################################################
# Patterns for compilation etc.
# Defeat default fortran rule.
% : %.f

% : %.f90

# Defeat the Modula-2 make booby-trap.
% : %.mod

###############################
%.o : %.f makefile ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

% : %.f90  makefile $(ACCISX) $(OBJECTS) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90 $(OBJECTS) $(LIBPATH) $(LIBRARIES)

% : %.f  makefile $(ACCISX) $(OBJECTS) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f $(OBJECTS) $(LIBPATH) $(LIBRARIES)

#########################################################################


