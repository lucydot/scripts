# ------------------------------------------------
# makefile de poche pour programme fortran
# ------------------------------------------------

CC = gfortran
CFLAGS= -g -fbounds-check -fbacktrace -Wall
LFLAGS=
# ------------------------------------------------

test.x: test_module.o crysmod.o
	$(CC) $^ -o $@

test_module.o: crysmod.o

crysmod: crysmod.o
	@echo compilation of $@.f90

%.o: %.f90
	$(CC) -c $(CFLAGS) -o $@ $< 

# ------------------------------------------------
# clean up
.PHONY: clean
clean:
	rm *.o *.mod test.x

# ------------------------------------------------
# $@  	Le nom de la cible
# $<  	Le nom de la première dépendance
# $^  	La liste des dépendances
# ------------------------------------------------

