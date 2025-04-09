executables = x_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o util.o constants.o activation.o random.o arch_fit.o arch_sim.o binomial.o lsq.o student_t.o pdf.o x.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

x_gfort.exe: kind.o util.o constants.o activation.o random.o arch_fit.o arch_sim.o binomial.o lsq.o student_t.o pdf.o x.o
	$(FC) -o x_gfort.exe kind.o util.o constants.o activation.o random.o arch_fit.o arch_sim.o binomial.o lsq.o student_t.o pdf.o x.o $(FFLAGS)

run: $(executables)
	./x_gfort.exe

clean:
	rm -f $(executables) $(obj)

