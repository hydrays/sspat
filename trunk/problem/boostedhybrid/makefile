RM := rm -rf

#FC = ifort
FC = gfortran
#SYSTEM := SinMolEnz
#SYSTEM := Kepler_Elston
#SYSTEM := gene_expression
#SYSTEM := schlogl
#SYSTEM := Hybrid_SinMolEnz
#SYSTEM := TwoGeneOscillator -- no scale seperation
#SYSTEM := promoter
#SYSTEM := oscillator
#SYSTEM := hsr
#SYSTEM := fge_modified
#SYSTEM := 0091
SYSTEM := system1
#SYSTEM := system3

SSA_OBJS += \
../../lib/nrtype.o \
../../lib/random.o \
./systems/$(SYSTEM).o \
./ssa.o 

HYBRID_OBJS += \
../../lib/nrtype.o \
../../lib/random.o \
./systems/$(SYSTEM).o \
./hybrid.o 

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking Fortran Compiler'
	$(FC) -g -O0 -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ssa hybrid

# Tool invocations
ssa: $(SSA_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking Fortran Linker'
	$(FC) -o "ssa" $(SSA_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

hybrid: $(HYBRID_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking Fortran Linker'
	$(FC) -o "hybrid" $(HYBRID_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(SSA_OBJS) $(HYBRID_OBJS) ssa hybrid *.mod *.o *~
	-@echo ' '
