RM := rm -rf

FC = ifort
#FC = gfortran


OBJS_HYBRID += \
../../lib/nrtype.o \
../../lib/random.o \
./paras.o \
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
all: hybrid

# Tool invocations
hybrid: $(OBJS_HYBRID)
	@echo 'Building target: $@'
	@echo 'Invoking Fortran Linker'
	$(FC) -o "hybrid" $(OBJS_HYBRID)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(OBJS_HYBRID) hybrid $(OUTPUT_FILE) *.mod *.o *~
	-@echo ' '
