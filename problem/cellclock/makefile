RM := rm -rf

FC = ifort
OBJS += \
../../lib/nrtype.o \
../../lib/random.o \
./setting.o \
./main.o 

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking Fortran Compiler'
	$(FC) -g -O0 -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: main

# Tool invocations
main: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking Fortran Linker'
	$(FC) -o "run" $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) run *.mod *.o *~
	-@echo ' '

