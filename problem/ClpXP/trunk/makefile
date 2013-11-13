RM := rm -rf
OUTPUT_FILE := run
OUTPUT_FILE3D := run3d

FC = gfortran -fopenmp -fdefault-real-8
#FC = ifort -openmp -r8

OBJS += \
./random.o \
./setting.o \
./main.o

OBJS3D += \
./random.o \
./setting3d.o \
./main3d.o

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
main: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Output file name:'
	@echo $(OUTPUT_FILE)
	@echo 'Invoking Fortran Linker'
	$(FC) -o $(OUTPUT_FILE) $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

3d: $(OBJS3D)
	@echo 'Building target: $@'
	@echo 'Output file name:'
	@echo $(OUTPUT_FILE3D)
	@echo 'Invoking Fortran Linker'
	$(FC) -o $(OUTPUT_FILE3D) $(OBJS3D)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) $(OBJS3D) $(OUTPUT_FILE) $(OUTPUT_FILE3D) *.mod *~ out/*.dat
	-@echo ' '
