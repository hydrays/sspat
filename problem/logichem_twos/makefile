RM := rm -rf

OUTPUT_FILE := run

# add inputs and outputs from these tool invocations to the build variables 
OBJS += \
../../lib/nrtype.o \
../../lib/random.o \
./setting.o \
./main.o


# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking Fortran Compiler'
	ifort -g -O0 -c -o "$@" "$<"
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
	ifort -o $(OUTPUT_FILE) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) $(OUTPUT_FILE) *.mod *~
	-@echo ' '
