RM := rm -rf
OUTPUT_FILE := run
PDE_OUTPUT_FILE := pde_run

FC = gfortran -fdefault-real-8
#FC = ifort -openmp -r8

OBJS += \
./setting.o \
./main.o

PDE_OBJS += \
./pde_setting.o \
./pde_solver.o

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f90
	@echo 'Building file: $<'
	@echo 'Invoking Fortran Compiler'
	$(FC) -g -O0 -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: md pde

# Tool invocations
md: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Output file name:'
	@echo $(OUTPUT_FILE)
	@echo 'Invoking Fortran Linker'
	$(FC) -o $(OUTPUT_FILE) $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

# Tool invocations
pde: $(PDE_OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Output file name:'
	@echo $(PDE_OUTPUT_FILE)
	@echo 'Invoking Fortran Linker'
	$(FC) -o $(PDE_OUTPUT_FILE) $(PDE_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(PDE_OBJS) $(PDE_OUTPUT_FILE) $(OBJS) $(OUTPUT_FILE) *.mod *~ out/*.dat
	-@echo ' '
