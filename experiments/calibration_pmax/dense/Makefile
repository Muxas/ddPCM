# Dense solver
EXE=../main.exe
# Inputs geometries
INPUTS=$(wildcard *.inp)
# Possible values of lmax (model degrees) to check
LMAXS:=6 8 10 12
# Size of Lebedev grids
NGRID:=590
# Outputs
OUTPUTS:=$(foreach lmax, $(LMAXS), $(foreach input, $(INPUTS:.inp=), \
	   $(input).$(lmax).out))

# Standard trick to ignore file "all" to build target "all"
.PHONY: all clean

# All jobs
all:	$(OUTPUTS)

# Set parameters of each dense test
%.out:	ARGV=$(subst ., ,$*)
%.out:	INPUT=$(firstword $(ARGV)).inp
%.out:	LMAX=$(word 2, $(ARGV))
%.out:	OUTPUT=$@

# Recipe to prepare dense result
%.out:
	@echo Running dense experiment with input $(INPUT) lmax=$(LMAX) \
	    ngrid=$(NGRID)
	@rm -f Input.txt
	@echo "0 ! print flag" > Input.txt
	@echo "1 ! number of OpenMP threads" >> Input.txt
	@echo "$(LMAX) ! model harmonics max.degree" >> Input.txt
	@echo "0 ! FMM harmonics max.degree" >> Input.txt
	@echo "$(NGRID) ! number of Lebedev grid points" >> Input.txt
	@echo "14 ! order of convergence tolerance" >> Input.txt
	@echo "0 ! whether or not to compute forces" >> Input.txt
	@echo "78.3553 ! dielectric permittivity" >> Input.txt
	@echo "0.1 ! regularization parameter" >> Input.txt
	@cat Input.txt > $(OUTPUT).tmp # save current parameters
	@echo "$(shell head -n 1 $(INPUT)) ! number of atoms" >> $(OUTPUT).tmp
	@cat $(INPUT) >> Input.txt # copy grid into input
	@# Run dense solver
	@$(EXE) >> $(OUTPUT).tmp
	@# Copy output
	@mv $(OUTPUT).tmp $(OUTPUT)

clean:
	rm -rf $(OUTPUTS)

