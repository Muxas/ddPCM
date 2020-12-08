# FMM solver
EXE=../main_fmm.exe
# Inputs geometries
INPUTS=$(wildcard *.inp)
# Possible values of lmax (model degrees) to check
LMAXS:=6 8 10 12
# Possible values of pmax (FMM degrees) to check
PMAXS:=0 2 4 6 8 10 12 15 18 21 24 27 30 35 40
# Size of Lebedev grids
NGRID:=590
# Outputs
OUTPUTS:=$(foreach lmax, $(LMAXS), $(foreach input, $(INPUTS:.inp=), \
	   $(foreach pmax, $(PMAXS), $(input).$(lmax).$(pmax).out)))

# Standard trick to ignore file "all" to build target "all"
.PHONY: all clean

# All jobs
all:	$(OUTPUTS)

# Set parameters of each FMM test
%.out:	ARGV=$(subst ., ,$*)
%.out:	INPUT=$(firstword $(ARGV)).inp
%.out:	LMAX=$(word 2, $(ARGV))
%.out:	PMAX=$(word 3, $(ARGV))
%.out:	OUTPUT=$@

# Recipe to prepare FMM result
%.out:
	@echo Running FMM experiment with input $(INPUT) lmax=$(LMAX) \
	    ngrid=$(NGRID) pmax=$(PMAX)
	@rm -f Input.txt
	@echo "0 ! print flag" > Input.txt
	@echo "1 ! number of OpenMP threads" >> Input.txt
	@echo "$(LMAX) ! model harmonics max.degree" >> Input.txt
	@echo "$(PMAX) ! FMM harmonics max.degree" >> Input.txt
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

