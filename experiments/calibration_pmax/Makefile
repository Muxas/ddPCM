# Dense solver
EXE=./main.exe
# FMM solver
EXE_FMM=./main_fmm.exe
# Inputs geometries
INPUTS=$(wildcard *.inp)
# Possible values of lmax (model degrees) to check
LMAXS=$(shell seq 0 20)
# Possible values of pmax (FMM degrees) to check
PMAXS=$(shell seq 0 40)
# Sizes of Lebedev grids
NGRIDS=6 14 26 38 50 74 86 110 146 170 194 230 266 302 350 434 590 770 \
       974 1202 1454 1730 2030 2354 2702 3074 3470 3890 4334 4802 5294 5810
# Outputs for dense solver
DENSE_DIR=dense_out/
DENSE_OUTS=$(foreach lmax, $(LMAXS), $(foreach input, $(INPUTS:.inp=), \
	   $(input).$(lmax).out))
# Outputs for FMM solver
FMM_DIR=fmm_out/
FMM_OUTS=$(foreach lmax, $(LMAXS), $(foreach input, $(INPUTS:.inp=), \
	   $(foreach pmax, $(PMAXS), $(input).$(lmax).$(pmax).fmm)))

# Standard trick to ignore file "all" to build target "all"
.PHONY: all dense fmm clean

# Ask to build all outputs
all:	dense fmm

# Dense job
dense:	$(DENSE_DIR) $(DENSE_OUTS)

# Dense output dir
$(DENSE_DIR):
	mkdir $(DENSE_DIR)

# FMM job
fmm:	$(FMM_DIR) $(FMM_OUTS)

# FMM output dir
$(FMM_DIR):
	mkdir $(FMM_DIR)

# Set parameters of each dense test
%.out:	ARGV=$(subst ., ,$*)
%.out:	INPUT=$(firstword $(ARGV)).inp
%.out:	LMAX=$(lastword $(ARGV))
%.out:	OUTPUT=$(DENSE_DIR)/$(firstword $(ARGV)).$(LMAX).out

# Recipe to prepare dense result
%.out:
	@echo Running dense experiment with input $(INPUT) and lmax=$(LMAX)
	@rm -f Input.txt
	@echo 0 > Input.txt # print flag
	@echo 1 >> Input.txt # number of OMP threads
	@echo $(LMAX) >> Input.txt # degree of harmonics of the model
	@echo $(word $(LMAX), $(NGRIDS)) >> Input.txt # number of grid points
	@echo 14 >> Input.txt # convergence threshold
	@echo 0 >> Input.txt # whether or not to compute forces
	@echo 78.3553 >> Input.txt # dielectric constant
	@echo 0.2 >> Input.txt # regularization parameter
	@cat $(INPUT) >> Input.txt # grid
	@# Run dense solver
	@echo $(EXE) > $(OUTPUT).tmp
	@# Copy output
	@mv $(OUTPUT).tmp $(OUTPUT)

# Set parameters of each FMM test
%.fmm:	ARGV=$(subst ., ,$*)
%.fmm:	INPUT=$(firstword $(ARGV)).inp
%.fmm:	LMAX=$(word 2,$(ARGV))
%.fmm:	PMAX=$(lastword $(ARGV))
%.fmm:	OUTPUT=$(FMM_DIR)/$(firstword $(ARGV)).$(LMAX).$(PMAX).out

# Recipe to prepare dense result
%.fmm:
	@echo Running FMM experiment with input $(INPUT) lmax=$(LMAX) and \
	    pmax=$(PMAX)
	@rm -f Input.txt
	@echo 0 > Input.txt # print flag
	@echo 1 >> Input.txt # number of OMP threads
	@echo $(LMAX) >> Input.txt # degree of harmonics of the model
	@echo $(PMAX) >> Input.txt # degree of harmonics of the FMM
	@echo $(word $(LMAX), $(NGRIDS)) >> Input.txt # number of grid points
	@echo 14 >> Input.txt # convergence threshold
	@echo 0 >> Input.txt # whether or not to compute forces
	@echo 78.3553 >> Input.txt # dielectric constant
	@echo 0.2 >> Input.txt # regularization parameter
	@cat $(INPUT) >> Input.txt # grid
	@# Run dense solver
	@echo $(EXE) > $(OUTPUT).tmp
	@# Copy output
	@mv $(OUTPUT).tmp $(OUTPUT)

clean:
	rm -rf $(DENSE_DIR) $(FMM_DIR)


