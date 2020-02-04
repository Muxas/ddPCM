EXE=./main.exe
EXE_FMM=./main_fmm.exe
INPUTS=$(wildcard *.inp)

OUTPUTS=$(patsubst %.inp, %.out, $(INPUTS))
OUTPUTS_FMM=$(patsubst %.inp, %.fmm, $(INPUTS))

.PHONY: all fmm dense clean

all: fmm dense

dense: $(OUTPUTS)

%.txt: %.inp
	@cat Input.txt.hdr $*.inp > $*.txt

%.out: %.txt
	@echo DENSE: $*
	@cp $*.txt Input.txt
	@$(EXE) > $*.out

fmm: $(OUTPUTS_FMM)

%.fmm: %.txt
	@echo FMM: $*
	@cp $*.txt Input.txt
	@$(EXE_FMM) > $*.fmm

clean:
	rm -f $(OUTPUTS) $(OUTPUTS_FMM) Input.txt
