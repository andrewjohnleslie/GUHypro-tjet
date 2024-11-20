
# Makefile for all

include		../Makefile.mk

TARGETS		= $(BinPath)/mixtureMex.mexa64 $(BinPath)/thermoKineticStateMex.mexa64 $(BinPath)/thermoStateMex.mexa64

all: $(TARGETS)

$(BinPath)/%.mexa64: %.cpp
	@echo 'Building target: $@'
	@echo 'Invoking: Matlab mex file builder'
	mex -v $(Include) -g CXXFLAGS='$(CXXFLAGS) $(CPPFLAGS)' $(Library) $<
	#@mv $(<:.cpp=.d) $(@:.mexa64=.d)
	@mv $(<:.cpp=.mexa64) $@
	#@sed 's@'$(<:.cpp=.o)'@'$@'@' $(@:.mexa64=.d) > $(@:.mexa64=.d.tmp) #fix problem with dependancy file
	#@rm $(@:.mexa64=.d)
	#@mv $(@:.mexa64=.d.tmp) $(@:.mexa64=.d)
	#@rm -f $(@:Mex.mexa64=.m)
	@ln -s $(ProjDirPath)/src/thermo/$(<:Mex.cpp=.m) $(BinPath)
	@echo 'Finished building target: $@'
	@echo ' '

#ifneq ($(strip $(TARGETS:.mexa64=.d)),)
#-include $(TARGETS:.mexa64=.d)
#endif


