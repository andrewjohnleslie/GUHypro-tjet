
# Makefile for all

include		../Makefile.mk

all: html/index.html

clean:
	rm -f -R html

html/index.html: Doxyfile $(addprefix ../, $(SRC))
	@echo 'Building target documentation'
	@echo 'Invoking: doxygen'
	doxygen $<
	@echo 'Finished building documentation'
	@echo ' '


