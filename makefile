# Include the MOAB configuration information so that
# all required flags and libs are populated correctly
include makefile.config

default: all

ALLEXAMPLES = malha_hexaedros

all: $(ALLEXAMPLES)

malha_hexaedros: malha_hexaedros.o
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}
