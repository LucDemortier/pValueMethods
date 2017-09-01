# Specify the source files, target files, the build directories,
# and the install directory
OUTPUTS = poissonPvalues gaussianPvalues pValueCombination
LIBCDF = libcdf.a
CDFDIR = cdflib

# Compilation flags
CFLAGS = -Wall

# Add Gnu Scientific Library and Math Library to link string
LDFLAGS = -lgsl -lm

# Default target
.PHONY: all
all: $(OUTPUTS)
$(OUTPUTS):
	$(CXX) -o $@ $@.cpp $(CDFDIR)/$(LIBCDF) $(LDFLAGS)

# Clean up directory
.PHONY: clean
clean:
	for file in $(OUTPUTS); do rm -f $$file; done
