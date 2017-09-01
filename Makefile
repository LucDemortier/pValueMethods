# Specify the target files and the libraries to link to.
OUTPUTS = poissonPvalues gaussianPvalues pValueCombination
LIBCDF = cdflib/libcdf.a
LDLIBS = -lgsl -lm

# Default target
.PHONY: all
all: $(OUTPUTS)
$(OUTPUTS): %: %.cpp
	$(CXX) -o $@ $< $(LIBCDF) $(LDLIBS)

# Clean up directory
.PHONY: clean
clean:
	for file in $(OUTPUTS); do rm -f $$file; done
