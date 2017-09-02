# Specify the target files and the libraries to link to.
OUTPUTS = poissonPvalues gaussianPvalues pValueCombination
CDFDIR = cdflib
LIBCDF = libcdf.a
LDLIBS = -lgsl -lm

# Default target
.PHONY: all
all: $(OUTPUTS)
$(OUTPUTS): %: %.cpp $(CDFDIR)/$(LIBCDF)
	$(CXX) -o $@ $< $(CDFDIR)/$(LIBCDF) $(LDLIBS)

# Target to create CDF library
.PHONY: libs
libs: $(CDFDIR)/$(LIBCDF)
$(CDFDIR)/$(LIBCDF):
	cd $(CDFDIR) && $(MAKE)
	cd $(CDFDIR) && $(MAKE) clean

# Clean up directory
.PHONY: clean
clean:
	for file in $(OUTPUTS); do rm -f $$file; done
