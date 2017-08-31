# Specify extensions of files to delete when cleaning
CLEANEXTS = o d
CLEANALLEXTS = o d a

# Compiler flags
CFLAGS = -O2 -Wall

# Specify the target file and the install directory
OUTPUTFILE = libcdf.a
INSTALLDIR = ./

# Default target
.PHONY: all
all: $(OUTPUTFILE)

# List of sources
SOURCES = $(wildcard *.cpp)

# Generate object files
%.o: %.cpp $(DEPS)
	$(CC) -MMD -c -o $@ $< $(CFLAGS)

# Include dependency make files, if any exist.
-include *.d

# Build library from object files
$(OUTPUTFILE): $(subst .cpp,.o,$(SOURCES))
	ar ru $@ $^
	ranlib $@

# No rule to build object files from source files is required;
# this is handled by make's database of implicit rules
.PHONY: install
install:
	mkdir -p $(INSTALLDIR)
	cp -p $(OUTPUTFILE) $(INSTALLDIR)

.PHONY: clean
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

.PHONY: cleanall
cleanall:
	for file in $(CLEANALLEXTS); do rm -f *.$$file; done
