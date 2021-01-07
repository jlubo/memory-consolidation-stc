#-------------------------------------------------------------------------
SOURCE=NetworkMain.cpp
PROGRAM=net
COMPILER=g++
CFLAGS=-std=c++11
CINCLUDE=-Wl,--format=binary -Wl,code.zip -Wl,plotFunctions.py -Wl,--format=default
LDFLAGS=-static
LDLIBS=-lboost_serialization
#-------------------------------------------------------------------------
OBJECTS=$(SOURCE:.cpp=.o)
.PHONY: all test clean
all: $(SOURCE)
	rm -f code.zip
	zip -q -D code.zip * -x *.out *.o
	$(COMPILER) $(CFLAGS) $(CINCLUDE) $< $(LDFLAGS) $(LDLIBS) -o $(PROGRAM).out
	rm -f $(PROGRAM) *.o
	rm code.zip
clean:
	rm -f $(PROGRAM) *.o
	rm code.zip
