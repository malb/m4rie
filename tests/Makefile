CXXFLAGS=-O2 -ggdb -I$(SAGE_ROOT)/local/include  -I.. 
LDFLAGS=-lm4rie -lm4ri -lgivaro -lntl -lgmpxx -lgmp -lm -lstdc++ -L/$(SAGE_ROOT)/local/lib

PRGS = test_elimination test_multiplication test_smallops

all:	$(PRGS)

.cc:	
	$(CXX) $(CXXFLAGS) $@.cc $(LDFLAGS) -o $@

clean:
	rm -f $(PRGS)