CXX = g++
CXXFLAGS = -O3 -Wall
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

all: ngsAssociation

.PHONY: clean

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp > $*.d

ngsAssociation: $(OBJ)
	$(CXX) $(CXXFLAGS) -o ngsAssociation *.o

clean:
	rm -f ngsAssociation *.o *.d

