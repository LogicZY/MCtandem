#makefile for c++ programs
#change the name of the executable to use for "any" project

EXECUTABLE = ../bin/MCtandem.exe
#EXECUTABLE = ../bin/p3.exe
LINKCC = $(CXX)

# Use openMP
ifeq ($(openmp),1)
LIB += -lgomp
endif

#CXXFLAGS denotes flags for the C++ compiler

CXX = icpc

#uncomment this line if you are using gcc 4.x
CXXFLAGS = -O2 -DGCC4_3
#CXXFLAGS = -O2 -DGCC4_3 -DX_P3
#CXXFLAGS = -O2 -DGCC4_3 -DPLUGGABLE_SCORING -DX_P3

#LDFLAGS = -lpthread -L/usr/lib -lm /usr/local/lib/libexpat.a -static-libstdc++
#LDFLAGS = -lpthread -L/usr/lib -lm /usr/local/lib/libexpat.a
LDFLAGS = -lpthread -L/system/local/lib -lm -lexpat
SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(EXECUTABLE)

#define the components of the program, and how to link them
#these components are defined as dependencies; that is they must be up-to-date before the code is linked

$(EXECUTABLE): $(DEPS) $(OBJS)
	$(LINKCC) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)

#specify the dep files depend on the cpp files

%.d: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@
	$(CXX) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@

clean:
	-rm $(OBJS) $(EXECUTABLE) $(DEPS) *~

explain:
	@echo "The following info represents the program:"
	@echo "Final exec name: $(EXECUTABLE)"
	@echo "Source files:       $(SRCS)"
	@echo "Object files:       $(OBJS)"
	@echo "Dep files:          $(DEPS)"

depend: $(DEPS)
	@echo "Deps are now up-to-date."
 	
-include $(DEPS)
