DEBUG ?= 0

# External libraries
BOOST_ROOT ?= /usr
LIBBAM_ROOT ?= /usr
BOOST_SUFFIX ?=
LDADD ?=

CXX = g++
CXX_FLAGS += -isystem ${BOOST_ROOT}/include -isystem ${LIBBAM_ROOT}/include -std=c++11 -Wall -Wextra -fPIC -Wstrict-aliasing
LD_FLAGS += -L${BOOST_ROOT}/lib -lboost_iostreams${BOOST_SUFFIX} -lboost_program_options${BOOST_SUFFIX} -lboost_system${BOOST_SUFFIX} -lboost_thread${BOOST_SUFFIX} -L${LIBBAM_ROOT}/lib -lbam -lz -lpthread -lm ${LDADD}

ifeq (${DEBUG}, 1)
    CXX_FLAGS += -g -O0 -fno-inline
else
    CXX_FLAGS += -O3 -DNDEBUG
endif

TARGETS = src/rnasequel

SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:.cpp=.o)

all: $(TARGETS)

src/rnasequel: $(OBJS)
	$(CXX) -o $@ $^ $(LD_FLAGS) 

%.o : %.cpp 
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

clean:
	rm -f $(TARGETS) $(TARGETS:=.o) $(OBJS)
