# global flags

BUILD_YAML:=1

.PHONY: clean all realclean test bench
.SECONDARY:

# compiler flags

# g++:
CFLAGS = -g -O0
CXXFLAGS = -std=c++11 -g -O0 -fopenmp
CXXFLAGS += -MMD -MP  # auto-generate dependencies

# pgc++:
#CXXFLAGS = -std=c++11 -g -O3 -mp
#CXXFLAGS += -MMD # auto-generate dependencies


# main targets

demos := demo_parse
tests := test_dummy
benches := 

all: $(demos) $(tests) $(benches)

# directories

top=..
GTEST_DIR=$(top)/gtest-1.7.0
YAML_DIR=$(top)/yaml-0.1.6
VPATH = $(top)

CPPFLAGS += -I$(top)

# need to make libyaml ourselves?
ifdef BUILD_YAML
CPPFLAGS += -isystem $(YAML_DIR)/include
LDFLAGS += libyaml.a
yaml_lib = libyaml.a

obj/yaml/%.o: $(YAML_DIR)/src/%.c
	$(CC) -c -o $@ $(CPPFLAGS) -DHAVE_CONFIG_H -I$(YAML_DIR)/src $(CFLAGS) $<

libyaml.a: $(patsubst %.c,obj/yaml/%.o,$(notdir $(wildcard $(YAML_DIR)/src/*.c)))
	$(AR) $(ARFLAGS) $@ $^

else
yaml_lib =
endif

# rdmini core routines

rdmini_sources := $(notdir $(wildcard $(top)/rdmini/*.cc))
rdmini_objects := $(patsubst %.cc,obj/%.o,$(rdmini_sources))
-include $(patsubst %.o,%.d,$(rdmini_objects))

obj/%.o: rdmini/%.cc
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<


# demo apps link with rdmini object files

demo_%: obj/demo_%.o $(rdmini_objects) $(yaml_lib)
	$(CXX) -o $@ $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS)

obj/demo_%.o: demo/demo_%.cc
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<


# tests link with rdmini object files and libgtestmain

GTEST_CPPFLAGS=-isystem $(GTEST_DIR)/include

test_%: obj/test_%.o $(rdmini_objects) $(yaml_lib) libgtestmain.a
	$(CXX) -o $@ $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS)

obj/test_%.o: test/test_%.cc
	$(CXX) -c -o $@ $(CPPFLAGS) $(GTEST_CPPFLAGS) $(CXXFLAGS) $<


# build libgtestmain.a

GTEST_HEADERS=$(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SOURCES=$(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

GTEST_CPPFLAGS=-I$(GTEST_DIR) -I$(GTEST_DIR)/include -DGTEST_HAS_PTHREAD=0

obj/gtest-all.o: $(GTEST_SOURCES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GTEST_CPPFLAGS) -o $@ -c $(GTEST_DIR)/src/gtest-all.cc

obj/gtest_main.o: $(GTEST_SOURCES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GTEST_CPPFLAGS) -o $@ -c $(GTEST_DIR)/src/gtest_main.cc

libgtestmain.a: obj/gtest-all.o obj/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^


# clean

clean:
	rm -f gtest-all.o gtest_main.o obj/*.o obj/yaml/*.o

realclean: clean
	rm -f $(demos) $(tests) $(benches) libgtestmain.a libyaml.a obj/*.d obj/yaml/*.d


