poliomv: main_Gillespie_multivillage.cpp Makefile
	g++ -O2 --std=c++11 -Wall --pedantic main_Gillespie_multivillage.cpp -o poliomv -lgsl -lgslcblas

polio: main_Gillespie.cpp Makefile
	g++ -O2 --std=c++11 -Wall --pedantic main_Gillespie.cpp -o polio -lgsl -lgslcblas

check_polio: polio original_checksums
	rm output/*.csv
	time -p ./polio
	md5sum -c original_checksums

# git submodule add -b master --depth 1 https://github.com/open-source-parsers/jsoncpp
json/include/nlohmann/json.hpp:
	git submodule update --init --recursive --remote

# vs
# git submodule add -b master --depth 1 https://github.com/nlohmann/json

-include local.mk

ARCHIVE ?= $(AR) -rv
CPP := g++ -O2 --std=c++11 -Wall --pedantic

JSONINC := -Ijson/include

# library object for parsing model parameters
Params.o: Params.cpp Params.h json/include/nlohmann/json.hpp
	$(CPP) -c $< -o $@ $(JSONINC)

# defines the polio model compartments,
# provides definitions for transitions,
# provides assorted analytical definitions (e.g., equilibrium fractions)
States.o: States.cpp States.h
	$(CPP) -c $< -o $@

# implements basic Gillespie algorithm:
# given event rates, select which event occurs when 
Gillespie.o: Gillespie.cpp Gillespie.h
	$(CPP) -c $< -o $@

Utils.o: Utils.cpp Utils.h
	$(CPP) -c $< -o $@

PolioEvents.o: PolioEvents.cpp PolioEvents.h
	$(CPP) -c $< -o $@

libsim.a: Params.o States.o Gillespie.o Utils.o PolioEvents.o
	$(ARCHIVE) $@ $^

test%.o: test%.cpp libsim.a
	$(CPP) -L. $< -o $@ -lgsl -lgslcblas -lsim

polio_cabp: main_Gillespie_cabp.cpp libsim.a
	$(CPP) -L. $< -o $@ -lgsl -lgslcblas -lsim

multiPatch: main_Gillespie_multivillage.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@ -lgsl -lgslcblas

risk: calculate_absolute_risk.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@

clean: *.o *.a
	rm $^