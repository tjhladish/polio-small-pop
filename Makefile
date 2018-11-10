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
CPP = g++ -O2 --std=c++11 -Wall --pedantic

JSONINC := -Ijson/include

Params.o: Params.cpp Params.h json/include/nlohmann/json.hpp
	$(CPP) -c $< -o $@ $(JSONINC)

States.o: States.cpp States.h
	$(CPP) -c $< -o $@

libsim.a: Params.o States.o
	$(ARCHIVE) $@ $^

testODE: main_testEqODESoln.cpp libsim.a
	g++ -L. -O2 --std=c++11 -Wall --pedantic $< -o $@ -lgsl -lgslcblas -lsim

multiPatch: main_Gillespie_multivillage.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@ -lgsl -lgslcblas

risk: calculate_absolute_risk.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@

clean: *.o *.a
	rm $^