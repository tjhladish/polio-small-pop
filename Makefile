poliomv: main_Gillespie_multivillage.cpp Makefile
	g++ -O2 --std=c++11 -Wall --pedantic main_Gillespie_multivillage.cpp -o poliomv -lgsl -lgslcblas

polio: main_Gillespie.cpp Makefile
	g++ -O2 --std=c++11 -Wall --pedantic main_Gillespie.cpp -o polio -lgsl -lgslcblas

check_polio: polio original_checksums
	rm output/*.csv
	time -p ./polio
	md5sum -c original_checksums

testODE: main_testEqODESoln.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@ -lgsl -lgslcblas

multiPatch: main_Gillespie_multivillage.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@ -lgsl -lgslcblas

risk: calculate_absolute_risk.cpp
	g++ -O2 --std=c++11 -Wall --pedantic $< -o $@ 