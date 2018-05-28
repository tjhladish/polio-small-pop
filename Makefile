polio: main_Gillespie.cpp Makefile
	g++ -O2 --std=c++11 -Wall --pedantic main_Gillespie.cpp -o polio -lgsl -lgslcblas

check_polio: polio original_checksums
	rm output/*.csv
	time ./polio
	md5sum -c original_checksums
