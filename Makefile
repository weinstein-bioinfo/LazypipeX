
# TESTS
#
# Phony targets: they name no file, so they always run.  Without .PHONY, make
# would silently do nothing here the day a file called `test` appears.
#
# `make test` is the reflex entry point and deliberately runs only the tiers that
# need no reference databases, so it stays fast enough to type without thinking.
# `make test-all` adds Tiers 2-4, which need the databases and take ~25 min.
# tests/run_tests.sh is the real interface: --tier, --keep-going, --list.
.PHONY: test test-all

test:
	tests/run_tests.sh --quick

test-all:
	tests/run_tests.sh

# UTILITIES USING bioio.h

retrieve_reads: cpp/retrieve_reads.cpp
	g++ -Wall -O3 -std=c++11 -Icpp cpp/retrieve_reads.cpp -o bin/retrieve_reads

# UTILITIES BASED ON SEQAN LIBRARY

get_contigs: cpp/get_contigs.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/get_contigs.cpp -o bin/get_contigs
filtfa: cpp/filtfa.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/filtfa.cpp -o bin/filtfa
filtfq: cpp/filtfq.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/filtfq.cpp -o bin/filtfq

# STATIC BUILDS:
sretrieve_reads: cpp/retrieve_reads.cpp
	g++ -Wall -O3 -std=c++11 -Icpp cpp/retrieve_reads.cpp -o bin/retrieve_reads  -static
sget_contigs: cpp/get_contigs.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/get_contigs.cpp -o bin/get_contigs  -static
sfiltfa: cpp/filtfa.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/filtfa.cpp -o bin/filtfa -static
sfiltfq: cpp/filtfq.cpp
	g++ -Wall -O3 -DNDEBUG -std=c++14 -Icpp -I${seqan}/include cpp/filtfq.cpp -o bin/filtfq -static
	
#-lrt -lpthread

# options
#-std=c++14, -std=gnu++14
	
#clean:
#	rm *.o
