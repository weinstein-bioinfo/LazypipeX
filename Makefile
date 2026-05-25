
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
