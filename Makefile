CPP=g++ --std=c++11 -Wall

all: BFandFriends

BFandFriends: hash.o BFandFriends.o
	$(CPP) -o $@ $^

hash.o: hash.cpp hash.hpp
	$(CPP) -c $<

BFandFriends.o: BFandFriends.cpp BloomFilter.hpp
	$(CPP) -c $<


clean :
	rm -rf BFandFriends *.o
