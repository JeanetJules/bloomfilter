#include <cstdint>
#include <iostream>


uint64_t xorshift64(uint64_t x);
void multihash(uint64_t x, uint64_t * hashes, uint64_t nb_hashes, uint64_t max_val);
