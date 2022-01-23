#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <time.h>


class BloomFilter {                       //initialisation de la classe.
public:                                   //en public (utilisable partout) : le constructeur & destructeur,
  BloomFilter( size_t s, size_t nf );     //les fonctions add_value() et is_present(), et la taille du Filtre.
  virtual ~BloomFilter();
  void add_value( uint64_t x );
  bool is_present( uint64_t x );
  const size_t size;
private:                                  //en privé, le tableau en lui-même, et le nombre de fonctions de hash
  uint8_t * array;
  const size_t nf;
};
