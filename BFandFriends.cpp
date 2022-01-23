#include "BloomFilter.hpp"
#include "hash.hpp"


static uint64_t hashX( uint64_t x, uint64_t nf ) {    //passe x (un k-mer intified) nf fois dans xorshift, fournie
  if( x == 0 )
    x--;                                              //si x == 0, on le transforme en -1, or les bits sont non-signés

  for( uint64_t i = 0; i < nf; i++ )
    x = xorshift64( x );
  return x;
}


void BloomFilter::add_value( uint64_t x ) {     //ajoute à la case x la valeur x TODO : remplacer x par 1
  x = hashX( x, nf ) % size;
  array[ x/8 ] |= (1 << (x % 8 ));
}


bool BloomFilter::is_present( uint64_t x ) {    //renvoie true si la case du Filtre de Bloom est à autre chose que 0.
  x = hashX( x, nf ) % size;

  return array[ x/8 ] & (1 << (x % 8 ));       //notre tableau fait n/8 octets de taille, donc on va au x/8è octet, et on ajoute le nombre de bits restants pour avoir x
}                                              //TODO : Workaround avec le 0 !!!


BloomFilter::BloomFilter( size_t n, size_t nbf ) : size{ n }, nf{ nbf }{            //cette syntaxe fait que tu initialises les objets avant la construction et donc t'as pas besoin de les assigner !
  array = new uint8_t( n/8 + 1 );                                                    //donc les const peuvent être assignées à la création du BF, et du coup big
}


BloomFilter::~BloomFilter() {              //libère la mémoire lorsque l'on a plus besoin d'un objet BloomFilter
  delete[] array;
}


char skipline(FILE * fichier){            //passe à la ligne suivante, permet d'éviter la description
  char ch;
  while (ch != '\n') ch = fgetc(fichier); //si on arrive dans cette fonction, c'est qu'on a lu un '>', donc on attend juste le prochain '\n'
  if (ch == EOF) return EOF;              //en cas d'erreur
  return fgetc(fichier);                  //on renvoie l'élément d'après TODO : vérifier si on ne skip pas une valeur !!!
}


char nextTrueLetter(FILE * fichier) {      //passe à la prochaine lettre sous la forme ACTG, ignore les 'N'
  char ch = fgetc(fichier);
  while(ch == 'N' || ch =='\n') ch = fgetc(fichier);  //vérifie si la valeur d'après est aussi un 'N' (fréquent), mais aussi s'il y a un saut de ligne, ça causait un bug
  return ch;
}

uint64_t chartoInt ( char c ){
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'T': return 2;
    case 'G': return 3;
  }
}

uint64_t kmerToInt( char * kmer, size_t k) {          //transforme un k-mer en une valeur en base 4
  uint64_t intifiedkmer = 0;
  for( size_t i = 0; i < k; i++ ) {
    intifiedkmer += chartoInt(kmer[i]) << ( 2 * i );           //TODO : vérifier si c'est bien la base 4 (on a 2 * i ???)
  }
  return intifiedkmer;
}


char readNextLetter(FILE * fichier){               //envoie le prochain caractère d'ADN dans le fichier FASTA
  char ch;
  if ((ch = fgetc(fichier)) != EOF){                                      //On vérifie que le fichier n'est pas terminé
    if (ch == '>') { ch = skipline(fichier); }                          //On récupère la première lettre du génôme en passant au-delà de la description
    if (ch == '\n' || ch == '\r') { ch = fgetc(fichier); }              //On ignore les retours à la ligne tant qu'on est pas dans une description
    if (ch == 'N' || ch == 'n') { ch = nextTrueLetter(fichier); }       //On ignore les 'N'

    return ch;
  }
  return EOF;                       //en cas d'erreur
}


void reverseKmer( char * kmer, char * reversed, size_t k) {   //remplit reversed avec le Reverse Complement du k-mer donné en paramètre
  for( size_t i = 0; i < k; i++ ) {
    reversed[i] = kmer[ k - i - 1 ];            //ON LIT LE BUFFER A L'ENVERS c'est pas évident : AACT se lit aussi AAGT, donc on inverse la lecture de kmer pour avoir la lecture de reversed
    switch( reversed[i] ) {
      case 'A': reversed[i] = 'T'; break;
      case 'C': reversed[i] = 'G'; break;
      case 'G': reversed[i] = 'C'; break;
      case 'T': reversed[i] = 'A'; break;
    }
  }
}


char * pickWhichKmer( char * kmer, char * reversed, size_t k ){          //compare le kmer et choisit lequel des deux on prend entre le reverse complement et le kmer directement
  for( size_t i = 0; i < k; i++ ) {
    if( kmer[i] > reversed [i] ) return reversed;
    if( kmer[i] < reversed [i] ) return kmer;
  } return kmer;
}


bool nextkmer( FILE * fichier, char * kmer, char * reversed, size_t k ){     //renvoie le prochain kmer dans le fichier, ou bien son reverse complement
  char ch;

    for (size_t i = 1; i < k; i++) { kmer[ i - 1 ] = kmer[i]; }
    if ((ch = readNextLetter(fichier)) == EOF ) {
      printf("fichier termine \n\n\n" );
      return false;
    }
    kmer[ k - 1 ] = ch;
    reverseKmer( kmer, reversed, k );
    kmer = pickWhichKmer(kmer, reversed, k );   //on vérifie si on prend un k-mer ou son reverse complement selon l'ordre lexicographique
    return true;
  }


int main(int argc, char * argv[]){

  if( argc != 6 ) {
    std::cout << "Invalid number of arguments, expected 5 and got" << (argc - 1) << std::endl;
    exit(1);
  }

  char * pathname = argv[1];
  size_t k = atoi( argv[2] );
  size_t bloom_size = atoi( argv[3] );
  size_t nf = atoi( argv[4] );
  size_t requests = atoi( argv[5] );


  FILE *fichier;

  char *kmer = new char(k);                                           //pour économiser de la mémoire, ne pas oublier de le passer en argument
  char *reversed = new char(k);                                       //pour économiser de la mémoire, ne pas oublier de le passer en argument

  fichier = fopen(pathname, "r");                         //on obtient la taille du fichier
  if(fichier == NULL) printf("No such file\n" );
  fseek(fichier, 0L, SEEK_END);
  size_t filesize = ftell(fichier);
  fclose(fichier);

  fichier = fopen(pathname, "r");                  //on ouvre le fichier ici, une seule fois

  BloomFilter bloom( bloom_size, nf );
  for( size_t i = 0; i < k; i++ ){ kmer[i] = readNextLetter(fichier); }   //remplit le k-mer initial
  reverseKmer(kmer, reversed, k );
  bloom.add_value( kmerToInt( pickWhichKmer( kmer, reversed, k ), k ));   //ajoute le k-mer initial au Bloom-Filter
//ajoute la valeur hashée nf fois du k-mer ou de son reverse complément dans le bloom filter


  printf("On teste avant le while \n\n\n\n\n" );
  while (true) {
    if( !nextkmer(fichier, kmer, reversed, k) )
      break;
    reverseKmer( kmer, reversed, k );
    bloom.add_value( kmerToInt( pickWhichKmer( kmer, reversed, k ), k ));
  }
//ajoute la valeur hashée nf fois du k-mer ou de son reverse complement dans le bloom filter
  printf("On teste avant le fclose\n\n\n\n\n\n" );
  srand(time(NULL));
  const uint64_t test_max = (1 << (2 * k) - 1);
  uint64_t random_kmer;
  printf("on teste avant le for\n\n\n\n");
  for( size_t i = 0; i < requests; i++ ) {
    random_kmer = rand() % test_max;
    if( bloom.is_present( random_kmer ) ) printf("%lu \n", random_kmer);
  }


  return 0;
}
