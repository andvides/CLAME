/*---------------------------------------------------------------
 *
 *   CLAME:"Clasificador Metagenómico"
 *
 *   CLAME is a binning software for metagenomic reads.
 *   It immplements a fm-index search algorithm for nucleotide 
 *   sequence alignment. Then it uses strongly connected component strategy
 *   to bin the similar sequences.
 *   
 *  
 ---------------------------------------------------------------*/
#include <sdsl/suffix_arrays.hpp>
#include "iostream"
#include<fstream>
#include <ctime>

#include <fcntl.h>
#include <sstream> 
#include <vector>
#include <omp.h>
#include <string>
#include <algorithm>
#include <iomanip>

#include<cstdint>
#include <cstdlib>
#include <tr1/unordered_map>
#include <vector>

#include <sys/resource.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//using std::unordered_map;
using std::tr1::unordered_map;

using namespace sdsl;
using namespace std;

struct Args {bool multiFasta; bool fastq; bool outputFile; bool size;};
struct Names {string multiFasta; string outputFile; string fm9;};
struct Parameters {int size;};

int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);
bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters);
bool readFasta(Names *names,string *fasta);
bool readFastQFile(Names *names,string* fasta);
bool alignemnt(Names *names,Parameters *parameters,string *fasta);
bool readFasta2(Names *names,string *fasta,Parameters *parameters);
bool readFastQFile2(Names *names,string* fasta,Parameters *parameters);

