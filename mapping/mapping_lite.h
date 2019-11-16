/*---------------------------------------------------------------
 *
 *   CLAME:"mapping"
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

struct Args {bool multiFasta; bool fastq; bool outputFile; bool numT; bool bases_Threshold; bool print; bool fm9; bool offsetFM9; bool sizeBin; bool w; bool list2Exclude; bool size;};
struct Names {string multiFasta; string outputFile; string fm9; string list2Exclude;};
struct Parameters {int offsetFM9; int numThreads; int query_size; bool enablePrint; bool loadFM9;int sizeBin; int w; int size;};

string reverse(string *str);
int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);
bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters);

bool readFasta(Names *names, vector<string> *bases);
bool readFastQFile(Names *names, vector<string> *bases);


bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases, vector<uint64_t> *index, vector<int>* MatrixList, int *reads2exclude);
void printResult(Names *names, Parameters *parameters,int numberOFreads, vector<int> *MatrixList,int *reads2exclude);
uint32_t indx2Loc(uint64_t locs, vector<uint64_t> *index);
bool readList2Exclude(Names *names, int *reads2exclude);



bool readFasta2(Names *names, vector<string> *bases,Parameters *parameters, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index,vector<uint64_t> *index);
bool readFastQFile2(Names *names, vector<string> *bases,Parameters *parameters, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index,vector<uint64_t> *index);
bool loadFM9(Names *names, Parameters *parameters, vector<uint64_t> *index, csa_wt<> *ptrfm_index);
bool alignemnt2(Names *names,Parameters *parameters,vector<string> *bases, vector<uint64_t> *index, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index,int offset);
void printResult2(Names *names, Parameters *parameters,int numberOFreads, vector<int> *MatrixList,int *reads2exclude, int offset, ofstream *myfileResult, ofstream *myfileLinks);

