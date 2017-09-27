/*---------------------------------------------------------------
 *
 *   CLAME:"Clasificador Metagenómico"
 *	
 *   CLAME is a binning software for metagenomic reads.
 *   It immplements a fm-index search algorithm for nucleotide 
 *   sequence alignment. Then it uses strongly connected component strategy
 *   to bin sequences with similar DNA composition.
 *  
 ---------------------------------------------------------------*/

Installation

1. To download and install the SDSL - Succinct Data Structure Library, use the instructions describe in the wep page.
https://github.com/simongog/sdsl-lite.

2. To compile CLAME, type make

3. Type clame -h to show the help


CLAME Versions

1. Type make little: Compile a lite version of CLAME, this version reduces the memory requirement to use reads' index instead of names.

2. Type make or make full: To compile a full version of CLAME, this version produces the final bin using the actuall name of the reads.
