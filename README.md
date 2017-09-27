/*---------------------------------------------------------------
//    ______  __          ___      .___  ___.  _______ 
//   /      ||  |        /   \     |   \/   | |   ____|
//  |  ,----'|  |       /  ^  \    |  \  /  | |  |__   
//  |  |     |  |      /  /_\  \   |  |\/|  | |   __|  
//  |  `----.|  `----./  _____  \  |  |  |  | |  |____ 
//   \______||_______/__/     \__\ |__|  |__| |_______|
//                                                     
/*---------------------------------------------------------------
 *
 *   CLAME:"Clasificador MetagenÃ³mico"
 *	
 *   CLAME is a binning software for metagenomic reads.
 *   It immplements a fm-index search algorithm for nucleotide 
 *   sequence alignment. Then it uses strongly connected component strategy
 *   to bin sequences with similar DNA composition.
 *  
 ---------------------------------------------------------------*/

---------------------------------------------------------------
Installation
---------------------------------------------------------------
1. To download and install the SDSL - Succinct Data Structure Library, use the instructions describe in the wep page.
https://github.com/simongog/sdsl-lite.

2. To compile CLAME, type make

3. Type clame -h to show the help

---------------------------------------------------------------
CLAME Versions
---------------------------------------------------------------
1. Type make little: Compile a lite version of CLAME, this version reduces the memory requirement to use reads' index instead of names.

2. Type make or make full: To compile a full version of CLAME, this version produces the final bin using the actuall name of the reads.


---------------------------------------------------------------
Authors
---------------------------------------------------------------
Benavides A(1), Alzate JF (2),(4) and Cabarcas F (1),(4)
1.	Grupo Sistemic, Departamento de Ingeniería Electrónica, Facultad de Ingenieria, Universidad de Antioquia.
2.	Centro Nacional de Secuenciacion Genomica-CNSG, Sede de Investigación Universitaria SIU, Universidad de Antioquia
3.	Escuela de Microbiología, Universidad de Antioquia
4.	Grupo de Parasitología, Departamento de Microbiología y Parasitología, Facultad de Medicina, Universidad de Antioquia

