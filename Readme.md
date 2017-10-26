---------------------------------------------------------------
    CLAME:"Clasificador Metagenomico"
---------------------------------------------------------------
Description:
    CLAME is a binning software for metagenomic reads.
    It immplements a fm-index search algorithm for nucleotide 
    sequence alignment. Then it uses strongly connected component strategy
    to bin sequences with similar DNA composition.

---------------------------------------------------------------
Installation
---------------------------------------------------------------
1. To download and install the SDSL - Succinct Data Structure Library, use the instructions describe in the wep page.
https://github.com/simongog/sdsl-lite.

2. Download CLAME source codes into the same directory used for the SDSL codes

3. To compile CLAME, type make

4. Type clame -h to show the help

---------------------------------------------------------------
CLAME Versions
---------------------------------------------------------------
1. Type make or make all: To Compile CLAME.

2. Type make or make full: To compile a debug version of CLAME.

---------------------------------------------------------------
Usage
---------------------------------------------------------------
./clame -b 70 -multiFasta test/Bancomini.fna -nt 4 -output bminiBins 

Output files
1. bminiBins_*.fastq: Output fastq
2. bminiBins.bins: All bins reported
3. bmini.fm9: FM-index output

---------------------------------------------------------------
Authors
---------------------------------------------------------------
Benavides A(1), Alzate JF (2),(3) and Cabarcas F (1),(3)
1.	Grupo Sistemic, Departamento de Ingeniería Electrónica, Facultad de Ingenieria, Universidad de Antioquia.
2.	Centro Nacional de Secuenciacion Genomica-CNSG, Sede de Investigación Universitaria SIU, Universidad de Antioquia
3.	Grupo de Parasitología, Departamento de Microbiología y Parasitología, Facultad de Medicina, Universidad de Antioquia

