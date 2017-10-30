/*---------------------------------------------------------------
 *
 *  CLAME:"Clasificador Metagenomico"
 *	
 *  CLAME is a binning software for metagenomic reads.
 *  It immplements a fm-index search algorithm for nucleotide 
 *  sequence alignment. Then it uses strongly connected component strategy
 *  to bin the similar sequences.
 *  
 *  Copyright (C) 2017 Benavides A.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 *   
 *  
 ---------------------------------------------------------------*/

#include "clame_lite.h"



int main(int argc,char *argv[])
{

    #ifdef debug
        cout<<"CLAME version_2.1 25/09/2017"<<endl;
        cout<<"Debug enable"<<endl;
    #endif
    //struct Args {bool multiFasta; bool fastq, bool outputFile; bool numT; bool bases_Threshold; bool print; bool fm9; bool lu; bool ld; bool sizeBin;};
    Args   args = {false,false,false,false,false,false,false,false,false,false};
    Names  names;
    Parameters parameters={0,10000,1,20,false,false,1000,false};
    names.outputFile="output";
    

    bool initOK=readArguments(argc,argv,&args,&names,&parameters);
    if (initOK)
    {
        parameters.fastq=args.fastq;
        parameters.loadFM9=args.fm9;
        bool runningError=false;
        
        vector<string> title;         //Array to hold the name of the reads
        vector<string> bases;         //Array to hold the bases of the reads
        vector<string> qual;          //Array to hold the bases of the reads' quality
        vector<uint32_t> index;       //Array to hold the bases->index assencion
        string fasta="";              //String to generate the FM-reference

        //reading the multiFasta file into vectors 
        if(parameters.fastq)
            runningError=readFastQFile(&names,&bases,&fasta,&index,&title,&qual);
        else
            runningError=readFasta(&names,&bases,&fasta,&index,&title);
        
        //aligment
        int numberOFreads=bases.size();
        int *queryList = new int[numberOFreads] ();
        vector <int> *MatrixList= new vector <int> [numberOFreads];
        runningError=alignemnt(&names, &parameters, &bases, &fasta, &index, queryList, MatrixList);
        
        //binning
        binningPrint(&names, &parameters, &title, queryList, MatrixList, numberOFreads, &bases, &qual);

    }
    else
        printerror(argv[0]);

    return 1;
} //end main




