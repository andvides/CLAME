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
#include "genFm9_lite.h"



int main(int argc,char *argv[])
{

    #ifdef debug
        cout<<"CLAME version_2.1 25/09/2017"<<endl;
        cout<<"Debug enable"<<endl;
    #endif
    //struct Args {bool multiFasta; bool fastq; bool outputFile; bool size;};                                                                         
    //struct Names {string multiFasta; string outputFile; string fm9;};
    //struct Parameters {int size;};   
    Args   args = {false,false,false,false};
    Names  names;
    Parameters parameters={1};
    names.outputFile="output";

    bool initOK=readArguments(argc,argv,&args,&names,&parameters);
    if (initOK)
    {
        bool runningError=false;
        string fasta=""; //String to generate the FM-reference

        if(!args.size) //genarate the FM9 from all dataset
        {
            //cout<<"ALL reads";
            if(args.fastq)
                runningError=readFastQFile(&names,&fasta);
            else
                runningError=readFasta(&names,&fasta);
            runningError=alignemnt(&names, &parameters, &fasta);
        }
        else //split the file to genarate the FM9 from all dataset
        {
            //cout<<"By group of  reads";
             if(args.fastq)
                runningError=readFastQFile2(&names,&fasta,&parameters);
            else
                runningError=readFasta2(&names,&fasta,&parameters);
        }

    }
    else
        printerror(argv[0]);

    return 1;
} //fin del main

