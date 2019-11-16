/*---------------------------------------------------------------
 *
 *   CLAME:"mapping"
 *
 *   
 *  
 ---------------------------------------------------------------*/
#include "mapping_lite.h"



int main(int argc,char *argv[])
{

    #ifdef debug
        cout<<"CLAME version_2.1 25/09/2017"<<endl;
        cout<<"Debug enable"<<endl;
    #endif
    //struct Args {bool multiFasta; bool fastq, bool outputFile; bool numT; bool bases_Threshold; bool print; bool fm9; bool offsetFM9; bool sizeBin;bool w;bool list2Exclude; bool size;};
    Args   args = {false,false,false,false,false,false,false,false,false,false,false,false};
    Names  names;
    //struct Parameters {int offsetFM9; int numThreads; int query_size; bool enablePrint; bool loadFM9;int sizeBin;bool w;int size};
    Parameters parameters={0,1,70,false,false,1000,0,0};
    names.outputFile="output";

    bool initOK=readArguments(argc,argv,&args,&names,&parameters);
    if (initOK)
    {
        parameters.loadFM9=args.fm9;
        parameters.enablePrint=args.print;
        bool runningError=false;
        vector<string> bases;         //Array to hold the bases of the reads
        vector<uint64_t> index;  
        if(args.size==false) //for the all dataset
        {
            //reading the multiFasta file into vectors 
            if(args.fastq)
                runningError=readFastQFile(&names,&bases);
            else
                runningError=readFasta(&names,&bases);
        
            //read file with list to exclude
            int numberOFreads=bases.size();
            int *reads2exclude = new int[numberOFreads] ();
            if(args.list2Exclude)
                runningError=readList2Exclude(&names,reads2exclude);

        
            //aligment
            vector <int> *MatrixList= new vector <int> [numberOFreads];
            runningError=alignemnt(&names, &parameters, &bases,  &index, MatrixList, reads2exclude);
        }
        else
        {
            csa_wt<> fm_index;
            csa_wt<> *ptrfm_index= &fm_index;
            loadFM9(&names, &parameters, &index, ptrfm_index);
            
            int numberOFreads=parameters.size;
            int *reads2exclude = new int[numberOFreads] ();
            readList2Exclude(&names,reads2exclude);
            vector <int> *MatrixList= new vector <int> [numberOFreads];
            if(args.fastq)
                runningError=readFastQFile2(&names,&bases,&parameters,MatrixList,reads2exclude,ptrfm_index,&index);
            else
                runningError=readFasta2(&names,&bases,&parameters,MatrixList,reads2exclude,ptrfm_index,&index);
        }
    }
    else
        printerror(argv[0]);

    return 1;
} //fin del main


