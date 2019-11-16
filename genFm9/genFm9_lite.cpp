/*---------------------------------------------------------------
 *
 *   CLAME:"Clasificador Metagenomico"
 *
 *   CLAME is a binning software for metagenomic reads.
 *   It immplements a fm-index search algorithm for nucleotide 
 *   sequence alignment. Then it uses strongly connected component strategy
 *   to bin the similar sequences.
 *   
 *  
 ---------------------------------------------------------------*/

#include "genFm9_lite.h"

bool readArguments(int argc,char *argv[], Args *args,Names *names, Parameters *parameters)
{
    bool   error = false, found1=false, mandatory=false;
    int    argum = 1;
    while (argum < argc && not error) 
    {
        found1 = false; 
        if ( IsParam(argv[argum],"-multiFasta") ) //multifasta file
        {
            argum ++;
            if (argum < argc)  
            {
                args->multiFasta = true;
                names->multiFasta = argv[argum];
                mandatory=true;
            } 
            else 
                error = true;
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-output") ) //output file
        {
            argum ++;
            if (argum < argc)  
            {
                args->outputFile = true;
                names->outputFile = argv[argum];
            } 
            else 
                error = true;
        
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-size") ) // number of reads to generate the FM9
        {
            argum ++;
            args->size = true;
            found1 = true;
            if (argum < argc)  
                parameters->size = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-fastq") ) 
        {
            argum ++;
            args->fastq = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-h") ) //to print the help
        {
            error = true;
            found1 = false; 
            continue;
        }
        if (not found1)
            error = true;
    }//end of the while

    if (error && not found1 ||  argc < 2  || not mandatory)
        return false;
    else
        return true;
}

bool readFasta2(Names *names, string *fasta, Parameters *parameters)
{
    //to print index
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".index";
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile.open(nameFile.c_str());
    //myfile2.open(nameFile2.c_str());

    ifstream infile;
    infile.open(names->multiFasta.c_str());
        
    if (infile.is_open()) 
    {
        string word, name="", read="";
        string auxNameFm9=names->outputFile;
        bool firstTime=true;
        double indexBases=0;
        int indexTitle=0;
        uint32_t n=0;
        int counReads=0;
        int bin=0;
        while (getline(infile,word))
        {
            if (word.find(">") != string::npos) 
            {
                if (firstTime)
                    firstTime=false;
                else
                {
                    counReads++;
                    myfile<<name<<"\t"<<indexTitle++<<endl; 
                    read+='$';
                    n = read.length();
                    *fasta+=read;
                    indexBases+=n;
                    read="";
                }
                int found = word.find(' ');
                name=word.substr (0,found);
            }
            else
                read=read+word;
            
            if(counReads==parameters->size)         //generate the FM9
            {
                counReads=0;
                names->outputFile=auxNameFm9+std::to_string(bin);
                bin++;
                alignemnt(names, parameters,fasta);
                *fasta="";

            }
            
        }
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        read+='$';
        n = read.length();
        *fasta+=read;
        indexBases+=n;
        read="";
        infile.close();
        myfile.close();
        
        //generate the FM9
        names->outputFile=auxNameFm9+std::to_string(bin);
        bin++;
        alignemnt(names, parameters,fasta);
        *fasta="";
                
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fasta file";
        #endif
        return true;
    }
}

bool readFastQFile2(Names *names, string* fasta, Parameters *parameters)
{
    //to print the index mandatory file
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".index";//ss.str();
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile.open(nameFile.c_str());

    int indexTitle=0;
    double indexBases=0;
    uint32_t n=0;

    //name file string word
    ifstream infile;
    infile.open(names->multiFasta.c_str());
    typedef enum {s0, s1,s2} STATES; //to read multifasta files
    STATES state=s0;

    if (infile.is_open()) 
    {
        string word, name="", read="", quality="";
        string auxNameFm9=names->outputFile;
        bool firstTime=true;
        int countLine=0;
        int counReads=0;
        int bin=0;
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                    {    
                        int found=word.find(' '); 
                        name=word.substr(0,found);
                        countLine=0; state=s1; 
                        myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                        counReads++;
                    }
                    break;
                case s1:
                    if (word[0]=='+') 
                        state=s2;
                    else
                        read=read+word, countLine++;
                    break;
                case s2:
                    //quality=quality+word;
                    if(--countLine==0)
                    {
                        read+='$';
                        n = read.length();
                        indexBases+=n;
                        *fasta+=read;
                        if(counReads==parameters->size)         //generate the FM9
                        {
                            counReads=0;
                            names->outputFile=auxNameFm9+std::to_string(bin);
                            bin++;
                            alignemnt(names, parameters,fasta);
                            *fasta="";
                            
                        }
                        read="";
                        state=s0;
                    }
                    break;
            }
        }
        
        if(counReads)         //generate the FM9 if there are reads without index
        {
            counReads=0;
            names->outputFile=auxNameFm9+std::to_string(bin);
            bin++;
            alignemnt(names, parameters,fasta);
            *fasta="";

        }  
        infile.close();
        myfile.close();
        //myfile2.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fastq file";
        #endif
        return true;
    }
}

bool readFasta(Names *names, string *fasta)
{
    //to print index
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".index";
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile.open(nameFile.c_str());
    //myfile2.open(nameFile2.c_str());

    ifstream infile;
    infile.open(names->multiFasta.c_str());
    
    if (infile.is_open()) 
    {
        string word, name="", read="";
        bool firstTime=true;
        double indexBases=0;
        int indexTitle=0;
        uint32_t n=0;
        //myfile2<<indexBases<<endl; 
        //index->push_back(indexBases);

        while (getline(infile,word))
        {
            if (word.find(">") != string::npos) 
            {
                if (firstTime)
                    firstTime=false;
                else
                {
                    myfile<<name<<"\t"<<indexTitle++<<endl; 
                    //title->push_back (name);
                    //bases->push_back(read);
                    read+='$';
                    n = read.length();
                    *fasta+=read;
                    indexBases+=n;
                    //myfile2<<indexBases-1<<endl; 
                    //index->push_back(indexBases);
                    read="";
                }
                int found = word.find(' ');
                name=word.substr (0,found);
            }
            else
                read=read+word;
        }
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        //title->push_back (name);
        //bases->push_back(read);
        read+='$';
        n = read.length();
        *fasta+=read;
        indexBases+=n;
        //index->push_back(indexBases);
        //myfile2<<indexBases-1<<endl; 
        read="";
        infile.close();
        
        myfile.close();
        //myfile2.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fasta file";
        #endif
        return true;
    }
}

bool readFastQFile(Names *names, string* fasta)
{
    //to print the index mandatory file
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".index";//ss.str();
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile.open(nameFile.c_str());
    //myfile2.open(nameFile2.c_str());

    int indexTitle=0;
    double indexBases=0;
    uint32_t n=0;

    //name file string word
    ifstream infile;
    infile.open(names->multiFasta.c_str());
    typedef enum {s0, s1,s2} STATES; //to read multifasta files
    STATES state=s0;

    if (infile.is_open()) 
    {
        string word, name="", read="", quality="";
        bool firstTime=true;
        int countLine=0;
        //myfile2<<indexBases<<endl; 
        //index->push_back(indexBases);
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                    {    
                        int found=word.find(' '); 
                        name=word.substr(0,found);
                        countLine=0; state=s1; 
                        myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                        //title->push_back (name);
                    }
                    break;
                case s1:
                    if (word[0]=='+') 
                        state=s2;
                    else
                        read=read+word, countLine++;
                    break;
                case s2:
                    //quality=quality+word;
                    if(--countLine==0)
                    {
                        //bases->push_back(read);
                        read+='$';
                        n = read.length();
                        indexBases+=n;
                        //index->push_back(indexBases);
                       //myfile2<<indexBases-1<<endl; 
                        *fasta+=read;
                        read="";
                        quality="";
                        state=s0;
                    }
                    break;
            }
        }
        infile.close();
        myfile.close();
        //myfile2.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "stage1: Error opening fastq file";
        #endif
        return true;
    }
}

bool alignemnt(Names *names,Parameters *parameters,string *fasta)
{
       
    //FM-Index
    csa_wt<> fm_index;
    string index_file = names->outputFile+".fm9";
    construct_im(fm_index,*fasta,1); 
    store_to_file(fm_index, index_file); // save it

}

int mstrlen(const char arg[])
{
  int out = 0;
  while(arg[out] != 0)
    out++;
  return out;
}

bool IsParam(char arg[],const char comp[])
{
  int len_comp = mstrlen(comp);
  bool equal = true;
  if (mstrlen(arg) == len_comp) 
  {
    for (int i = 0; i < len_comp; i++) 
      if (arg[i] != comp[i]) 
        equal = false;
  } 
  else 
    equal = false;
  
  return equal;
}

void printerror(const char arg[])
{
    cout << "genFM9 program version 2019-10-04" << endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h\t\t\t(Help)" << endl;
    cout << "  -fastq input file is in a fastq format  " << endl;
    cout << "  -multiFasta\t\tFILE  with all the reads " << endl;
    cout << "  -output name for the output-file  if print option was selected (default output)" << endl;
    cout << "  -size number of reads to generate the FM9 (default all dataset)" << endl;
    cout << ""<< endl;

}

