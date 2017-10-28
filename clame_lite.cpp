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
        if ( IsParam(argv[argum],"-fm9") ) //FM) reference
        {
            argum ++;
            if (argum < argc)  
            {
                args->fm9 = true;
                names->fm9 = argv[argum];
            } 
            else 
                error = true;
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-lu") ) //superior cut
        {
            argum ++;
            if (argum < argc)  
            {
                args->lu = true;
            parameters->lu = atoi(argv[argum]);
            } 
            else 
                error = true;
         
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-ld") ) //inferior cut
        {
            argum ++;
            if (argum < argc)  
            {
                args->ld = true;
                parameters->ld = atoi(argv[argum]);
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
        if ( IsParam(argv[argum],"-nt") ) // number of cpus
        {
            argum ++;
            args->numT = true;
            found1 = true;
            if (argum < argc)  
                parameters->numThreads = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-b") ) //minimal base lenght
        {
            argum ++;
            args->bases_Threshold = true;
            found1 = true;
            if (argum < argc)  
                parameters->query_size = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-print") ) //print to file
        {
            argum ++;
            parameters->enablePrint = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-sizeBin") ) //minimal reads set to report a bin
        {
            argum ++;
            if (argum < argc)  
            {
                args->sizeBin = true;
                parameters->sizeBin = atoi(argv[argum]);
            } 
            else 
                error = true;
                       
            argum ++;
            found1 = true;
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

bool readFasta(Names *names, vector<string> *bases,string *fasta,vector<uint32_t> *index)
{
    
    //to print index
    //mandatory file
    ofstream myfile;
    string nameFile=names->outputFile+".index";
    myfile.open(nameFile.c_str());

    ifstream infile;
    infile.open(names->multiFasta.c_str());
    
    if (infile.is_open()) 
    {
        string word, name="", read="";
        bool firstTime=true;
        int indexBases=0;
        int indexTitle=0;
        uint32_t n=0;
        index->push_back(indexBases);

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
                    bases->push_back(read);
                    read+='$';
                    n = read.length();
                    *fasta+=read;
                    indexBases+=n;
                    index->push_back(indexBases);
                    read="";
                }
                name=word;
            }
            else
                read=read+word;
        }
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        //title->push_back (name);
        bases->push_back(read);
        read+='$';
        n = read.length();
        *fasta+=read;
        indexBases+=n;
        index->push_back(indexBases);
        read="";
        infile.close();
        
        myfile.close();
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

bool readFastQFile(Names *names,  vector<string>* bases, string* fasta, vector<uint32_t> *index)
{
    //to print the index mandatory file
    ofstream myfile;
    string nameFile=names->outputFile+".index";//ss.str();
    myfile.open(nameFile.c_str());

    int indexTitle=0;
    int indexBases=0;
    uint32_t n=0;

    //name file string word
    ifstream infile;
    infile.open(names->multiFasta.c_str());
    typedef enum {s0, s1,s2} STATES; //to read multifasta files
    STATES state=s0;

    if (infile.is_open()) 
    {
        string word, name="", read="";
        bool firstTime=true;
        int countLine=0;
        index->push_back(indexBases);
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                        name=word, countLine=0, state=s1, myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                    break;
                case s1:
                    if (word[0]=='+') 
                        name=word, state=s2;
                    else
                        read=read+word, countLine++;
                    break;
                case s2:
                    if(--countLine==0)
                    {
                        bases->push_back(read);
                        read+='$';
                        n = read.length();
                        indexBases+=n;
                        index->push_back(indexBases);
                        *fasta+=read;
                        read="";
                        state=s0;
                    }
                    break;
            }
        }
        infile.close();
        myfile.close();
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

bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases,string *fasta, vector<uint32_t> *index, int* queryList, vector<int>* MatrixList)
{
       
    std::vector<string>::iterator ptrBases= bases->begin();
   
    //FM-Index
    csa_wt<> fm_index;
    string index_file =""; 
    if(parameters->loadFM9)//load from file
    {

        index_file   = names->fm9;
        if (!load_from_file(fm_index, index_file)) 
        {
            #ifdef debug
                cout << "No index "<<index_file<< "error" << endl;
            #endif
            return true;
        }
        else
        {
            index->clear();
            index->push_back(0);
            string querySearch="$";
            auto locations = locate(fm_index, querySearch.begin(), querySearch.begin()+1);
            sort(locations.begin(), locations.end());
           
            for (auto it=locations.begin(); it<locations.end(); it++)
                index->push_back(*it);
        }
    }
    else//buidl suffix tree
    {
    
        index_file   = names->outputFile+".fm9";
        construct_im(fm_index,*fasta,1); 
        *fasta="";
        if(parameters->enablePrint)
            store_to_file(fm_index, index_file); // save it
    }
    
    //search
    int cpus=parameters->numThreads;
    int numberOFreads=bases->size();
    bool runningError=false;
    
    if(runningError==false)
    {
        #pragma omp parallel num_threads(cpus)
        {
            #pragma omp for schedule(runtime) 
            for(int i=0; i<numberOFreads;i++)//process all parts except the last one because it can get a diffent size
            {
                uint32_t subjectID=0;
                string query=*(ptrBases+i);
                if (query.length()>=parameters->query_size)
                {
                    string reverseQuery=reverse(query);
                    
                    string queryF1=query.substr(0,parameters->query_size);//forward init
                    string queryF2=query.substr((query.length()-parameters->query_size),parameters->query_size);;//forward end
                    string queryR1=reverseQuery.substr(0,parameters->query_size);//reverse init
                    string queryR2=reverseQuery.substr((reverseQuery.length()-parameters->query_size),parameters->query_size);//reverse end
    
                    string Queries[4]={queryF1,queryF2,queryR1,queryR2};
                    string orientations[4]={"5a'","5b'","3a'","3b'"};
	
                    queryList[i]=0;// Generate the nodes of the tree and init as unmarked read
                    for (int q=0;q<4;q++)
                    {
                        string querySearch=Queries[q];
                        uint32_t sizeq = querySearch.length();
                        auto locations = locate(fm_index, querySearch.begin(), querySearch.begin()+sizeq);
                        for (auto it=locations.begin(); it<locations.end(); it++)
                        {
                            subjectID=indx2Loc(*it,index);
                            #pragma omp critical
                            {
                                MatrixList[i].push_back(subjectID);
                                MatrixList[subjectID].push_back(i);
                            }
                        }
                    }
                }
               
           
            }
        }//end pragma
        
        //sort and deleted duplicated items to get the number of links
        for(int k=0;k<numberOFreads;k++)
        {
            std::sort (MatrixList[k].begin(), MatrixList[k].end());           
            std::vector<int>::iterator it=std::unique (MatrixList[k].begin(), MatrixList[k].end()); 
            MatrixList[k].resize( std::distance(MatrixList[k].begin(),it) );
        }
        
        //print the result
        if(parameters->enablePrint)
            printResult(names,numberOFreads,MatrixList);
            
        return runningError;
    }
}

void printResult(Names *names, int numberOFreads, vector<int> *MatrixList)
{
    ofstream myfile, myfile2;
    string nameFile=names->outputFile+".result";//ss.str();
    string nameFile2=names->outputFile+".links";//ss.str();
    myfile.open(nameFile.c_str());
    myfile2.open(nameFile2.c_str());

    for(int k=0;k<numberOFreads;k++)
    {
        //print links
        myfile2<<k<<"\t"<<MatrixList[k].size()<<endl;
        
        //print result
        myfile<<k<<"\t";
        for (auto it=MatrixList[k].begin(); it<MatrixList[k].end(); it++)
            myfile<<*it<<"\t";
        myfile<<"\n";
    }
    myfile.close();
    myfile2.close();
}

void binningPrint(Names *names, Parameters *parameters, vector<string> *title, int *queryList, vector<int> *MatrixList, int numberOFreads, vector<string> *bases)
{
    int cpus=parameters->numThreads;
    std::vector<string>::iterator ptrBases= bases->begin();
   
    //1.Number of links by read
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime)
        for(int k=0; k<numberOFreads;k++) 
            if(MatrixList[k].size()<parameters->ld || MatrixList[k].size()>parameters->lu )
                queryList[k]=1; //marked as visited 
    }
    
    
    //2. Genera los bins
    ofstream myfile;
    string nameFile=names->outputFile+".binning";
    myfile.open(nameFile.c_str());
    
                
    int *stack = new int[numberOFreads] ();
    int *put=stack;
    int *get=stack;
    int numBin=0; //0: no visietd 1:visited
    for(int i=0; i<numberOFreads;i++)
    {
        if(queryList[i]==0) //not visited
        {
            queryList[i]=1;
            *(put++)=i;
            while(get!=put)
            {
                int q=*get++;
                for( vector<int>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                {
                    if(queryList[*j]==0) //not visited
                    {
                        queryList[*j]=1;
                        *(put++)=*j;
                    }
                }
            }

            //PRINT BIN
            int size=get-stack;
            if(size>parameters->sizeBin)
            {
                float links1=0.0, mean1=0.0, std1=0.0;
                double *localLinks = new double[size]; 
                int p=0;
                
                ofstream myfile2;
                stringstream ss;
                ss << numBin;
                string strNumBin = ss.str();
                string nameFile2=names->outputFile+"_"+strNumBin+".fasta";
                myfile2.open(nameFile2.c_str());

                myfile<<">Bin:"<<numBin<<":\t"<<size<<endl;
                while(get!=stack)
                {
                    --get;
                    myfile<<*get<<"\t"<<numBin<<endl;
                    
                    string base=*(ptrBases+*(get));
                    myfile2<<">Read_"<<*get<<endl;
                    myfile2<<base<<endl;
                }
                numBin++;
                myfile2.close();

            }
            put=stack; //pointers restart
            get=stack; 

		}
	}
    myfile.close();
}

uint32_t indx2Loc(uint32_t locs, vector<uint32_t> *index)
{
	std::vector<uint32_t>::iterator up;
	up= std::upper_bound (index->begin(), index->end(), locs); 
	return (up-1-index->begin());
}

string reverse(string str)
{
// This table is used to transform nucleotide letters into reverse. //A=T, C=G, G=C, T=A, Other=N
  static const char  nt_table_reverse[128] = {
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
        };
  
        string rs="";
        for (std::string::reverse_iterator rit=str.rbegin(); rit!=str.rend(); ++rit)
        {
            char base=nt_table_reverse[*rit]; 
            rs=rs+base;
        }
        //cout<<str<<"\t"<<rs<<endl;
        return rs;

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
    cout << "CLAME:'Clasificador Metagenomico'" << endl;
    cout << "Bin name and reads into of this bin"<<endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h\t\t\t(Help)" << endl;
    cout << "  -b minimum number of bases to take an alignment (default 20) " << endl;
    cout << "  -fm9 Load fm9 file  " << endl;
    cout << "  -fastq input file is in a fastq format  " << endl;
    cout << "  -ld minimun number of links (default 0) " << endl;
    cout << "  -lu maximun number of links (default 10000) " << endl;
    cout << "  -multiFasta\t\tFILE  with all the reads " << endl;
    cout << "  -nt number of threads to use (default 1) " << endl;
    cout << "  -output name for the output-file  if print option was selected (default output)" << endl;
    cout << "  -print enable print output to file (default false) " << endl;
    cout << "  -sizeBin minimum number of reads to report a bin (default 1000) " << endl;
    cout << ""<< endl;

}

