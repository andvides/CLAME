/*---------------------------------------------------------------
 *
 *   CLAME:"mapping"
 *      
 *   
 *  
 ---------------------------------------------------------------*/

#include "mapping_lite.h"

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
        if ( IsParam(argv[argum],"-offsetFM9") ) //superior cut
        {
            argum ++;
            if (argum < argc)  
            {
                args->offsetFM9 = true;
            parameters->offsetFM9 = atoi(argv[argum]);
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
        if ( IsParam(argv[argum],"-w") ) //minimal reads set to report a bin
        {
            argum ++;
            if (argum < argc)  
            {
                args->w = true;
                parameters->w = atoi(argv[argum]);
            } 
            else 
                error = true;
                       
            argum ++;
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
        if ( IsParam(argv[argum],"-list2Exclude") ) //to exclude reads in list file
        {
            argum ++;
            if (argum < argc)  
            {
                args->list2Exclude = true;
                names->list2Exclude = argv[argum];
            } 
            else 
                error = true;
                       
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-print") ) 
        {
            argum ++;
            args->print = true;
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

bool readList2Exclude(Names *names, int *reads2exclude)
{
    ifstream infile;
    infile.open(names->list2Exclude.c_str());
    if (infile.is_open()) 
    {
        string word;
        int i=0;
        while (getline(infile,word))
        {
            reads2exclude[atoi(word.c_str())]=2; //to indecate that it should be excluded
            //cout<<word<<endl;
        }
        infile.close();
        return false;
    }
    else
    {
        #ifdef debug
            std::cout << "Error opening readList2Exclude file";
        #endif
        return true;
    }
}

bool readFasta2(Names *names, vector<string> *bases,Parameters *parameters, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index, vector<uint64_t> *index)
{
    ofstream myfileResult, myfileLinks;
    string nameLinks=names->outputFile+".links";//ss.str();
    myfileLinks.open(nameLinks.c_str());//, ios::out | ios::app);ios::out | ios::app | ios::binary);
    
    string nameResult=names->outputFile+".resultb";//ss.str();
    myfileResult.open(nameResult.c_str(), ios::out | ios::binary);
    
    //to print index
    ofstream myfile;
    string nameFile=names->outputFile+".index";
    myfile.open(nameFile.c_str());

    ifstream infile;
    infile.open(names->multiFasta.c_str());
    if (infile.is_open()) 
    {
        string word, name="", read="";
        bool firstTime=true;
        int indexTitle=0;
        int countReads=0;
        int block=0;
        int offset=0;
        while (getline(infile,word))
        {
            if (word.find(">") != string::npos) 
            {
                if (firstTime)
                    firstTime=false;
                else
                {
                    countReads++;
                    myfile<<name<<"\t"<<indexTitle++<<endl; 
                    bases->push_back(read);
                }
                int found = word.find(' ');
                name=word.substr (0,found);
            }
            else
                read=read+word;
            
            if(countReads==parameters->size)         //generate the FM9
            {
                offset=block*parameters->size;
                alignemnt2(names, parameters, bases,  index, MatrixList, reads2exclude, ptrfm_index,offset);
                printResult2(names, parameters,countReads, MatrixList, reads2exclude, offset,  &myfileResult, &myfileLinks);
                block++;
                bases->clear();
                countReads=0;
            }
        }
        
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        bases->push_back(read);
        
        if(countReads)         //generate the FM9
        {
            offset=block*parameters->size;
            parameters->size=countReads;
            alignemnt2(names, parameters, bases,  index, MatrixList, reads2exclude, ptrfm_index,offset);
            printResult2(names, parameters,countReads, MatrixList, reads2exclude, offset,  &myfileResult, &myfileLinks);
            block++;
            bases->clear();
            countReads=0;
        }
        
       
        infile.close();
        myfile.close();
        myfileResult.close();
        myfileLinks.close();
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

bool readFastQFile2(Names *names, vector<string> *bases,Parameters *parameters, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index, vector<uint64_t> *index)
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
    
    ofstream myfileResult, myfileLinks;
    string nameLinks=names->outputFile+".links";//ss.str();
    myfileLinks.open(nameLinks.c_str());//, ios::out | ios::app);ios::out | ios::app | ios::binary);
    
    string nameResult=names->outputFile+".resultb";//ss.str();
    myfileResult.open(nameResult.c_str(), ios::out  | ios::binary);
    


    if (infile.is_open()) 
    {
        string word, name="", read="", quality="";
        bool firstTime=true;
        int countLine=0;
        int countReads=0;
        int block=0;
        int offset=0;
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                    {    
                        int found=word.find(' '); 
                        name=word.substr(0,found);
                        countLine=0; state=s1; myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                        countReads++;
                    }
                    break;
                case s1:
                    if (word[0]=='+') 
                        state=s2;
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
                        if(countReads==parameters->size)         
                        {
                            offset=block*parameters->size;
                            alignemnt2(names, parameters, bases,  index, MatrixList, reads2exclude, ptrfm_index,offset);
                            printResult2(names, parameters,countReads, MatrixList, reads2exclude, offset, &myfileResult, &myfileLinks);
                            block++;
                            bases->clear();
                            countReads=0;
                        }
                        read="";
                        state=s0;
                    }
                    break;
            }
        }
        
        if(countReads)         
        {
            offset=block*parameters->size;
            parameters->size=countReads;
            alignemnt2(names, parameters, bases,  index, MatrixList, reads2exclude, ptrfm_index,offset);
            printResult2(names, parameters,countReads, MatrixList, reads2exclude, offset, &myfileResult, &myfileLinks);
            block++;
            bases->clear();
            countReads=0;
        }  
        infile.close();
        myfile.close();
        myfileResult.close();
        myfileLinks.close();
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
bool readFasta(Names *names, vector<string> *bases)
{
    
    //to print index
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

        while (getline(infile,word))
        {
            if (word.find(">") != string::npos) 
            {
                if (firstTime)
                    firstTime=false;
                else
                {
                    myfile<<name<<"\t"<<indexTitle++<<endl; 
                    bases->push_back(read);
                    read+='$';
                    n = read.length();
                    indexBases+=n;
                    read="";
                }
                int found = word.find(' ');
                name=word.substr (0,found);
            }
            else
                read=read+word;
        }
        myfile<<name<<"\t"<<indexTitle++<<endl; 
        bases->push_back(read);
        read+='$';
        n = read.length();
        indexBases+=n;
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

bool readFastQFile(Names *names,  vector<string>* bases)
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
        string word, name="", read="", quality="";
        bool firstTime=true;
        int countLine=0;
        while (getline(infile,word))
        {
            switch ( state ) 
            {
                case s0:
                    if (word[0]=='@') 
                    {    
                        int found=word.find(' '); 
                        name=word.substr(0,found);
                        countLine=0; state=s1; myfile<<name<<"\t"<<indexTitle++<<endl; //print Read0 0...;
                    }
                    break;
                case s1:
                    if (word[0]=='+') 
                        state=s2;
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

bool loadFM9(Names *names, Parameters *parameters, vector<uint64_t> *index, csa_wt<> *ptrfm_index)
{
    ofstream  myfile2;
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile2.open(nameFile2.c_str());
       
   
    //FM-Index
    string index_file =""; 
    if(parameters->loadFM9)//load from file
    {

        index_file   = names->fm9;
        if (!load_from_file(*ptrfm_index, index_file)) 
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
            auto locations = locate(*ptrfm_index, querySearch.begin(), querySearch.begin()+1);
            sort(locations.begin(), locations.end());
            int i=0;
            for (auto it=locations.begin(); it<locations.end(); it++)
            {
                myfile2<<i++<<"\t"<<*it<<endl; 
                index->push_back(*it);
 
            }
        }
    }
    else//buidl suffix tree
    {
    
            #ifdef debug
                cout << "No index file found" << endl;
            #endif
            return true;
    }
    myfile2.close();//exit(10);
    cout << "LOAD OK" << endl;
}
bool alignemnt2(Names *names,Parameters *parameters,vector<string> *bases, vector<uint64_t> *index, vector<int>* MatrixList, int *reads2exclude, csa_wt<> *ptrfm_index, int offset)
{
    //search
    std::vector<string>::iterator ptrBases= bases->begin();
    bool runningError=false;
    int cpus=parameters->numThreads;
    int numberOFreads=parameters->size;
    int offsetFM9=parameters->offsetFM9;
    
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime) 
        for(int i=0; i<numberOFreads;i++)//process all parts except the last one because it can get a diffent size
        {
                
            if (reads2exclude[i+offset]!=2)
            {
                uint32_t subjectID=0;
                string query=*(ptrBases+i);
                int w=parameters->w;
                if (query.length()>=w+parameters->query_size)
                {
                    string reverseQuery=reverse(&query);
                    
                    string queryF1=query.substr(w,parameters->query_size);//forward init
                    string queryF2=query.substr((query.length()-w-parameters->query_size),parameters->query_size);;//forward end
                    string queryR1=reverseQuery.substr(w,parameters->query_size);//reverse init
                    string queryR2=reverseQuery.substr((reverseQuery.length()-w-parameters->query_size),parameters->query_size);//reverse end
    
                    string Queries[4]={queryF1,queryF2,queryR1,queryR2};
                    string orientations[4]={"5a'","5b'","3a'","3b'"};

                    for (int q=0;q<4;q++)
                    {
                        string querySearch=Queries[q];
                        uint32_t sizeq = querySearch.length();
                        auto locations = locate(*ptrfm_index, querySearch.begin(), querySearch.begin()+sizeq);
                        for (auto it=locations.begin(); it<locations.end(); it++)
                        {
                            subjectID=indx2Loc(*it,index);
                            if (reads2exclude[subjectID]!=2)//(it2 == reads2exclude->end())
                                MatrixList[i].push_back(subjectID+offsetFM9);
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
    cout << "Aligment OK" << endl;

    return runningError;
      
}
bool alignemnt(Names *names,Parameters *parameters,vector<string> *bases, vector<uint64_t> *index, vector<int>* MatrixList, int *reads2exclude)
{
    ofstream  myfile2;
    string nameFile2=names->outputFile+".baseIndex";//ss.str();
    myfile2.open(nameFile2.c_str());
       
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
            int i=0;
            for (auto it=locations.begin(); it<locations.end(); it++)
            {
                myfile2<<i++<<"\t"<<*it<<endl; 
                index->push_back(*it);
                /*#ifdef debug
                    cout<<i++<<"\t"<<*it<<endl;
                #endif*/
            }
        }
    }
    else//buidl suffix tree
    {
    
            #ifdef debug
                cout << "No index file found" << endl;
            #endif
            return true;
    }
    myfile2.close();//exit(10);
    cout << "LOAD OK" << endl;
    
    //search
    int cpus=parameters->numThreads;
    int numberOFreads=bases->size();
    int numBlocks=numberOFreads/1;
    bool runningError=false;
    int offsetFM9=parameters->offsetFM9;
    if(runningError==false)
    {
        #pragma omp parallel num_threads(cpus)
        {
            #pragma omp for schedule(runtime) 
            for(int i=0; i<numberOFreads;i++)//process all parts except the last one because it can get a diffent size
            {
                //std::vector<int>::iterator it0 = std::find(reads2exclude->begin(), reads2exclude->end(), i);
                if (reads2exclude[i]!=2)
                {
                    uint32_t subjectID=0;
                    string query=*(ptrBases+i);
                    int w=parameters->w;
                    if (query.length()>=w+parameters->query_size)
                    {
                        string reverseQuery=reverse(&query);
                    
                        string queryF1=query.substr(w,parameters->query_size);//forward init
                        string queryF2=query.substr((query.length()-w-parameters->query_size),parameters->query_size);;//forward end
                        string queryR1=reverseQuery.substr(w,parameters->query_size);//reverse init
                        string queryR2=reverseQuery.substr((reverseQuery.length()-w-parameters->query_size),parameters->query_size);//reverse end
    
                        string Queries[4]={queryF1,queryF2,queryR1,queryR2};
                        string orientations[4]={"5a'","5b'","3a'","3b'"};

                        for (int q=0;q<4;q++)
                        {
                            string querySearch=Queries[q];
                            uint32_t sizeq = querySearch.length();
                            auto locations = locate(fm_index, querySearch.begin(), querySearch.begin()+sizeq);
                            for (auto it=locations.begin(); it<locations.end(); it++)
                            {
                                subjectID=indx2Loc(*it,index);
                                //std::vector<int>::iterator it2 = std::find(reads2exclude->begin(), reads2exclude->end(), subjectID);
                                if (reads2exclude[subjectID]!=2)//(it2 == reads2exclude->end())
                                    MatrixList[i].push_back(subjectID+offsetFM9);
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
        cout << "Aligment OK" << endl;

        //print the result
        printResult(names,parameters,numberOFreads,MatrixList,reads2exclude);
        return runningError;
    }
}

void printResult2(Names *names, Parameters *parameters,int numberOFreads, vector<int> *MatrixList,int *reads2exclude, int offset, ofstream * myfileResult, ofstream * myfileLinks)
{
    cout << "print start" << endl;

    //ofstream myfile, myfile2;
    //string nameFile2=names->outputFile+".links";//ss.str();
    //myfile2.open(nameFile2.c_str(), ios::app);//, ios::out | ios::app);ios::out | ios::app | ios::binary);
    //if (parameters->enablePrint)
    //{
        //string nameFile=names->outputFile+".resultb";//ss.str();
        //myfile.open(nameFile.c_str(), ios::out | ios::app | ios::binary);
    //}

    for(int k=0;k<numberOFreads;k++)
    {
        if (reads2exclude[k]!=2)//(it0 == reads2exclude->end())
        {
                //print links
                *myfileLinks<<k+offset<<"\t"<<MatrixList[k].size()<<endl;
                //print result
                int i=k+offset;
                myfileResult->write(reinterpret_cast<const char*>(&i), sizeof(i));//myfile<<k+offset<<"\t";
                for (auto it=MatrixList[k].begin(); it<MatrixList[k].end(); it++)
                { 
                  int j=*it;
                  myfileResult->write(reinterpret_cast<const char*>(&j), sizeof(j));//myfile<<*it<<"\t";
                }
                int j=4294967295;
                myfileResult->write(reinterpret_cast<const char*>(&j), sizeof(j));//myfile<<"\n";
                MatrixList[k].clear();
        }
    }
    //myfile2.close();
    //if (parameters->enablePrint)
        //myfile.close();
}
void printResult(Names *names, Parameters *parameters,int numberOFreads, vector<int> *MatrixList,int *reads2exclude)
{
    cout << "print start" << endl;

    ofstream myfile, myfile2;
    string nameFile2=names->outputFile+".links";//ss.str();
    //myfile2.open(nameFile2.c_str(), ios::out | ios::app | ios::binary);
    myfile2.open(nameFile2.c_str());//, ios::out | ios::app);
    if (parameters->enablePrint)
    {
        //ofstream myfile;
        string nameFile=names->outputFile+".result";//ss.str();
        myfile.open(nameFile.c_str());
    }

    for(int k=0;k<numberOFreads;k++)
    {
        //if(MatrixList[k].size()>=parameters->ld && MatrixList[k].size()<=parameters->lu )
        //std::vector<int>::iterator it0 = std::find(reads2exclude->begin(), reads2exclude->end(), k);
        if (reads2exclude[k]!=2)//(it0 == reads2exclude->end())
        {
                //print links
                //myfile2<<*(ptrTitle+k)<<"\t"<<MatrixList[k].size()<<endl;
                myfile2<<k<<"\t"<<MatrixList[k].size()<<endl;
                //print result
                //if (parameters->enablePrint)
                {
                        myfile<<k<<"\t";
                        for (auto it=MatrixList[k].begin(); it<MatrixList[k].end(); it++)
                                myfile<<*it<<"\t";
                        myfile<<"\n";
                }
        }
    }
    myfile2.close();
    if (parameters->enablePrint)
        myfile.close();
}

uint32_t indx2Loc(uint64_t locs, vector<uint64_t> *index)
{
        std::vector<uint64_t>::iterator up;
        uint64_t a=locs+1;
        up= std::upper_bound (index->begin(), index->end(), locs); 
        return (up-1-index->begin());
}


string reverse(string *str)
{
        string rs="";
        for (std::string::reverse_iterator rit=str->rbegin(); rit!=str->rend(); ++rit)
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
    cout << "CLAME: mapping to reference" << endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h\t\t\t(Help)" << endl;
    cout << "  -b minimum number of bases to take an alignment (default 20) " << endl;
    cout << "  -fm9 Load fm9 file  " << endl;
    cout << "  -fastq input file is in a fastq format  " << endl;
    cout << "  -list2Exclude file with sequeces to exclude of the Aligment" <<endl;
    cout << "  -multiFasta\t\tFILE  with all the reads " << endl;
    cout << "  -nt number of threads to use (default 1) " << endl;
    cout << "  -offsetFM9 when use several FM9 indexes (default 0) " << endl;
    cout << "  -output name for the output-file  if print option was selected (default output)" << endl;
    cout << "  -print print the result file (default false)" <<endl;
    cout << "  -size size for the block aligment (default all dataset)" <<endl;
    cout << "  -w windows offset to start the alignemnt (default 0) "<<endl;
    cout << ""<< endl;
}

