/*  Binnig tool: Version que toma los links totales buenos para cortar el grafo(using FM_Index results)"
 *  Se usa and en lugar de or
 *
 *        -h --help
 *        -matrixAlignments  cleaning read FILE after blast
 *                    */


#include<iostream>
#include<cstdint>
#include<fstream>
#include<string>
#include <cstdlib>
#include <tr1/unordered_map>
#include <ctime>
#include <fcntl.h>
#include <string.h>
#include <sstream> 
#include <algorithm>
#include <omp.h>
#include <vector>
#include <math.h>
#define bufferSize 10000

//using std::unordered_map;
using std::tr1::unordered_map;
using namespace std;

struct Args {bool matrixAlignments; bool nt; bool links_Threshold; bool numberOFreads; bool tol; bool index; bool outputFile;};
struct Names {string matrixAlignments; string index; string outputFile;};
struct Parameters {int cpus; int numberOFreads; int linksThresholdDW; uint32_t linksThresholdUP; int sizeBin; float tol;bool dM;};
struct StaticsVals {float mean; float std; float std1; float cv; float median; float mode;float mad; float p; int min; int max; float tol; float cut;};

int mstrlen(const char arg[]);
void printerror(const char arg[]);
bool IsParam(char arg[],const char comp[]);

void binning(Names& names, Parameters *parameters, string *titles, int* queryList, vector<uint32_t>* MatrixList, string *lines);
void result2Matrix(int cpus,int sizeLoad,string *lines,int* queryList, vector<uint32_t>* MatrixList);
//void prinBin(Names& names, string *titles, int numBin, int *stack, int *get);
void loadIndex(Names& names, string *titles);//, vector<string> *title)

void readResultFile(Names& names, Parameters *parameters, vector<uint32_t>* MatrixList);
void filterReadsbyLinks(Parameters *parameters, vector<uint32_t>* MatrixList, uint32_t* queryList);
void grahTraversal(Names& names,Parameters *parameters, vector<uint32_t>* MatrixList, uint32_t* queryList, string *titles);
void static_parameters(StaticsVals *statics_vals,int numBin, uint32_t *stack, uint32_t size, vector<uint32_t> *MatrixList, uint32_t *queryList);
void printBin(Names& names, int numBin, uint32_t *stack, uint32_t *get, string *titles);

int main(int argc,char *argv[])
{
    //{int cpus; int numberOFreads; int linksThresholdDW; linksThresholdUP; int sizeBin; float tol;}
    Parameters parameters={1,0,0,100000000,1,0.5,false};
    bool   error = false, found1;
    Args   args = {false, false, false, false, false, false, false};
    Names  names;
    int    argum = 1;

    names.outputFile = "bin";
    
    while (argum < argc && not error) 
    {
        found1 = false; 
        if ( IsParam(argv[argum],"-rt") ) 
        {
            argum ++;
            if (argum < argc)  //verifica si el parametro existe
            {
                args.matrixAlignments = true;
                names.matrixAlignments = argv[argum];
                //cout << "matrix OK" << endl;
            } 
            else 
                error = true;
           argum ++;
           found1 = true;
           continue;
        }
        if ( IsParam(argv[argum],"-i") ) 
        {
            argum ++;
            if (argum < argc)  //verifica si el parametro existe
            {
                args.index = true;
                names.index = argv[argum];
                //cout << "matrix OK" << endl;
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
                args.outputFile = true;
                names.outputFile = argv[argum];
            } 
            else 
                error = true;
        
            argum ++;
            found1 = true;
            continue;
        }
        if (IsParam(argv[argum],"-n") ) 
        {
            argum ++;
            args.numberOFreads = true;
            found1 = true;
            if (argum < argc)  //verifica si el parametro existe
                parameters.numberOFreads = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-nt") ) 
        {
            argum ++;
            args.nt = true;
            found1 = true;
            if (argum < argc)  //verifica si el parametro existe
                parameters.cpus = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-lu") )
        {
            argum ++;
            args.links_Threshold = true;
            found1 = true;
            if (argum < argc)  //verifica si el parametro existe
                parameters.linksThresholdUP = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-ld") )
        {
            argum ++;
            args.links_Threshold = true;
            found1 = true;
            if (argum < argc)  //verifica si el parametro existe
                parameters.linksThresholdDW = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-sizeBin") )
        {
            argum ++;
            args.links_Threshold = true;
            found1 = true;
            if (argum < argc)  //verifica si el parametro existe
                parameters.sizeBin = atoi(argv[argum]);
            argum ++;
            continue;
        }
        if ( IsParam(argv[argum],"-tol") ) //MAD tol
        {
            argum ++;
            if (argum < argc)  
            {
                args.tol = true;
                std::string::size_type sz;     
                parameters.tol = std::stof (argv[argum],&sz);
           } 
            else 
                error = true;
         
            argum ++;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-dM") ) 
        {
            argum ++;
            parameters.dM = true;
            found1 = true;
            continue;
        }
        if ( IsParam(argv[argum],"-h") ) 
        {
            error = true;
            found1 = true;
            continue;
        }
        if (not found1)
        {
                        error = true;
        }
    }//termina el while


    found1 = false;
    if (not error && args.matrixAlignments) 
    {
        //string *lines = new string [bufferSize];
        string *titles = new string[parameters.numberOFreads];
        uint32_t *queryList = new uint32_t[parameters.numberOFreads] ();
        vector <uint32_t> *MatrixList= new vector <uint32_t> [parameters.numberOFreads];

	loadIndex(names, titles);
        readResultFile(names, &parameters, MatrixList);
        filterReadsbyLinks(&parameters, MatrixList, queryList);
        grahTraversal(names, &parameters, MatrixList, queryList,titles);

        //binning(names, &parameters, titles, queryList, MatrixList,lines);
        found1 = true;
        error = true;
    } 
  
    if (error && not found1)
        printerror(argv[0]);
  
    return 1;
} //fin del main

void loadIndex(Names& names, string *titles)//, vector<string> *title)
{
    //input file
    ifstream infile;
    infile.open(names.index.c_str());
    string line, name="";
    int i=0;
    if (infile.is_open()) 
    {
        while (getline(infile,line))
        {
            istringstream iss(line);
            getline(iss, name, '\t'); 
            titles[i++]=name;
            //cout<<name<<endl;
        }
    }
}

void readResultFile(Names& names, Parameters *parameters, vector<uint32_t>* MatrixList)
{
    //1. carga la matrix del alineamiento a mem
    int numberOFreads=parameters->numberOFreads; 
    int state=0;
    uint32_t  subject;
    uint32_t query;
        
    std::stringstream ss;  
    ss<<names.matrixAlignments;
    
    while(ss.good() )
    {
        string nameFile;
        getline(ss, nameFile, ',' );
        cout<<"Reading .result "<<nameFile<<"\n";
    
        ifstream myFile (nameFile, ios::in | ios::binary);

        while(myFile.read((char *)&subject,sizeof(subject)))
        {
            if (state==0)
            {
                query=subject;
                state=1;
            }
            else
            {
                if (subject!=4294967295)//0xFFFFFFFF)
                    MatrixList[query].push_back(subject);
                else
                    state=0;
            }
        }
        myFile.close();
    }
    
    /*for (int q=0; q<numberOFreads;q++)
    {
        cout<<q<<"\t";
        for( vector<uint32_t>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
        {
                
            cout<<*j<<"\t";
        }
        cout<<"\n";
    }*/
    
   
}

void filterReadsbyLinks(Parameters *parameters, vector<uint32_t>* MatrixList, uint32_t* queryList)
{
    int cpus=parameters->cpus;
    int numberOFreads=parameters->numberOFreads; 
    int ld=parameters->linksThresholdDW;
    int lu=parameters->linksThresholdUP;
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime)
        for(uint32_t k=0; k<numberOFreads;k++) 
            if(MatrixList[k].size()<= ld || MatrixList[k].size()>= lu)
                queryList[k]=3; //marked as  visited 
    }
}

void grahTraversal(Names& names, Parameters *parameters, vector<uint32_t>* MatrixList, uint32_t* queryList, string *titles)
{
    int cpus=parameters->cpus;
    uint32_t numberOFreads=parameters->numberOFreads; 
    bool dM=parameters->dM;
    int sizeBin=parameters->sizeBin;

    uint32_t *stack = new uint32_t[numberOFreads] ();
    uint32_t *put=stack;
    uint32_t *get=stack;
    int numBin=0; 

    
    //StaticsVals{float mean; float std; float std1; float cv; float median; float mode;float mad; float p; int min; int max; float tol; float cut;};
    StaticsVals statics_vals={0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,2147483647,0.5,3.0};
    statics_vals.tol=parameters->tol;
    
    for(uint32_t i=0; i<numberOFreads;i++)
    {
        if(queryList[i]==0) //not visited
        {
            queryList[i]=1;
            *(put++)=i;
            while(get!=put)
            {
                uint32_t q=*get++;
                for( vector<uint32_t>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                {
                    if(queryList[*j]==0) //not visited
                    {
                        queryList[*j]=1;
                        *(put++)=*j;
                    }
                }
            }
            
            //Outliers for the BIN
            uint32_t size=get-stack;
            while(dM==false && size>=sizeBin && (abs(statics_vals.p-1.0)>statics_vals.tol) && statics_vals.cut>0 && statics_vals.p>1.0)
            {
               static_parameters(&statics_vals,numBin,stack,size,MatrixList,queryList);
               for(uint32_t k=0; k<size;k++) 
               {
                    uint32_t s=*(stack+k);
                    if(MatrixList[s].size()>statics_vals.max || MatrixList[s].size()<statics_vals.min)
                        queryList[s]=0; //marked as  non taked queryList[s]=3
                    else
                        queryList[s]=1;
                }
                
                //Redefine the binning
                put=stack;
                get=stack;
                for(uint32_t k=0; k<size;k++)
                {
                    uint32_t s=*(stack+k);
                    if(queryList[s]==1) //no binned
                    {
                        queryList[s]=2;
                        *(put++)=s;
                        while(get!=put)
                        {
                            uint32_t q=*get++;
                            for( vector<uint32_t>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                            {
                                if(queryList[*j]==1) //not binned
                                {
                                    queryList[*j]=2;
                                    *(put++)=*j;
                                }
                            }
                        }
                    }
                }
                //cout<<"\tFin redefine"<<endl;
                size=get-stack;
                statics_vals.cut=statics_vals.cut-0.25;
            }
            //cout<<"salio outliers"<<endl;
            
            //cout<<"imprimir"<<endl;
            if (size>=sizeBin)
            {    
                cout<<">bin_"<<numBin<<"\t"<<size<<"\t"<<statics_vals.mean<<"\t"<<statics_vals.std<<"\t"<<statics_vals.median<<"\t"<<statics_vals.mad<<"\t"<<statics_vals.p<<"\t"<<statics_vals.cut<<"\t"<<statics_vals.min<<"\t"<<statics_vals.max<<endl;
                
                printBin(names, numBin, stack, get, titles);
                //while(get!=stack)
                    //cout<<*(--get)<<endl;
                numBin++;
               
            }    
            //cout<<"fin imprimir"<<endl;
            statics_vals.p=3.0;
            statics_vals.cut=3.0;
            put=stack; //pointers restart
            get=stack; 
        }
    }
    cout<<"Binning end"<<endl;
}

void static_parameters(StaticsVals *statics_vals,int numBin, uint32_t *stack, uint32_t size, vector<uint32_t> *MatrixList, uint32_t *queryList) 
{
    float mean=0.0;
    float std,std1=0.0;
    float cv;      
    float median;
    float mode;
    
    vector <float> v,v2;
    vector <int> v3;
    
    //mean
    //cout<<"mean"<<endl;
    float links=0.0;
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        links+=edges;
        v.push_back(edges);
    }
    mean=(float)(links/size);

    //std
    //cout<<"std"<<endl;
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        std1+=pow(edges - mean, 2);
    }
    std=sqrt(std1/size);   
    
    //variation coefficient
    cv=100*std/mean;          
    
    //median
    sort(v.begin(), v.end());
    median = size % 2 ? v[size/2] : (v[size/2-1] + v[size/2]/2);
    
    //mad
    //cout<<"mad"<<endl;
    for (int i=0; i<size; i++)
    {
        float aux=pow(v[i]-median, 2);
        v2.push_back(sqrt(aux));
        //cout<<"aux:"<<aux<<endl;
    }
    sort(v2.begin(), v2.end());
    float median2 = size % 2 ? v2[size/2] : (v2[size/2-1] + v2[size / 2] / 2);
    float mad = 1.4826*median2;   
    
    //cutoff
    float cutoff, aux;
    //cout<<"cutoff"<<endl;
    for (int i=0; i<size; i++)
    {
        int edges=MatrixList[*(stack+i)].size();
        float aux=sqrt(pow(edges-median, 2));  
        cutoff=0;
        if (mad>0)
            cutoff= aux/mad;
        //cout<<">bin"<<numBin<<"\t"<<size<<"\t"<<*(stack+i)<<"\t"<<edges<<"\t"<<cutoff<<endl;
        if(cutoff<statics_vals->cut)
             v3.push_back(edges);
    }
    
    //sort(v3.begin(), v3.end());
    //cout<<"max min:"<<v2[size/2-1]<<"->"<<median<<"->"<<size<<"->"<<v2.size()<<"->"<<median2<<"->"<<cutoff<<"->"<<v3.size()<<endl;
    int max=v3[0]; int min=v3[0];
    for (std::vector<int>::iterator it = v3.begin() ; it != v3.end(); ++it)
    {
        if(*it>max)
            max=*it;
        
        if (*it<min)
            min=*it;
            
    }
    //cout<<"END max min"<<endl;

    float p=1.0;
    if (mean>0)
      p=3*std/mean;
    //cout<<"p="<<p<<endl;

    statics_vals->mean=mean;
    statics_vals->std=std;
    statics_vals->std1=std1;
    statics_vals->cv=cv;      
    statics_vals->median=median;
    statics_vals->mode=mode;
    statics_vals->mad=mad;
    statics_vals->p=p;
    statics_vals->min=min;
    statics_vals->max=max;
    
    /*cout<<">bin_"<<numBin<<"\t"<<size<<"\t"<<mean<<"\t"<<std<<"\t"<<median<<"\t"<<mad<<"\t"<<p<<"\t"<<min<<"\t"<<max<<endl;
    for( int p=0;p<size;p++)
        cout<<*(stack+p)<<endl;*/
    
}
void printBin(Names& names,  int numBin, uint32_t *stack, uint32_t *get, string *titles)
{
    //binned reads' list to be exclude in future stages
    ofstream myfileEx;
    stringstream ss;
    string nameFile=names.outputFile+".exclude";
    myfileEx.open(nameFile.c_str(), std::fstream::app);
    
    
    //output file
    ofstream myfile;
    stringstream ss2;
    ss2 << numBin;
    string strNumBin = ss2.str();
    string nameFile2=names.outputFile+"_"+strNumBin+".list";
    myfile.open(nameFile2.c_str());

    while(get!=stack)
    {
        --get;
        myfileEx<<*get<<endl;
        //myfile<<*get<<endl;
        myfile<<titles[*get]<<endl;
    }
    myfileEx.close();
    myfile.close();
}
/*
void result2Matrix(int cpus,int sizeLoad,string *lines,int* queryList, vector<uint32_t>* MatrixList)
{
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime)
        for (int i=0; i<sizeLoad; i++)
        {
            string item;
            istringstream iss(lines[i]);
            getline(iss, item, '\t'); //iss >> item;
            uint32_t query=atoi(item.c_str());
            queryList[query]=0;
            while(getline(iss, item, '\t'))//while(iss >> item)
            {
                uint32_t subject=atoi(item.c_str());
                MatrixList[query].push_back(subject);
                //MatrixList[subject].push_back(query);
                //cout<< query <<","<<subject<<endl;
            }
            //cout<<"------------"<<endl;
        }
    }
   
}
void binning(Names& names, Parameters *parameters, string *titles, int* queryList, vector<uint32_t>* MatrixList, string *lines)
{
    int cpus=parameters->cpus;
    int ld=parameters->linksThresholdDW;
    int sizeBin=parameters->sizeBin;
    int numberOFreads=parameters->numberOFreads; 
    bool dM=parameters->dM;
    int *stack = new int[numberOFreads] ();
    int *put=stack;
    int *get=stack;
    int numBin=0; 
    //StaticsVals{float mean; float std; float std1; float cv; float median; float mode;float mad; float p; int min; int max; float tol; float cut;};
    StaticsVals statics_vals={0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0,2147483647,0.5,3.0};
    statics_vals.tol=parameters->tol;
    
    //1. carga la matrix del alineamiento a mem
    //string lines [numberOFreads];
    string item;
    ifstream in;
    std::ostringstream ss;
    ss<<names.matrixAlignments;
    string nameFile=ss.str();
    in.open(nameFile.c_str());

    //cout<<"Reading .result "<<names.matrixAlignments<<"\n";
    int i=0; int j=0;
    bool full=false;
    while (getline(in,lines[i++]).good()) 
    {
        full=false;
        if (i==bufferSize)
        {
            i=0;
            j++;
            result2Matrix(cpus,bufferSize,lines,queryList,MatrixList);
            full=true;
        }
    }
    if  (full==false)
        result2Matrix(cpus,i-1,lines,queryList,MatrixList);

    //cout<<"File loaded "<<names.matrixAlignments<<"\t"<<j*bufferSize+i-1<<"\n";
    in.close();  

    
    //2.Load index
    loadIndex(names, titles);

    //3. Genera los bins
    //cout<<"Binning "<<names.matrixAlignments<<"\n";
    //cout<<">bin\tsize\tmean\tstd\tmedian\tMAD\tp=3std/mean\tdist2mean\tOutliers\n";

    //2.1.Number of links by read
    //cout<<"removiendo nodos con pocos links"<<endl;
    #pragma omp parallel num_threads(cpus)
    {
        #pragma omp for schedule(runtime)
        for(uint32_t k=0; k<numberOFreads;k++) 
            if(MatrixList[k].size()<= ld)
                queryList[k]=3; //marked as  visited 
    }
    //cout<<"End links"<<endl;
      
    for(int i=0; i<numberOFreads;i++)
    {
        if(queryList[i]==0) //not visited
        {
            //cout<<"En el grafo:"<<i<<endl;
            queryList[i]=1;
            *(put++)=i;
            while(get!=put)
            {
                uint32_t q=*get++;
                for( vector<uint32_t>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                {
                    if(queryList[*j]==0) //not visited
                    {
                        queryList[*j]=1;
                        *(put++)=*j;
                    }
                }
            }
            //cout<<"Salio del grafo"<<endl;
            
            //Outliers for the BIN
            //cout<<"En outliers"<<endl;
            uint32_t size=get-stack;
            while(dM==false && size>=sizeBin && (abs(statics_vals.p-1.0)>statics_vals.tol) && statics_vals.cut>0 && statics_vals.p>1.0)
            {
               //cout<<"\tEn stadistics"<<endl;
               static_parameters(&statics_vals,numBin,stack,size,MatrixList,queryList);
               for(uint32_t k=0; k<size;k++) 
               {
                    uint32_t s=*(stack+k);
                    if(MatrixList[s].size()>statics_vals.max || MatrixList[s].size()<statics_vals.min)
                        queryList[s]=0; //marked as  non taked queryList[s]=3
                    else
                        queryList[s]=1;
                }
               //cout<<"\tfin stadistics"<<endl;
                
                //Redefine the binning
                //cout<<"\tredefine"<<endl;
                put=stack;
                get=stack;
                for(uint32_t k=0; k<size;k++)
                {
                    uint32_t s=*(stack+k);
                    if(queryList[s]==1) //no binned
                    {
                        queryList[s]=2;
                        *(put++)=s;
                        while(get!=put)
                        {
                            uint32_t q=*get++;
                            for( vector<uint32_t>::iterator j=MatrixList[q].begin(); j!=MatrixList[q].end(); ++j)
                            {
                                if(queryList[*j]==1) //not binned
                                {
                                    queryList[*j]=2;
                                    *(put++)=*j;
                                }
                            }
                        }
                    }
                }
                //cout<<"\tFin redefine"<<endl;
                size=get-stack;
                statics_vals.cut=statics_vals.cut-0.25;
            }
            //cout<<"salio outliers"<<endl;
            
            //cout<<"imprimir"<<endl;
            if (size>=sizeBin)
            {    
                cout<<">bin_"<<numBin<<"\t"<<size<<"\t"<<statics_vals.mean<<"\t"<<statics_vals.std<<"\t"<<statics_vals.median<<"\t"<<statics_vals.mad<<"\t"<<statics_vals.p<<"\t"<<statics_vals.cut<<"\t"<<statics_vals.min<<"\t"<<statics_vals.max<<endl;
                
                prinBin(names, titles, numBin, stack, get);

                //while(get!=stack)
                    //cout<<*(--get)<<endl;
                numBin++;
            }    
            //cout<<"fin imprimir"<<endl;
            statics_vals.p=3.0;
            statics_vals.cut=3.0;
            put=stack; //pointers restart
            get=stack; 
        }
    }
    cout<<"Binning end"<<endl;
}

void static_parameters(StaticsVals *statics_vals,int numBin, int *stack, int size, vector<uint32_t> *MatrixList, int *queryList) 
{
    float mean=0.0;
    float std,std1=0.0;
    float cv;      
    float median;
    float mode;
    
    vector <float> v,v2;
    vector <int> v3;
    
    //mean
    //cout<<"mean"<<endl;
    float links=0.0;
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        links+=edges;
        v.push_back(edges);
    }
    mean=(float)(links/size);

    //std
    //cout<<"std"<<endl;
    for( int p=0;p<size;p++)
    {
        int edges=MatrixList[*(stack+p)].size();
        std1+=pow(edges - mean, 2);
    }
    std=sqrt(std1/size);   
    
    //variation coefficient
    cv=100*std/mean;          
    
    //median
    sort(v.begin(), v.end());
    median = size % 2 ? v[size/2] : (v[size/2-1] + v[size/2]/2);
    
    //mad
    //cout<<"mad"<<endl;
    for (int i=0; i<size; i++)
    {
        float aux=pow(v[i]-median, 2);
        v2.push_back(sqrt(aux));
        //cout<<"aux:"<<aux<<endl;
    }
    sort(v2.begin(), v2.end());
    float median2 = size % 2 ? v2[size/2] : (v2[size/2-1] + v2[size / 2] / 2);
    float mad = 1.4826*median2;   
    
    //cutoff
    float cutoff, aux;
    //cout<<"cutoff"<<endl;
    for (int i=0; i<size; i++)
    {
        int edges=MatrixList[*(stack+i)].size();
        float aux=sqrt(pow(edges-median, 2));  
        cutoff=0;
        if (mad>0)
            cutoff= aux/mad;
        //cout<<">bin"<<numBin<<"\t"<<size<<"\t"<<*(stack+i)<<"\t"<<edges<<"\t"<<cutoff<<endl;
        if(cutoff<statics_vals->cut)
             v3.push_back(edges);
    }
    
    //sort(v3.begin(), v3.end());
    //cout<<"max min:"<<v2[size/2-1]<<"->"<<median<<"->"<<size<<"->"<<v2.size()<<"->"<<median2<<"->"<<cutoff<<"->"<<v3.size()<<endl;
    int max=v3[0]; int min=v3[0];
    for (std::vector<int>::iterator it = v3.begin() ; it != v3.end(); ++it)
    {
        if(*it>max)
            max=*it;
        
        if (*it<min)
            min=*it;
            
    }
    //cout<<"END max min"<<endl;

    float p=1.0;
    if (mean>0)
      p=3*std/mean;
    //cout<<"p="<<p<<endl;

    statics_vals->mean=mean;
    statics_vals->std=std;
    statics_vals->std1=std1;
    statics_vals->cv=cv;      
    statics_vals->median=median;
    statics_vals->mode=mode;
    statics_vals->mad=mad;
    statics_vals->p=p;
    statics_vals->min=min;
    statics_vals->max=max;
    
    cout<<">bin_"<<numBin<<"\t"<<size<<"\t"<<mean<<"\t"<<std<<"\t"<<median<<"\t"<<mad<<"\t"<<p<<"\t"<<min<<"\t"<<max<<endl;
    for( int p=0;p<size;p++)
        cout<<*(stack+p)<<endl;
    
}

void loadIndex(Names& names, string *titles)//, vector<string> *title)
{
    //input file
    ifstream infile;
    infile.open(names.index.c_str());
    string line, name="";
    int i=0;
    if (infile.is_open()) 
    {
        while (getline(infile,line))
        {
            istringstream iss(line);
            getline(iss, name, '\t'); 
            titles[i++]=name;
            //cout<<name<<endl;
        }
    }
}
void prinBin(Names& names, string *titles, int numBin, int *stack, int *get)
{
    //binned reads' list to be exclude in future stages
    ofstream myfileEx;
    stringstream ss;
    string nameFile=names.outputFile+".exclude";
    myfileEx.open(nameFile.c_str(), std::fstream::app);
    
    
    //output file
    ofstream myfile;
    stringstream ss2;
    ss2 << numBin;
    string strNumBin = ss2.str();
    string nameFile2=names.outputFile+"_"+strNumBin+".list";
    myfile.open(nameFile2.c_str());

    while(get!=stack)
    {
        --get;
        myfileEx<<*get<<endl;
        myfile<<titles[*get]<<endl;
    }
    myfileEx.close();
    myfile.close();
}
*/
int mstrlen(const char arg[]) {
  int out = 0;
  while(arg[out] != 0)
    out++;
  return out;
}

bool IsParam(char arg[],const char comp[])
{
  int len_comp = mstrlen(comp);
  bool equal = true;
  if (mstrlen(arg) == len_comp) {
    for (int i = 0; i < len_comp; i++) {
      if (arg[i] != comp[i]) {
        equal = false;
      }
    }
  } else {
    equal = false;
  }
  return equal;
}



void printerror(const char arg[])
{
    cout << "Binnig tool: It uses a number of bases threshold and links profile to generate the bins" << endl;
    cout << "Binnig tool: Version 5, Using total well-score links, whose agree with the link threshold (and onlybigger than)" << endl;
    cout << endl;
    cout << arg << endl;
    cout << "  -h               (Help)" << endl;  
    cout << "  -dM Disable MAD processs (default enabled)"<< endl;  
    cout << "  -i  Index file with reads name" << endl;  
    cout << "  -n  number of reads" << endl;  
    cout << "  -nt number of threads to use" << endl;  
    cout << "  -lu number of links to cut by UP threshold" << endl;  
    cout << "  -ld number of links to cut by Down threshold" << endl;  
    cout << "  -rt result FILE (comma delimted for several files)" << endl;
    cout << "  -sizeBin number of reads to report a bin (default 1)" << endl;
    cout << "  -tol MAD error tolerance (default 0.5) " << endl;
    cout << ""<< endl;
}


