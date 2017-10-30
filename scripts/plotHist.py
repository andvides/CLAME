#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os,argparse,math
from subprocess import Popen, PIPE
from pylab import *


if sys.argv.count("-f")and sys.argv.count("-c") :  #chequea que la opcion -src haya sido introducida

    column=int (sys.argv[sys.argv.index("-c")+1])

    #Graph options
    if sys.argv.count("-lu"):
        bins=int (sys.argv[sys.argv.index("-lu")+1])
    else:
        bins=1000;

    if sys.argv.count("-ld"):
        binsINF=int (sys.argv[sys.argv.index("-ld")+1])
    else:
        binsINF=0;

    if sys.argv.count("-title"):
        title=(sys.argv[sys.argv.index("-title")+1])
    else:
        title=" "

    if sys.argv.count("-xlabel"):
        xlabel=(sys.argv[sys.argv.index("-xlabel")+1])
    else:
        xlabel="Links"

    if sys.argv.count("-ylabel"):
        ylabel=(sys.argv[sys.argv.index("-ylabel")+1])
    else:
        ylabel="Number of reads"     

    if sys.argv.count("-format"):
        format=str(sys.argv[sys.argv.index("-format")+1])
    else:
        format='png'
        
    if sys.argv.count("-label1"):
        label1=sys.argv[sys.argv.index("-label1")+1]
    else:
        label1=' ';    

    #Main plot
    hist=[0]*(bins+1)
    filename=sys.argv[sys.argv.index("-f")+1]
    fsrc1 = open(sys.argv[sys.argv.index("-f")+1],'r') 
    for line in fsrc1:
        words=line.split("\n") 
        words=words[0] 
        words=words.split("\t")
        data=int(words[column])
        if data<bins and data>binsINF:
            hist[data]+=1
    fsrc1.close()
    
    if sys.argv.count("-line"):
        p1,=plt.plot(range(len(hist)),hist,color='r')
    else:
        p1=plt.bar(range(len(hist)),hist,color='r')

    #Auxiliar plot
    if sys.argv.count("-f2"):
        column=int (sys.argv[sys.argv.index("-c2")+1])
        fsrc1 = open(sys.argv[sys.argv.index("-f2")+1],'r') #abre el archivo
        hist2=[0]*(bins+1)
        for line in fsrc1:
            words=line.split("\n") 
            words=words[0] 
            words=words.split("\t")
            data=int(words[column])
            if data<bins and data>binsINF:
                hist2[data]+=1
        fsrc1.close()
        
        if sys.argv.count("-line"):
            p2,=plt.plot(range(len(hist2)),hist2,color='b')
        else:
            p2=plt.bar(range(len(hist2)),hist2,color='b')

        if sys.argv.count("-label2"):
            label2=sys.argv[sys.argv.index("-label2")+1]
        else:
            label2=' ';

    #Auxiliar plot
    if sys.argv.count("-f3"):
        column=int (sys.argv[sys.argv.index("-c3")+1])
        fsrc1 = open(sys.argv[sys.argv.index("-f3")+1],'r') #abre el archivo
        hist3=[0]*(bins+1)
        for line in fsrc1:
            words=line.split("\n") 
            words=words[0] 
            words=words.split("\t")
            data=int(words[column])
            if data<bins and data>binsINF:
                hist3[data]+=1
        fsrc1.close()
        
        if sys.argv.count("-line"):
            p3,=plt.plot(range(len(hist3)),hist3,color='g')
        else:
            p3=plt.bar(range(len(hist3)),hist3,color='g')

        if sys.argv.count("-label3"):
            label3=sys.argv[sys.argv.index("-label3")+1]
        else:
            label3=' ';

    #Plot parameters
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if sys.argv.count("-f3"):
        plt.legend([p1, p2, p3], [label1, label2, label3])
    elif sys.argv.count("-f2"):
        plt.legend([p1, p2], [label1, label2])
    else:
        plt.legend([p1], [label1])

    if sys.argv.count("-log"):
        plt.yscale('log')
    grid(True)
    #plt.show()
    fig_name = filename.split('.')
    fig_name =fig_name[0]+'.'+format
    plt.savefig( fig_name, dpi=300,  format=format )
    plt.close() 

else:
    print "Source files not specified"; 
    print "-f file with list"
    print "-c column to plot "
    print "-f2 another file with list"
    print "-c2 column to plot the thrid graph "
    print "-f3 another file with list"
    print "-c3 column to plot the thrid graph "
    print "-log plot in log scale"
    print "-lu maximun bins (default 1000)"
    print "-ld maximun bins (default 0)"
    print "-title Title label"
    print "-xlabel X label"          
    print "-ylabel Y label"          
    print "-line plot stile bar or line (default bar)"       
    print "-label1 legend name for the figure 1"       
    print "-label2 legend name for the figure 2"       
    print "-label3 legend name for the figure 3"       
    print "-format output format (default png)"
    sys.exit()

