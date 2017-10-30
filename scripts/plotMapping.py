#!/usr/bin/env python
import os,argparse,math
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from subprocess import Popen, PIPE

def add_depth(namelist1,positionNorm,depth):
        if positionNorm in namelist1:
                namelist1[positionNorm]+=depth
        else:
                namelist1[positionNorm]=depth
def plotNorm(args):    
        yvalues=[0,100]
        lim=100
        norma=1000   
        labelA = args.label1
        fileA = open(args.file1,'r') 
        namelist1= dict();
        oldID="";
        position=0;
        i=1
        firstTime=True
        firstPlot=False
        for line in fileA: #File with original Crypto depth
                word=line.split("\n") 
                word=word[0].split()
                Id=word[0]
                if (firstTime):
                        firstTime=False
                        firstPlot=True
                        oldID=Id
                elif(Id!=oldID):
                        oldID=Id
                        plotFig(namelist1,i,norma,labelA,'r',firstPlot)
                        plt.yticks(yvalues)
                        i+=1
                        grid(True)
                        firstPlot=False
                position=int(word[1])
                lenght=position
                positionNorm=lenght/norma
                depth=int(word[2])
                add_depth(namelist1,positionNorm,depth)

        plotFig(namelist1,i,norma,labelA,'r',False)
        plt.yticks(yvalues)
        fileA.close()
        
        fileA = open(args.file2,'r') #File with original binning depth
        labelA = args.label2
        namelist1= dict();
        oldID="";
        position=0;
        i=1
        firstTime=True
        for line in fileA:
                word=line.split("\n") 
                word=word[0].split()
                Id=word[0]
                if (firstTime):
                        firstTime=False
                        oldID=Id
                elif(Id!=oldID):
                        oldID=Id
                        plotFig(namelist1,i,norma,labelA,'b', firstPlot)
                        plt.yticks(yvalues)
                        i+=1

                position=int(word[1])
                lenght=position
                positionNorm=lenght/norma
                depth=int(word[2])
                add_depth(namelist1,positionNorm,depth)

        plotFig(namelist1,i,norma,labelA,'b', False)
        plt.yticks(yvalues)
        fileA.close()
        
        grid(True)
        plt.title(args.title)
        plt.xlabel("Kbp")
        plt.subplots_adjust(hspace=1,wspace=0.35)
        
      
        fig_name = args.output
        plt.savefig( fig_name, dpi=300, format='pdf')
        return

def plotFig(namelist1,i, norma,labelA, col, firstPlot):  
        depthArray =[]; base =[];  
        plt.subplot(8, 1,1*i) 
        for item in namelist1:
                base.append(item)
                depthArray.append(namelist1[item]/norma)
                fig1=plt.plot(base,depthArray, color=col,linewidth=1.0)
        h = plt.ylabel('CH_'+str(i))
        if firstPlot:
            ax = plt.gca()
            fakeLine1 = plt.Line2D([0,0],[0,1], color='red', linestyle='-')
            fakeLine2 = plt.Line2D([0,0],[0,1], color='blue', linestyle='-')
            ax.legend([fakeLine1,fakeLine2], ["All reads C. hominis", "Reads from bins"], loc=4, ncol=2, bbox_to_anchor=(1.0,  1.0), fontsize='small')
        

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='bam plotting depth')
  parser.add_argument('-f1', '--file1', help='File name with depth', type = str)
  parser.add_argument('-f2', '--file2', help='File name with depth', type = str)
  parser.add_argument('-lb1', '--label1',  help='label 1',  type = str, default=" ")
  parser.add_argument('-lb2', '--label2',  help='label 2',  type = str, default=" ")                       
  parser.add_argument('-t', '--title',  help='title',  type = str, default=" ")                       
  parser.add_argument('-o', '--output',  help='save name',  type = str, default=" ")                       
  args = parser.parse_args()


  if args.file1:
    plotNorm(args)

