#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os,argparse,math
from subprocess import Popen, PIPE
from pylab import *



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot histogram of total edges')

    parser.add_argument('-f', '--file',     help='Links files, comma separeted for several files', required=True, type = str)
    parser.add_argument('-c', '--columns',  help='Columns to plot, comma separeted for several fieles', required=True, type = str)
    parser.add_argument('-e', '--edge',     help='Edge-threshold filters', default='0,10000',type = str)
    parser.add_argument('-t', '--title',    help='Title plot name', default='', type = str)
    parser.add_argument('-a', '--axes',     help='Axes labels', default='Total edges,Number of Reads', type = str)
    parser.add_argument('-l', '--labels',   help='Legend labels', type = str)


    args = parser.parse_args()
    
    if args.file and args.columns:
        files=args.file.split(',')
        columns=args.columns.split(',')
        edges=args.edge.split(',')
        bins=int(edges[1])
        binsINF=int(edges[0])
        
        title=args.title
        axes_names=args.axes.split(',')
        xlabel=axes_names[0]
        ylabel=axes_names[1]
        colors=['r','g','b','m','c','y','gray','orange']
        
        #Main plot
        i=0
        p=[]
        for fname in files:
            hist=[0]*(bins+1)
            filename=fname
            fsrc1 = open(filename,'r') 
            for line in fsrc1:
                words=line.split('\n') 
                words=words[0] 
                words=words.split('\t')
                data=int(words[int(columns[i])])
                if data<bins and data>binsINF:
                    hist[data]+=1
            fsrc1.close()
            color=colors[i%(len(colors))]
            if args.labels:
                labels=args.labels.split(',')
                plt.plot(range(len(hist)),hist,color=color,label=labels[i])
            else:
                plt.plot(range(len(hist)),hist,color=color)    
            i+=1
        
        #Plot parameters
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()
        grid(True)
        plt.show()
