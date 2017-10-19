#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 18:43:30 2017

@author: kingd
"""
import pandas as pd
import numpy as np
f=pd.read_csv('clus8.csv',sep='\t')

df_genes= pd.DataFrame(columns=['Hugo_Symbol','Chromosome','Start_Position','Ref','Alt'])
t1=pd.read_table('MMRF_1049_1.txt',sep='\t')
t2=pd.read_table('MMRF_1049_2.txt',sep='\t')
t3=pd.read_table('MMRF_1049_3.txt',sep='\t')
#clus=open('Clus6Genes.txt','w')
#genes=pd.DataFrame(columns='Genes')\
test=open('GenesClus8.txt','w')
glist=[]
info=f['index']
f=f.drop('Remove',1)
t4=t1.merge(t2,how='outer')
tall=t4.merge(t3,how='outer')
for line in info:


    nline=line.split('_')
    #print(line)
    chrom=nline[0]
    #print(line)
    chrom=nline[0];ref=nline[3]
    pos=nline[1];alt=nline[3]

    for i, d in tall.groupby('Chromosome'):
        if chrom == i:
            b=np.int64(pos)
            a=d['Start_Position'][d['Start_Position']==b]
            if len(a)==1:
                gene=d['Hugo_Symbol'][d['Start_Position']==b].item()
            
                vaf1=f['1'][f['index']==line].item()
                vaf2=f['2'][f['index']==line].item()
                vaf3=f['3'][f['index']==line].item()
                vaf1=str(vaf1)
                vaf2=str(vaf2)
                vaf3=str(vaf3)
                wline= gene+'\t'+chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+vaf1+'\t' \
                +vaf2+'\t'+vaf3+'\n'
            
                test.write(wline)
                #print(wline)            
           # print('finall######FGHFHGFUFGy',d[d['Start_Position']==b])
            if len(a)>1:
                for i in a:
                    print('KKKK',a)
                    print('KKKK',i)
                    print(d['Hugo_Symbol'][d['Start_Position']==b])
         

"""                
    for i, d in t2.groupby('Chromosome'):
        if chrom == i:
            b=np.int64(pos)
            a=d['Start_Position'][d['Start_Position']==b]
            if len(a)>=1:
                print(nline)
                print('finall######FGHFHGFUFGy',d[d['Start_Position']==b])
                gene=d['Hugo_Symbol'][d['Start_Position']==b].item()
                vaf1=f['1'][f['index']==line].item()
                vaf2=f['2'][f['index']==line].item()
                vaf3=f['3'][f['index']==line].item()
                vaf1=str(vaf1)
                vaf2=str(vaf2)
                vaf3=str(vaf3)
                wline= gene+'\t'+chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+vaf1+'\t' \
                +vaf2+'\t'+vaf3+'\n'
                
                test.write(wline)
                print(wline)

    for i, d in t3.groupby('Chromosome'):
        if chrom == i:
            b=np.int64(pos)
            a=d['Start_Position'][d['Start_Position']==b]
            if len(a)>=1:
                print(nline)
                print('finall######FGHFHGFUFGy',d[d['Start_Position']==b])
                gene=d['Hugo_Symbol'][d['Start_Position']==b].item()
                vaf1=f['1'][f['index']==line].item()
                vaf2=f['2'][f['index']==line].item()
                vaf3=f['3'][f['index']==line].item()
                vaf1=str(vaf1)
                vaf2=str(vaf2)
                vaf3=str(vaf3)
                wline= gene+'\t'+chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+vaf1+'\t' \
                +vaf2+'\t'+vaf3+'\n'
                
                test.write(wline)
                print(wline) 
                """
test.close()
            