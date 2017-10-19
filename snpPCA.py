#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 17:45:09 2017

@author: kingd
"""

from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import visuals as vs
import numpy as np
import matplotlib.cm as cm
from IPython.display import display
from sklearn.mixture import GaussianMixture as GMM
from sklearn.metrics import silhouette_score
from pylab import figure
f1=open('snpvaf1.txt','r')
f2=open('snpvaf2.txt','r')
f3=open('snpvaf3.txt','r')

snpdict={}
snpcounts={}
tpvafs={}
for line in f1:

    line=line.split('\t')
    #print(line)
    chrom=line[0]
    pos=line[1]
    ref=line[2]
    alt=line[3]
    vaf=line[5].replace('\n','')
    ref_to_alt=ref+' -> '+alt
    recurrentvafs=chrom+'_'+pos+'_'+ref+'_'+alt
    if ref_to_alt in snpdict:
        snpcounts[ref_to_alt]+=1
        snpdict[ref_to_alt].append(vaf)
        
    else:
        snpcounts[ref_to_alt]=1        
        snpdict[ref_to_alt]=[vaf]
    tpvafs[recurrentvafs]=[vaf,0,0]

for line2 in f2:
    line2=line2.split('\t')
    #print(line2)
    chrom2=line2[0]
    pos2=line2[1]
    ref2=line2[2]
    alt2=line2[3]
    vaf2=line2[5].replace('\n','')
    recurrentvafs2=chrom2+'_'+pos2+'_'+ref2+'_'+alt2
   # print(recurrentvafs2)
    if recurrentvafs2 in tpvafs:
        l=tpvafs[recurrentvafs2]
        l[1]=vaf2
        tpvafs[recurrentvafs2]=l
       # print(tpvafs[recurrentvafs2])
    else:
        tpvafs[recurrentvafs2]=[0,vaf2,0]
       # print(tpvafs[recurrentvafs2])

for line3 in f3:
    #print(line3)
    line3=line3.split('\t')
   # print(line3)
    chrom3=line3[0]
    pos3=line3[1]
    ref3=line3[2]
    alt3=line3[3]
    vaf3=line3[5].replace('\n','')
    recurrentvafs3=chrom3+'_'+pos3+'_'+ref3+'_'+alt3
    if recurrentvafs3 in tpvafs:
        l3=tpvafs[recurrentvafs3]
        #print("KKKK",l3)
        l3[2]=vaf3
       # print(l3)
        tpvafs[recurrentvafs3]=l3
    else:
        tpvafs[recurrentvafs3]=[0,0,vaf3]

df=pd.DataFrame.from_items(tpvafs.items(), orient='index',columns=['1','2','3'])
display(df.describe())
df = df.astype(float)
display(df.describe())
pd.plotting.scatter_matrix(df, alpha = 0.2, figsize = (14,8), diagonal = 'kde')


pca = PCA(n_components=3) # random_state only available from 0.18.0 onwards
pca.fit(df)

pca_results = vs.pca_results(df, pca)
#vs.biplot(df, reduced_data, pca)
#pca = PCA(n_components=2)
#pca.fit(df)

reduced_data = pca.transform(df)
#reduced_data = pd.DataFrame(reduced_data, columns = ['Dimension 1', 'Dimension 2'])
#components=[2,3,4,5,6,7,8,9,10,11,12,14,20]
components=[9]
All_Scores=[]
for i in range(len(components)):
    clusterer = GMM(n_components=components[i],init_params='random').fit(reduced_data)

    # TODO: Predict the cluster for each data point
    preds = clusterer.predict(reduced_data)

    # TODO: Find the cluster centers
    centers = clusterer.means_

    # TODO: Predict the cluster for each transformed sample data point
   # sample_preds = clusterer.predict(pca_samples)

    # TODO: Calculate the mean silhouette coefficient for the number of clusters chosen
    score = silhouette_score(reduced_data,preds)
    All_Scores.append('n = '+str(components[i]) +": "+str(score))
print(All_Scores)
#df['Cluster']=preds
m = pd.DataFrame(reduced_data, columns = ['Dimension 1', 'Dimension 2','Dimension 3'])
m['cluster']=preds
#vs.cluster_results(rd, preds, centers)

 
fig = plt.figure(1, figsize=(4, 3))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
#m=rand(3,3) # m is an array of (x,y,z) coordinate triplets
cmap = cm.get_cmap('gist_rainbow')
fig = figure()
ax = Axes3D(fig)

	#for i, cluster in plot_data.groupby('Cluster'):   
	 #   cluster.plot(ax = ax, kind = 'scatter', x = 'Dimension 1', y = 'Dimension 2', z='Dimension 3', \
	  #               color = cmap((i)*1.0/(len(centers)-1)), label = 'Cluster %i'%(i), s=30);
for i, cluster in m.groupby('cluster'):
    #plot each point + it's index as text above
    ax.scatter(cluster["Dimension 1"],cluster["Dimension 2"], \
               cluster["Dimension 3"],color = cmap((i)*1.0/(len(centers)-1)),\
               label = 'Cluster %i'%(i), s=30)
 #ax.text(m["Dimension 1"][i],m["Dimension 2"][i],m["Dimension 3"][i],  '%s' % (str(i)), size=20, zorder=1,  
 #color='k') 
for i, c in enumerate(centers):
    
    ax.scatter(c[0], c[1], c[2], color = 'white', edgecolors = 'black', \
	               alpha = 1, linewidth = 2, marker = '.', s=200);
    ax.scatter(c[0], c[1], c[2], marker='$%d$'%(i), alpha = 1, s=100);
ax.set_xlabel('Dimension 1')
ax.set_ylabel('Dimension 2')
ax.set_zlabel('Dimension 3')
ax.set_title("Cluster Learning on PCA Data - Centroids Marked by Number\nTransformed Sample Data Marked by Black Circle");
plt.show()
#for name, label in [('Setosa', 0), ('Versicolour', 1), ('Virginica', 2)]:

#df.drop('VAR', 1,inplace=True)
"""
sv = sorted(snpcounts)
#svaf=pd.DataFrame(snpcounts)       
import matplotlib.pyplot as plt

#D = {u'Label1':26, u'Label2': 17, u'Label3':30}

plt.bar(range(len(snpcounts)), snpcounts.values(), align='center')
plt.xticks(range(len(snpcounts)), snpcounts.keys(),rotation=90)

plt.show()
"""
df['cluster']=preds
clus0=df[df['cluster']==0]
clus1=df[df['cluster']==1]
clus2=df[df['cluster']==2]
clus3=df[df['cluster']==3]
clus4=df[df['cluster']==4]
clus5=df[df['cluster']==5]
clus6=df[df['cluster']==6]
clus7=df[df['cluster']==7]
clus8=df[df['cluster']==8]




#clus0['info'] = clus0.index
clus0.reset_index(level=0, inplace=True)
clus1.reset_index(level=0, inplace=True)
clus2.reset_index(level=0, inplace=True)
clus3.reset_index(level=0, inplace=True)
clus4.reset_index(level=0, inplace=True)
clus5.reset_index(level=0, inplace=True)
clus6.reset_index(level=0, inplace=True)
clus7.reset_index(level=0, inplace=True)
clus8.reset_index(level=0, inplace=True)

clus0.to_csv('clus0.csv',sep='\t',index_label='Remove')
clus1.to_csv('clus1.csv',sep='\t',index_label='Remove')
clus2.to_csv('clus2.csv',sep='\t',index_label='Remove')
clus3.to_csv('clus3.csv',sep='\t',index_label='Remove')
clus4.to_csv('clus4.csv',sep='\t',index_label='Remove')
clus5.to_csv('clus5.csv',sep='\t',index_label='Remove')
clus6.to_csv('clus6.csv',sep='\t',index_label='Remove')
clus7.to_csv('clus7.csv',sep='\t',index_label='Remove')
clus8.to_csv('clus8.csv',sep='\t',index_label='Remove')

