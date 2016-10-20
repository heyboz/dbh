# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 12:50:28 2016

@author: Zheyi
"""

import numpy as np
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt

def branch_unique(dbhdf_species):  
    #only the bigest branch for each tag
    dbh_unique = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
    dbh_unique.index = range(len(dbh_unique))
    dbh_unique.index.name = 'No'
    return dbh_unique#dataframe

def trees_distance(dbh_unique):
    #distance matrix
    X = np.array([ dbh_unique['gx'], dbh_unique['gy'] ]).T
    D = scipy.spatial.distance_matrix(X, X)
    return D#array

def cluster_size(dbh_unique, D, dbh_threshold, distance_threshold):
    bigers = dbh_unique[dbh_unique['dbh'] > dbh_threshold]
    clusters = {}
    cluster_size = []
    for i in bigers.index:
        cluster = D[i][D[i] < distance_threshold]
        clusters[i] = cluster
<<<<<<< HEAD
        cluster_size[i] = len(cluster)
   
    n, bins, patches = plt.hist(cluster_size, bins=30, color='steelblue')
    plt.set_xlabel('Cluster Size')
    plt.set_ylabel('Frequency')
    plt.set_title('Cluster Size Distribution_'+sp)
    figname1 = str(dbh_threshold)+'dbh_'+'Cluster_Size_'+sp+'.png'
    plt.savefig('E:/Outputs/Cluster/threshold0.9/ClusterSize/' + figname1)
    plt.close()
    
    n, bins, patches = plt.hist(dbh_unique['dbh'], bins=30, color='darkolivegreen')
    plt.set_xlabel('DBH')
    plt.set_ylabel('Frequency')
    plt.set_title('DBH Distribution_'+sp+'(biggest branches only)')
    figname2 = str(dbh_threshold)+'dbh_'+'DBH_'+sp+'.png'
    plt.savefig('E:/Outputs/Cluster/threshold0.9/DBH/' + figname2)
    plt.close()
    
    sizes = (dbh_unique['dbh']**2) / (4 * 3.14159)
    colors = 1 * (dbh_unique['dbh']>dbh_threshold)
    plt.scatter(dbh_unique['gx'], dbh_unique['gy'], s=sizes, c=colors, alpha=0.7)  
    plt.set_title(sp+'_Plot')
    figname3 = str(dbh_threshold)+'dbh_'+'Location_'+sp+'.png'
    plt.savefig('E:/Outputs/Cluster/threshold0.9/Location/' + figname3)
=======
        cluster_size.append(len(cluster))

    plt.figure(figsize=(4,3))
    n, bins, patches = plt.hist(cluster_size, bins=30)
    plt.xlabel('Cluster Size')
    plt.ylabel('Frequency')
    plt.title('Cluster Size Distribution_'+sp)
    figname1 = sp+'_Cluster_Size'+'_dbh'+str(dbh_threshold)+'_distance'+str(distance_threshold)+'.png'
    plt.savefig('E:/Outputs/Cluster/threshold/ClusterSize/' + figname1)
    plt.close()
    
    plt.figure(figsize=(4,3))    
    n, bins, patches = plt.hist(dbh_unique['dbh'], bins=30)
    plt.xlabel('DBH')
    plt.ylabel('Frequency')
    plt.title('Histogram of DBH_'+sp+'(biggest branches only)')
    plt.savefig('E:/Outputs/Cluster/threshold/DBH/' + sp+'_DBH'+'.png')
    plt.close()
    
    plt.figure(figsize=(16,12))
    sizes = (dbh_unique['dbh']**2) / (4 * 3.14159)
    colors = 1 * (dbh_unique['dbh']>dbh_threshold)
    plt.scatter(dbh_unique['gx'], dbh_unique['gy'], s=sizes, c=colors, alpha=0.7)   
    plt.savefig('E:/Outputs/Cluster/threshold/Location/' + sp+'_Location'+'.png')
>>>>>>> origin/master
    plt.close()
    #return clusters
    
#read data
df = pd.read_csv('HSD_plot_qiaomu.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
#sorted by species abundance
species, counts = np.unique(df1['sp.code'], return_counts=True)
species_c = pd.Series(counts, species)
species_c = species_c.sort_values()
species_sorted = species_c.index

for sp in ['NEOPHA','CRYCON']:
#for sp in species_sorted:
#for sp in sp_maxdbh_sorted:
    dbhdf = df1[df1['sp.code'] == sp]
    dbh_unique = branch_unique(dbhdf)
    D = trees_distance(dbh_unique)
    dbh_threshold = dbhdf['dbh'].quantile(.9)
    #dbh_threshold = 20
    distance_threshold = dbh_threshold * (5/2)
    clusters = cluster_size(dbh_unique, D, dbh_threshold, distance_threshold)
    print sp
