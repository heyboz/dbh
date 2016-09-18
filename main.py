# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 12:50:28 2016

@author: Zheyi
"""

import numpy as np
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg') 
#import networkx as nx

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
        cluster_size.append(len(cluster))

    n, bins, patches = plt.hist(cluster_size, bins=30)
    plt.xlabel('Cluster Size')
    plt.ylabel('Frequency')
    plt.title('Cluster Size Distribution_'+sp)
    figname1 = 'Cluster_Size_'+sp+'_dbh'+str(dbh_threshold)+'_distance'+str(distance_threshold)+'.png'
    plt.savefig('ClusterSize/qiaomu/' + figname1)
    plt.clf()
    
    n, bins, patches = plt.hist(dbh_unique['dbh'], bins=30)
    plt.xlabel('DBH')
    plt.ylabel('Frequency')
    plt.title('Histogram of DBH_'+sp+'(biggest branches only)')
    figname2 = 'DBH_'+sp+'.png'
    plt.savefig('DBH/qiaomu/' + figname2)
    plt.clf()
    return clusters
    
#read data
df = pd.read_csv('HSD_plot_qiaomu.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
species = np.unique(df1['sp.code'])
for sp in species:
    dbhdf = df1[df1['sp.code'] == sp]
    dbh_unique = branch_unique(dbhdf)
    D = trees_distance(dbh_unique)
    clusters = cluster_size(dbh_unique, D, 20, 50)
    print sp