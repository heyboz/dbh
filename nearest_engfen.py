# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 18:00:34 2016

@author: Administrator
"""
from __future__ import division
import numpy as np
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt
import random

#read data
df = pd.read_csv('HSD_plot_qiaomu.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
#sorted by species aboundance
species, counts = np.unique(df1['sp.code'], return_counts=True)
species_c = pd.Series(counts, species)
species_c = species_c.sort_values(ascending = False)
species_sorted = species_c.index

def nullmodel(rawdata, dbh_adultRatio, dbh_distribution=False, spacial=False):
    if dbh_adultRatio != 0:
        a = int(len(rawdata)*dbh_adultRatio)
        dbhr = [1]*len(rawdata)    
        for i in range(a):
            dbhr[i] = random.uniform(10, max(rawdata['dbh']))
        for i in range(a+1, len(rawdata)):
            dbhr[i] = random.uniform(0, 10)
        rawdata['dbh'] = dbhr
    if dbh_distribution == True:
        dbhd = rawdata['dbh'][:]
        random.shuffle(dbhd)
        rawdata['dbh'] = dbhd
    if spacial == True:
        xr = [1]*len(rawdata)
        yr = [0]*len(rawdata)
        for i in range(len(rawdata)):
            xr[i] = random.uniform(0, 1000)
            yr[i] = random.uniform(0, 500)
        rawdata['gx'] = xr
        rawdata['gy'] = yr
    return rawdata

def Cluster_Size(data):
    dbhrange = range(int(min(data['dbh'])), int(max(data['dbh'])), 1)
    norm_mean_cluster_size = [1] * len(dbhrange)
    #norm_mean_d = [1] * len(dbhrange)     
    for i in range(len(dbhrange)):
        adults = data[data['dbh'] > dbhrange[i]]
        children = data[data['dbh'] <= dbhrange[i]]
        print 'adults: ', len(adults), 'children: ', len(children), 'species: ', dbhrange[i]
        
        A = np.array([ np.array(adults['gx']), np.array(adults['gy']) ]).T
        C = np.array([ np.array(children['gx']), np.array(children['gy']) ]).T
        D = scipy.spatial.distance_matrix(C, A)#axis1: distances from all adults of the focal child
        mom = np.argmin(D, axis=1)
        #find the nearest adult from the focal child and as its mom
        #index is the index of child in 'children', the value is the index of mom in 'adults'
        dist_from_mom = np.min(D,axis=1)#store distance information in a array with the same structure as the mom array
        
        #establish clusters in dictionary
        index_mom, counts = np.unique(mom, return_counts = True)
        Clusters = {}
        Distances = {}
        for j in range(len(children)):
            Clusters.setdefault(int(mom[j]), []).append(j)
            #key is the index of mom in 'adults', value is the index of child in 'children'
            Distances.setdefault(int(mom[j]), []).append(dist_from_mom[j])
            #key is the index of mom in 'adults', value is the distance between that child from mom
        
        cluster_size = [len(Clusters[index_mom[k]]) for k in range(len(index_mom))]#index is the index of mom in 'adults'
        norm_mean_cluster_size[i] = sum([x*x for x in cluster_size])/(len(dbh_unique)*len(dbh_unique))
        #cluster_size_list = [[len(Clusters[index_mom[m]])]*len(Clusters[index_mom[m]]) for m in range(len(index_mom))]
        #c = [y for x in cluster_size_list for y in x]
        #dist = [Distances[index_mom[n]] for n in range(len(index_mom))]    
        #d = [y for x in dist for y in x]
        #norm_mean_d[i] = sum(d)/((max(d)-min(d))*len(d))

    return cluster_size, norm_mean_cluster_size    
    
for sp in species_sorted[:20]:
    dbhdf = df1[df1['sp.code'] == sp]#ENGFEN ALTCHI NEOPHA
        
    #biggest branch only
    dbh_unique0 = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
    dbh_unique0.index = range(len(dbh_unique0))
    dbh_unique0.index.name = 'No' 
        
    dbh_unique = nullmodel(rawdata=dbh_unique0, dbh_adultRatio=0, dbh_distribution=True, spacial=True)
    cs_data, nmcs_data = Cluster_Size(dbh_unique0)
        

    plt.plot(dbhrange, norm_mean_cluster_size)
    plt.xlabel('dbh')
    plt.ylabel('normalized mean cluster size')
    plt.savefig('E:/Outputs/Cluster/nearest/dbh-percolation' + sp)
    plt.close()
    