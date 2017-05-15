# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 20:53:46 2017

@author: Administrator
"""

from __future__ import division
import numpy as np
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt
import os
os.chdir('E:/BaiduYunSych/Github/ClusterSize')

#read data
df = pd.read_csv('HSDdata/hsd_plot_qiaomu.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
#sorted by species aboundance
species, counts = np.unique(df1['sp.code'], return_counts=True)
species_c = pd.Series(counts, species)
species_c = species_c.sort_values(ascending = False)
species_sorted = species_c.index

def Biggest_Branch_Only(data):
    dbh_unique = data.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
    dbh_unique.index = range(len(dbh_unique))
    dbh_unique.index.name = 'No' 
    return dbh_unique
                        
output = pd.DataFrame()  
for sp in species_sorted:
    sp = 'DIOMOR'
    dbhdf = df1[df1['sp.code'] == sp]
    dbh_unique0 = Biggest_Branch_Only(dbhdf)
    
    data = dbh_unique0
    dbh_threshold = 10
    
    adults = data[data['dbh'] > dbh_threshold]
    children = data[data['dbh'] <= dbh_threshold]
    print 'adults: ', len(adults), 'children: ', len(children), 'species: ', sp, 'dbh threshold: ', dbh_threshold
    if (len(adults) * len(children)) != 0 and len(adults) > 50 and len(children) > 50:
        A = np.array([ np.array(adults['gx']), np.array(adults['gy']) ]).T
        C = np.array([ np.array(children['gx']), np.array(children['gy']) ]).T
        D = scipy.spatial.distance_matrix(C, A)#axis1: distances from all adults of the focal child, 脚标j
        mom = np.argmin(D, axis=1)
        #find the nearest adult from the focal child and as its mom
        #index is the index of child in 'children', the value is the index of mom in 'adults'
        dist_from_mom = np.min(D,axis=1)#store distance information in a array with the same structure as the mom array
        
        #establish clusters in dictionary
        index_mom, counts = np.unique(mom, return_counts = True)
        Clusters = {}
        Distances = {}
        for j in range(len(children)):
            if j not in index_mom:
                Clusters[j] = 'NA'
                Distances[j] = 'NA'
            Clusters.setdefault(int(mom[j]), []).append(j)
            #key is the index of mom in 'adults', value is the index of child in 'children'
            Distances.setdefault(int(mom[j]), []).append(dist_from_mom[j])
            #key is the index of mom in 'adults', value is the distance between that child from mom
        cs = [0] * len(adults)
        for ia in range(len(adults)): cs[ia] = adults.iat[ia,2]
        
        
        #for ii, im in enumerate(index_mom): cs[im] = counts[ii]
        
        
        
        for im in index_mom:
        ###climax
#            a = Distances[im]
#            if len(a) > 20:
#                plt.hist(a)
#                plt.savefig('Outputs/nearest/trial_analysis/climaxD/diomor_'+ str(im), dpi=100)
#                plt.close()
        ###BasalArea
            BA = 0
            for ic in Clusters[im]:
                BA = BA + children.iat[ic,2]#此处BA为dbh
            cs[im] = cs[im] + BA
        
        
        
        
        
        na = ['NA']*(2518-len(cs))
        cs = cs + na
        output[sp] = cs

output.to_csv('Outputs/nearest/all/cs_scanner/cs_BA0.csv')