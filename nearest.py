# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 18:00:34 2016

@author: Administrator
"""

import numpy as np
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt

#read data
df = pd.read_csv('HSD_plot_qiaomu.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
#sorted by species aboundance
species, counts = np.unique(df1['sp.code'], return_counts=True)
species_c = pd.Series(counts, species)
species_c = species_c.sort_values()
species_sorted = species_c.index
#sorted by maximal dbh
#sp_maxdbh = df1.groupby('sp.code', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
#sp_maxdbh = sp_maxdbh.drop_duplicates(subset="sp.code")
#sp_maxdbh_sorted = sp_maxdbh.sort(columns='dbh', axis=0)['sp.code']

for sp in species_sorted:
    dbhdf = df1[df1['sp.code'] == sp]
    #biggest branch only
    dbh_unique = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
    dbh_unique.index = range(len(dbh_unique))
    dbh_unique.index.name = 'No'    
   
    #find nearest mom
    dbh_threshold = 10
    adults = dbh_unique[dbh_unique['dbh'] > dbh_threshold]
    children = dbh_unique[dbh_unique['dbh'] <= dbh_threshold]
    print 'adults: ', len(adults), 'children: ', len(children), 'species: ', sp
    if len(adults) == 0 or len(children) == 0:
        continue
    A = np.array([ np.array(adults['gx']), np.array(adults['gy']) ]).T
    C = np.array([ np.array(children['gx']), np.array(children['gy']) ]).T
    D = scipy.spatial.distance_matrix(C, A)#axis1: distances from all adults of the focal child
    
    mom = np.argmin(D, axis=1)#find the nearest adult from the focal child and as its mom 
    dist_from_mom = np.min(D,axis=1)#store distance information in a array with the same structure as the mom array
    
    #establish clusters in dictionary
    index_mom, cluster_size = np.unique(mom, return_counts = True)
    Clusters = {}
    Distances = {}
    for i in range(len(children)):
        Clusters.setdefault(int(mom[i]), []).append(i)
        Distances.setdefault(int(mom[i]), []).append(dist_from_mom[i])
        
    ###calculate cluster size distribution
    ###visiulize cluster structure
    ###point out mothers and children in HSD map if it possible
    
    #show cluster size distribution
    Fig = plt.figure(figsize = (8,13.5))
    Ax1 = Fig.add_subplot(311)
    Ax2 = Fig.add_subplot(312)
    Ax3 = Fig.add_subplot(313)
    
    n, bins, patches = Ax1.hist(cluster_size, bins=50)
    Ax1.set_xlabel('Cluster Size')
    Ax1.set_ylabel('Frequency')
    Ax1.set_title('Cluster Size Distribution_'+sp)

    n, bins, patches = Ax2.hist(dbh_unique['dbh'], bins=50)
    Ax2.set_xlabel('DBH')
    Ax2.set_ylabel('Frequency')
    Ax2.set_title('Histogram of DBH_'+sp+'(biggest branches only)')

    sizes = (dbh_unique['dbh']**2) / (4 * 3.14159)
    colors = 1 * (dbh_unique['dbh']>dbh_threshold)
    Ax3.scatter(dbh_unique['gx'], dbh_unique['gy'], s=sizes, c=colors, alpha=0.7)
    Ax3.set_title('Trees Map')

    figname1 = sp+'_Cluster_Size'+'_dbh'+str(dbh_threshold)+'.png'
    Fig.savefig('E:/Outputs/Cluster/nearest/' + figname1)
    plt.close()