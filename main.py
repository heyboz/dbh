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
from matplotlib.mlab import griddata
import os
os.chdir('E:/BaiduYunSych/Github/ClusterSize')

#read data
df = pd.read_csv('HSDdata/hsd_plot_qiaomu.csv')
#dioecious = pd.read_csv('HSD_dioecious.csv')
envir = pd.read_csv('HSDdata/hsd_envir20_data.csv')
envir['xy'] = envir['x'] * 1000 + envir['y']
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

def Nullmodel(rawdata, dbh_adultRatio, dbh_distribution=False, spatial=False):
    #when dbh_distribution true, spatial false, keep spatial but randomize dbh
    #when dbh_distribution false, spatial true, randomize everything
    newdata = pd.DataFrame.copy(rawdata)    
    if dbh_adultRatio != 0:
        a = int(len(newdata)*dbh_adultRatio)
        dbhr = [1]*len(newdata)    
        for i in range(a):
            dbhr[i] = random.uniform(10, max(newdata['dbh']))
        for i in range(a+1, len(newdata)):
            dbhr[i] = random.uniform(0, 10)
        newdata['dbh'] = dbhr
    if dbh_distribution == True:
        dbhd = newdata['dbh'][:]
        random.shuffle(dbhd)
        newdata['dbh'] = dbhd
    if spatial == True:
        xr = [1]*len(newdata)
        yr = [0]*len(newdata)
        for i in range(len(newdata)):
            xr[i] = random.uniform(0, 1000)
            yr[i] = random.uniform(0, 500)
        newdata['gx'] = xr
        newdata['gy'] = yr
    return newdata

def Cluster_Size(data, dbh_threshold):
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
    cluster_size = [0] * len(adults)
    for ii, im in enumerate(index_mom): cluster_size[im] = counts[ii]
    #cluster_size = [len(Clusters[index_mom[k]]) for k in range(len(index_mom))]#index is the index of mom in 'adults'
    #cluster_size_list = [[len(Clusters[index_mom[m]])]*len(Clusters[index_mom[m]]) for m in range(len(index_mom))]
    #c = [y for x in cluster_size_list for y in x]
    #dist = [Distances[index_mom[n]] for n in range(len(index_mom))]    
    #d = [y for x in dist for y in x]
    #norm_mean_d[i] = sum(d)/((max(d)-min(d))*len(d))    
    
    return cluster_size

def DBH_Percolation(data):
    dbhrange = range(int(min(data['dbh']))-1, int(max(data['dbh']))+1, 1)
    norm_mean_cluster_size = [1] * len(dbhrange)
    #norm_mean_d = [1] * len(dbhrange)     
    for i in range(len(dbhrange)):
        cluster_size= Cluster_Size(data, dbhrange[i])
        norm_mean_cluster_size[i] = sum([x*x for x in cluster_size])/(len(data)*len(data))
    return norm_mean_cluster_size, dbhrange

def NC(data, dbh_threshold, mode='intra_adult', datall=False):
    adults = data[data['dbh'] > dbh_threshold]
    A = np.array([ np.array(adults['gx']), np.array(adults['gy']) ]).T
    #T = np.array([ np.array(data['gx']), np.array(data['gy']) ]).T
    nc = [0] * len(adults)
    if mode == 'intra_adult':
        Dnc = scipy.spatial.distance_matrix(A, A)
        for m in range(len(adults)):
            for n in range(len(adults)):
                if Dnc[m,n] < 20 and Dnc[m,n] > 0:
                    dbh = data.iat[n, 2]
                    base_area = (dbh/2)*(dbh/2)*3.1415926
                    nc[m] = nc[m] + (base_area/Dnc[m,n])
    elif mode == 'intra_child':
        children = data[data['dbh'] <= dbh_threshold]
        C = np.array([ np.array(children['gx']), np.array(children['gy']) ]).T
        Dnc = scipy.spatial.distance_matrix(A, C)
        for m in range(len(adults)):
            for n in range(len(children)):
                if Dnc[m,n] < 20 and Dnc[m,n] > 0:
                    dbh = data.iat[n, 2]
                    base_area = (dbh/2)*(dbh/2)*3.1415926
                    nc[m] = nc[m] + (base_area/Dnc[m,n])
    elif mode == 'inter':
        dbhall = Biggest_Branch_Only(df1)
        dbhinter = dbhall[dbhall['sp.code'] != sp]
        Tt = np.array([ np.array(dbhinter['gx']), np.array(dbhinter['gy']) ]).T
        Dnct = scipy.spatial.distance_matrix(A, Tt)
        for m in range(len(adults)):
            for n in range(len(dbhinter)):
                if Dnct[m,n] < 20 and Dnct[m,n] > 0:
                    dbh = dbhinter.iat[n, 2]
                    base_area = (dbh/2)*(dbh/2)*3.1415926
                    nc[m] = nc[m] + (base_area/Dnct[m,n])
    else: print 'invalid mode'
    return nc

def Quadrat_Center(data, dbh_threshold):
    xy = []
    adults = data[data['dbh'] > dbh_threshold]
    for i, gx in enumerate(adults['gx']):
        qx = int(gx/20)*20+10
        if qx > 990: qx = 990
        qy = int(adults.iat[i,4]/20)*20+10
        if qy > 490: qy = 490
        xy.append(qx * 1000 + qy)
    return xy
    
#maxdbhlist = []
#for sp in species_sorted:
#    dbhdf = df1[df1['sp.code'] == sp]
#        
#    #biggest branch only
#    dbh_unique0 = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])
#    dbh_unique0.index = range(len(dbh_unique0))
#    dbh_unique0.index.name = 'No' 
#    
#    maxdbh = max(dbh_unique0['dbh'])
#    text = "species: %s, maxdbh: %s \n" % (sp, maxdbh)
#    maxdbhlist.append(maxdbh)
#    f = open('E:/Outputs/Cluster/nearest/maxdbh.txt', 'a')
#    f.writelines(text)
#f.close()
#plt.figure(figsize=(12,9))
#plt.hist(maxdbhlist,bins=20,color='steelblue')
#plt.xlabel('maximum dbh')
#plt.ylabel('frequency')
#plt.savefig('E:/Outputs/Cluster/nearest/dbh-percolation/maxdbh', dpi=100)

species_demo = ['MACBRE','QUELIT','LITACU',
                'DISRAC','ARTSTY','DIOMOR',
                'SYMLAU','ORMPAC','CINPAR']
for sp in species_sorted:
#for sp in species_demo:
    #sp = 'IXORET'
    dbhdf = df1[df1['sp.code'] == sp]
    
    dbh_unique0 = Biggest_Branch_Only(dbhdf)
    
    dbh_threshold = 10
    adults = dbh_unique0[dbh_unique0['dbh'] > dbh_threshold]
    children = dbh_unique0[dbh_unique0['dbh'] <= dbh_threshold]
    print 'adults: ', len(adults), 'children: ', len(children), 'species: ', sp, 'dbh threshold: ', dbh_threshold
    
    if (len(adults) * len(children)) != 0 and len(adults) > 50 and len(children) > 50:
    
    #dbh_unique1 = Nullmodel(rawdata=dbh_unique0, dbh_adultRatio=0, dbh_distribution=False, spatial=True)
    #dbh_unique2 = Nullmodel(rawdata=dbh_unique0, dbh_adultRatio=0, dbh_distribution=True, spatial=False)
    
        #cs0 = Cluster_Size(dbh_unique0, 10)
    #cs1 = Cluster_Size(dbh_unique1, 10)
    #cs2 = Cluster_Size(dbh_unique2, 10)
    
#####calculating nc intraspecies   
        nc0a = NC(dbh_unique0, 10, mode = 'intra_adult')
        nc0c = NC(dbh_unique0, 10, mode = 'intra_child')
    #nc1 = NC(dbh_unique1, 10)
    #nc2 = NC(dbh_unique2, 10)
    
        nc0t = NC(dbh_unique0, 10, mode='inter')
    #nc1t = NC(dbh_unique1, 10, mode='inter')
    #nc2t = NC(dbh_unique2, 10, mode='inter')
    
######dbh percolation
#    nmcs_data, dbhrange = DBH_Percolation(dbh_unique0)
#    nmcs_spatial, dbhrange = DBH_Percolation(dbh_unique1)
#    nmcs_dbh, dbhrange = DBH_Percolation(dbh_unique2)    
#    
#    fig, ax1 = plt.subplots()
#    fig.set_size_inches(16, 9)
#    ax2 = ax1.twinx()
#    ax2.hist(dbh_unique0['dbh'], bins=round(max(dbh_unique0['dbh'])), histtype='stepfilled', facecolor='silver', alpha=0.5)
#    ax1.plot(dbhrange, nmcs_data, '-', label='empirical', linewidth=3)
#    ax1.plot(dbhrange, nmcs_dbh, ':', label='without parental effect', linewidth=3)
#    ax1.plot(dbhrange, nmcs_spatial, '--', label='without spatial aggregation', linewidth=3)
#    ax1.legend(loc='higher left')
#    ax1.set_xlabel('dbh')
#    ax1.set_ylabel('normalized mean cluster size')
#    ax2.set_ylabel('frequency')
#    plt.tight_layout(pad=8, w_pad=2, h_pad=2)
#    plt.savefig('E:/Outputs/Cluster/nearest/dbh-percolation/' + sp, dpi=100)
#    plt.close()

#####contour
#        contour = pd.read_csv('hsd_envir20_data.csv')
#        x = contour['x']
#        y = contour['y']
#        z = contour['elev']
#        xi = np.linspace(10,990,50)
#        yi = np.linspace(10,490,25)
#        zi = griddata(x,y,z,xi,yi,interp='linear')
#        
#        plt.figure(figsize=(20,10))
#        CS = plt.contour(xi,yi,zi,linewidths=2,colors='k')
#        plt.clabel(CS, inline=1, fontsize=25)
#        sizes = (dbh_unique0['dbh'])**1.5
#        #colors = 1 * (dbh_unique['dbh']<=dbh_threshold)
#        colors = 1 * (dbh_unique0['dbh']<=10)
#        plt.scatter(dbh_unique0['gx'], dbh_unique0['gy'], s=sizes, c=colors, alpha=0.9, edgecolors='none', cmap='Spectral')
#        
#        plt.xlim(0,1000)
#        plt.ylim(0,500)
#        plt.savefig('E:/Outputs/Cluster/nearest/contour/all/' + sp, dpi=100)

####qudart center
        xy0 = Quadrat_Center(dbh_unique0, 10)
    #x1, y1 = Quadrat_Center(dbh_unique1, 10)
    #x2, y2 = Quadrat_Center(dbh_unique2, 10)

####output data
        for mi in [0]:#[0,1,2]
            exec'output = dbh_unique%s[dbh_unique%s[\'dbh\'] > 10]' % (mi, mi)
            #exec'output[\'cluster.size\'] = cs%s' % mi
            exec'output[\'nc_intra_adult\'] = nc%sa' % mi
            exec'output[\'nc_intra_child\'] = nc%sc' % mi
            exec'output[\'nc_inter\'] = nc%st' % mi
            exec'output[\'xy\'] = xy%s' % mi
            output1 = pd.merge(output, envir, on = 'xy')
            output1.to_csv('Outputs/nearest/all/nc/' + sp + '.csv')