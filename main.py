# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 12:50:28 2016

@author: Administrator
"""

import numpy as np
import scipy as sp
import pandas as pd
import scipy.spatial
import matplotlib.pyplot as plt

#read data
df = pd.read_csv('HSD_plot.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
dbh_ADF = df1[df1['sp.code'] == 'MYTLAO']
dbhdf = dbh_ADF

#only the bigest branch for each tag
dbh_unique = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])

#distance matrix
X = np.array([ dbh_unique['gx'], dbh_unique['gy'] ]).T
D = sp.spatial.distance_matrix(X, X)

#define big trees according to dbh
n, bins, patches = plt.hist(dbh_unique['dbh'], bins=30)
plt.xlabel('DBH')
plt.ylabel('frequency')
plt.title('Histogram of DBH(the bigest branches only)')
plt.show()

n, bins, patches = plt.hist(dbh_ADF['dbh'], bins=30)
plt.xlabel('DBH')
plt.ylabel('frequency')
plt.title('Histogram of DBH')
plt.show()

