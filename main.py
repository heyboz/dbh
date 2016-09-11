# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 12:50:28 2016

@author: Administrator
"""

import numpy as np
import scipy as sp
import pandas as pd
import scipy.spatial

df = pd.read_csv('HSD_plot.csv')
df1 = df[['sp.code', 'tag', 'dbh', 'gx', 'gy']]
dbh_ADF = df1[df1['sp.code'] == 'ALTCHI']
dbhdf = dbh_ADF[:10]

dbh_unique = dbhdf.groupby('tag', sort = False).apply(lambda t: t[t.dbh==t.dbh.max()])

X = np.array([ dbh_unique['gx'], dbh_unique['gy'] ]).T
D = sp.spatial.distance_matrix(X, X)


