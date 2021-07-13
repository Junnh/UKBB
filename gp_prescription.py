#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import time
import datetime
from datetime import datetime

# read gp prescription data
ps = pd.read_csv('gp_scripts.txt', sep='\t', encoding='latin_1')

# drop null, error issue_date
null_time_index = ps[ps['issue_date'].isnull()].index
error_time_index = ps[ps['issue_date'].isin(['07/07/2037', '01/01/1901', '02/02/1902', '03/03/1903'])].index
ps.drop(np.union1d(null_time_index, error_time_index), axis=0, inplace=True)

# read Read2 drug code data
r2_drug = pd.read_csv('read2_drug.csv')
r2_drug.drop([67612, 67613], axis=0, inplace=True)
r2_drug.rename(columns={'read_code':'read_2', 'term_description':'drug_term_description'}, inplace=True)
r2_drug.drop('status_flag', axis=1, inplace=True)
r2_drug.head()

# join prescription-read2 drug
ps_read = pd.merge(ps, r2_drug, on='read_2', how='left')

# drop null read2 code, sort by time
null_drug = ps_read[ps_read['drug_name'].isnull()]
no_info_index = null_drug[null_drug['read_2'].isnull()].index
ps_read.drop(no_info_index, axis=0, inplace=True)
ps_read.rename(columns={'issue_date':'event_dt'}, inplace=True)
ps_read['event_dt'] = pd.to_datetime(ps_read['event_dt'])
ps_read.sort_values(by='event_dt', inplace=True)

# ### make DB
# con = sqlite3.connect("/path/gp.db")
# sqlite3.Connection
# cursor = con.cursor()
# read_gp_processed.to_sql('gp_prescription', con) # table name = gp_prescription

