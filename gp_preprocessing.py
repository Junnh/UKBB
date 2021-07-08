#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import time
import datetime
from datetime import datetime
import sqlite3


### read files
gp = pd.read_csv('gp_clinical.txt', sep='\t', encoding='cp949')

r3_i10 = pd.read_csv('r3_icd10.csv') # Read3-ICD10 mapping
r3_i10 = r3_i10.drop([116374, 116375, 116376], axis=0)
r3_i10 = r3_i10.rename(columns={'read_code':'read_3'})

r3_i9 = pd.read_csv('r3_icd9.csv') # Read3-ICD9 mapping
r3_i9 = r3_i9.drop([67152, 67153, 67154], axis=0)
r3_i9 = r3_i9.rename(columns={'read_code':'read_3'})

r2_i10 = pd.read_csv('r2_icd10.csv') # Read2-ICD10 mapping
r2_i10 = r2_i10.drop([36664, 36665, 36666], axis=0)
r2_i10 = r2_i10.rename(columns={'read_code':'read_2'})

r2_i9 = pd.read_csv('r2_icd9.csv') # Read2-ICD9 mapping
r2_i9 = r2_i9.drop([35661, 35662, 35663], axis=0)
r2_i9 = r2_i9.rename(columns={'read_code':'read_2'})

r3 = pd.read_csv('read3.csv') # Read3 description
r3 = r3.drop([0, 1, 2, 395652, 395653], axis=0)
r3 = r3.rename(columns={'read_code':'read_3'})

r2 = pd.read_csv('read2.csv') # Read2 descripiton
r2 = r2.drop([101881, 101882], axis=0)
r2 = r2.rename(columns={'read_code':'read_2'})


### save index before preprocessing
gp_r3 = gp['read_3'].dropna()
gp_r3 = gp.loc[gp_r3.index]
gp_r3_ids = gp_r3['eid'].values

gp_r2 = gp['read_2'].dropna()
gp_r2 = gp.loc[gp_r2.index]
gp_r2_ids = gp_r2['eid'].values

# inter = np.intersect1d(gp_r3_ids, gp_r2_ids)
gp_r3_index = gp_r3.index
gp_r2_index = gp_r2.index


### Join read code, gp file
# filtering read code
pure_r3 = r3[r3['description_type'] == 'P']
pure_r3 = pure_r3[pure_r3['status'] == 'C']
pure_r2 = r2[r2['term_code'] == '00']

# filtering read-icd mapping
r3_i10_processed = r3_i10[r3_i10['mapping_status'].isin(['E', 'G'])]
r2_i10_processed = r2_i10[r2_i10['icd10_code_def'] == 1]

# merge read, icd 
gp_r3_join = pd.merge(gp_r3, pure_r3, on='read_3', how='left')
gp_r2_join = pd.merge(gp_r2, pure_r2, on='read_2', how='left')
gp_r3_i10_join = pd.merge(gp_r3_join, r3_i10_processed, on='read_3', how='left')
gp_r2_i10_join = pd.merge(gp_r2_join, r2_i10_processed, on='read_2', how='left')
# reindexing
gp_r3_i10_join.index = gp_r3_index
gp_r2_i10_join.index = gp_r2_index
# concat read2, read3 / drop unnecessary columns / index sorting
read_gp = pd.concat([gp_r3_i10_join, gp_r2_i10_join], axis=0)
read_gp_processed = read_gp.drop(['icd10_code_def', 'add_code_flag', 'element_num', 'block_num', 'term_code', 'description_type', 'status'], axis=1)
read_gp_processed.sort_index(axis=0, inplace=True)
# convert date type
read_gp_processed['event_dt'] = pd.to_datetime(read_gp_processed['event_dt'])


### make DB and check
con = sqlite3.connect("/home/lrainsoul/GP/gp.db")
sqlite3.Connection
cursor = con.cursor()
read_gp_processed.to_sql('gp_clinical', con) # table name = gp_clinical

## sql check example
cursor.execute("SELECT * FROM gp_clinical WHERE eid=1004070 ORDER BY event_dt")
rows = cursor.fetchall() 
for row in rows: 
    print(row)
con.close()


### Read code description counts
# make Read3_sorted_counts.csv
read_counts = read_gp_processed['read_3'].value_counts()
read_counts_df = pd.DataFrame(read_counts)
read_counts_df['code'] = read_counts_df.index
read_counts_df.columns = ['counts', 'read_3']
merged = pd.merge(read_counts_df, pure_r3, on='read_3', how='left')
merged.drop(['description_type', 'status'], axis=1, inplace=True)
merged.to_csv('Read3_sorted_counts.csv', index=None)

# make Read2_sorted_counts.csv
read2_counts = read_gp_processed['read_2'].value_counts()
read2_counts_df = pd.DataFrame(read2_counts)
read2_counts_df['code'] = read2_counts_df.index
read2_counts_df.columns = ['counts', 'read_2']
merged2 = pd.merge(read2_counts_df, pure_r2, on='read_2', how='left')
merged2.drop(['term_code'], axis=1, inplace=True)
merged2.to_csv('sorted2_counts.csv', index=None)

### Type2 diabetes

# READ 2
# Non-insulin dependent diabetes mellitus - C109.
# T2DB - C10F.
# Pre-existing type 2 diabetes mellitus in pregnancy - L180B
# Insulin treated Type II diabetes mellitus - C109J
# Insulin treated Type II diabetes mellitus - C10FJ

# READ 3
# Insulin treated Type 2 diabetes mellitus - X40J6

r2_t2db_id = read_gp_processed[read_gp_processed['read_2'].isin(['C10F.', 'C109.'])]['eid'].unique()
r3_t2db_id = read_gp_processed[read_gp_processed['read_3'].isin(['X40J5'])]['eid'].unique()
t2db_all_id = np.union1d(r2_t2db_id, r3_t2db_id)
t2db = read_gp_processed[read_gp_processed['eid'].isin(t2db_all_id)]
np.sum(t2db['eid'].value_counts() >= 100) # record check, count > 100 record 

## complication
# read 3 (Non-insulin dependent~)
# C1090 renal complication
# C1091 ophthalmic complication
# C1092 neurological complication
# C1093 multiple complication
# C1094 with ulcer
# C1095 gangrene
# C1096 retinopathy
# C1097 poor control

# XaELQ without complication
# XaFmA diabetic cataract
# XalfG t2db on insulin
# Xalfl t2db on diet only
# ...

# read 2
# C10F0 retinal
# C10F1 ophthalamic
# C10F2 neurlogical
# C10F3 multiple
# C10F4 ulcer
# C10F5 gangree
# C10F6 retinopathy
# C10F7 poor control 
# C10F9 without complication
# C10FE diabetic cataract
# ...
# (Non-insulin dependent~)
# C1090 renal
# C1091 ophthalamic 
# ...

# ophthalmic complication
oph3_id = t2db[t2db['read_3'].isin(['C1091'])]['eid'].unique()
print('read 3:', len(oph3_id)) # result=6
oph2_id = t2db[t2db['read_2'].isin(['C10F1', 'C1091'])]['eid'].unique()
print('read 2:', len(oph2_id)) # result=0

# diabetic cataract
catar3_id = t2db[t2db['read_3'].isin(['XaFmA'])]['eid'].unique()
print('read 3:', len(catar3_id)) # result=6
catar2_id = t2db[t2db['read_2'].isin(['C10FE', 'C109E'])]['eid'].unique()
print('read 2:', len(catar2_id)) # result=4

