#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import math

final_t2 = pd.read_csv('filtered_final_t2.csv')
counts = pd.read_csv('sorted_counts_t2db.csv')
counts_r2 = pd.read_csv('sorted_counts_t2db_read2.csv')
hesin = pd.read_csv('hesin_diag_main.csv')
uniq = final_t2['eid'].unique()
hesin_t2 = hesin[hesin['eid'].isin(uniq)]
hesin_t2.drop(['diag_icd9_nb', 'diag_icd10_nb'], axis=1, inplace=True)

str_filter = final_t2[final_t2['value1'].isin(['^'])].index
final_t2.drop(str_filter, axis=0, inplace=True)

# final_t2.drop(2813483, axis=0, inplace=True)
temp = final_t2.copy() # for all data
# temp = final_t2[final_t2['eid'].isin([3452033])]

temp_bp = temp[(temp['read_2'].isin(['246..', '22A..', '42W4.', '42W5.', '44J3.', '22K..', '44Q..', '42J..', '42M..', '44P5.', '44P6.',
                                    '451E.', '42H..', '44P..'])) 
               | (temp['read_3'].isin(['2469.', '246A.', '22A..', 'XaERp', 'XaPbt', 'XE2q5', '22K..', 'XE2q9', '42J..', '42M..',
                                      '44P5.', '44P6.', 'XaK8y', 'XaIdY', 'XE2eD']))]
temp_bp.replace('0.000', np.NaN, inplace=True)

# col = systolic / diastolic / weight / hba1c(%) / creatinine level / BMI / triglycerides / neutrophil / lymphocyte / HDL / LDL
#        / eGFR / WBC / total Cholesterol / 
# + diabetes check / time delta /
# after + // diabetes duration / 
cols = [['246..', '2469.'], ['246..', '246A.'], ['22A..', '22A..'], ['42W4.', 'XaERp'], ['44J3.', 'XE2q5'], ['22K..', '22K..'],
       ['44Q..', 'XE2q9'], ['42J..', '42J..'], ['42M..', '42M..'], ['44P5.', '44P5.'], ['44P6.', '44P6.'],
       ['451E.', 'XaK8y'], [ '42H..', 'XaIdY'], ['44P..','XE2eD']]

#2 DCCT(%) = a*IFCC + b
# read3 XaPbt XaERp(%)
# read2 42W5. 42W4.(%)

# ## read3
# hb_r3 = temp_bp[temp_bp.loc[:, 'read_3'] == 'XaPbt']['value1']
# # drop outliers > 240, keep >= 20
# ii = hb_r3[hb_r3.astype(np.float) > 240].index
# hb_r3.drop(ii, axis=0, inplace=True)
# index = hb_r3[hb_r3.astype(np.float) >= 20].index
# val = hb_r3[hb_r3.astype(np.float) >= 20].values
# # convert IFCC to DCCT
# temp_bp.loc[index, 'value1'] = val.astype(np.float)*0.09148+2.152
# temp_bp['read_3'].replace('XaPbt', 'XaERp', inplace=True)


# ## read2
# hb_r2 = temp_bp[temp_bp.loc[:, 'read_2'] == '42W5.']
# # check value2, value3
# hb_val2_r2 = hb_r2[hb_r2['value1'] == 'OPR003']['value2']
# hb_val1_r2 = hb_r2.drop(hb_val2_r2.index, axis=0)['value1']
# # drop outliers > 240, keep < 20
# ii2 = hb_val1_r2[hb_val1_r2.astype(np.float) > 240].index
# hb_val1_r2.drop(ii2, axis=0, inplace=True)
# index2 = hb_val1_r2[hb_val1_r2.astype(np.float) >= 20].index
# val2 = hb_val1_r2[hb_val1_r2.astype(np.float) >= 20].values

# v2_ii2 = hb_val2_r2[hb_val2_r2.astype(np.float) > 240].index
# hb_val2_r2.drop(v2_ii2, axis=0, inplace=True)
# v2_index2 = hb_val2_r2[hb_val2_r2.astype(np.float) >= 20].index
# v2_val2 = hb_val2_r2[hb_val2_r2.astype(np.float) >= 20].values

# # convert IFCC to DCCT
# temp_bp.loc[index2, 'value1'] = val2.astype(np.float)*0.09148+2.152
# temp_bp.loc[v2_index2, 'value2'] = v2_val2.astype(np.float)*0.09148+2.152
# temp_bp.loc[v2_index2, 'value1'] = temp_bp.loc[v2_index2, 'value2']
# temp_bp['read_2'].replace('42W5.', '42W4.', inplace=True)

#1 drop null values
for c in cols:
    # read 2
    filter_na2 = temp_bp[temp_bp['read_2'].isin([c[0]])]
    if filter_na2.shape[0] > 0:
        r2 = 1
        filter_na2_id = filter_na2[(filter_na2['value1'].isnull() == True) & (filter_na2['value2'].isnull() == True) &
                                (filter_na2['value3'].isnull() == True)].index
        print(filter_na2_id)
    else:
        r2 = 0

    # read 3
    filter_na3 = temp_bp[temp_bp['read_3'].isin([c[1]])]
    if filter_na3.shape[0] > 0:
        r3 = 1
        filter_na3_id = filter_na3[(filter_na3['value1'].isnull() == True) & (filter_na3['value2'].isnull() == True) &
                                (filter_na3['value3'].isnull() == True)].index
    else:
        r3 = 0

    if r2*r3:
        temp_bp.drop(np.union1d(filter_na2_id, filter_na3_id), axis=0, inplace=True)
    elif not r3:
        temp_bp.drop(filter_na2_id, axis=0, inplace=True)
    elif not r2:
        temp_bp.drop(filter_na3_id, axis=0, inplace=True)


#3 make col (read code, disease, complication)
id_array = temp_bp['eid'].unique()
input_list = []
time_list = []
t2_time_list = []
ID = 0
time = True
process = 0

for i in id_array:
    slice_i = temp_bp[temp_bp['eid'].isin([i])]
    time_array = slice_i['event_dt'].unique()
    
    # save diabetes time
    slice_hes = hesin_t2[hesin_t2['eid'].isin([i])]
    if slice_hes[slice_hes['diag_icd10'].isin(['E119', 'E110', 'E111', 'E112', 'E113', 'E114', 
                                                         'E115', 'E116', 'E117', 'E118'])].shape[0] > 0:
        t2_time = slice_hes[slice_hes['diag_icd10'].isin(['E119', 'E110', 'E111', 'E112', 'E113', 'E114', 
                                                         'E115', 'E116', 'E117', 'E118'])]['event_dt'].iloc[0] # without complication
        t2_time_list.append(t2_time)
    
#     # save ophthalmic complication time
#     if slice_hes[slice_hes['diag_icd10'].isin(['E113'])].shape[0] > 0:
#         comp_time = slice_hes[slice_hes['diag_icd10'].isin(['E113'])]['event_dt'].iloc[0]    
#     else:
#         comp_time = 0
    
    ### for other complications
    if slice_hes[slice_hes['diag_icd10'].isin(['E110', 'E111', 'E112', 'E113', 'E114', 
                                                         'E115', 'E116', 'E117', 'E118'])].shape[0] > 0:
        all_comp_time = slice_hes[slice_hes['diag_icd10'].isin(['E110', 'E111', 'E112', 'E113', 'E114', 
                                                         'E115', 'E116', 'E117', 'E118'])]['event_dt'].iloc[0]    
    else:
        all_comp_time = 0    
    
    # append timedelta
    if time:
        L = np.zeros((len(time_array), len(cols)+3))
    else:
        L = np.zeros((len(time_array), len(cols)+2))
    turn = 0
    for t in time_array:
        slice_t = slice_i[slice_i['event_dt'].isin([t])]
        col = 0
        for c in cols:
            # get mean value
            filter_r2 = slice_t[slice_t['read_2'].isin([c[0]])]
            if filter_r2.shape[0] > 0:
                v_r2 = 1
                if c[1] in ['246A.']:# bp reading, 
                    dup_r2 = filter_r2['value2'].values
                else:
                    dup_r2 = filter_r2['value1'].values
                    
                if np.sum(filter_r2['value1'].isin(['OPR005', 'OPR004', 'OPR003', 'OPR002'])) > 0:
                    mask_opr = filter_r2['value1'].isin(['OPR005', 'OPR004', 'OPR003', 'OPR002'])
#                     index_opr = filter_r2[filter_r2['value1'].isin(['OPR003'])].index
                    opr_sample = filter_r2[mask_opr]
                    filter_r2.loc[mask_opr, 'value1'] = opr_sample['value2'].values
                    dup_r2 = filter_r2['value1'].values
            else:
                v_r2 = 0

            filter_r3 = slice_t[slice_t['read_3'].isin([c[1]])]
            if filter_r3.shape[0] > 0:
                v_r3 = 1
                dup_r3 = filter_r3['value1'].values
            else:
                v_r3 = 0

            if v_r2*v_r3:
                dup_r2_r3 = np.union1d(dup_r2, dup_r3)
                mean_value = np.mean(dup_r2_r3.astype(np.float))
            elif not (v_r2 + v_r3):
                mean_value = np.NaN
            elif not v_r2:
                mean_value = np.mean(dup_r3.astype(np.float))
            elif not v_r3:
                mean_value = np.mean(dup_r2.astype(np.float))

            L[turn, col] = mean_value

            col += 1

        # diabetes check
        if pd.to_datetime(t2_time) < pd.to_datetime(t):
            L[turn, -3] = 1
        else:
            L[turn, -3] = 0

#         # ophthalmic complication
#         if comp_time == 0:
#             L[turn, -3] = 0
#         else:
#             if pd.to_datetime(comp_time) < pd.to_datetime(t):
#                 L[turn, -3] = 1
#             else:
#                 L[turn, -3] = 0
            
        # all complications
        if all_comp_time == 0:
            L[turn, -2] = 0
        else:
            if pd.to_datetime(all_comp_time) < pd.to_datetime(t):
                L[turn, -2] = 1
            else:
                L[turn, -2] = 0        
    
        # timedelta
        if time:
            if turn == 0:
                L[turn, -1] = 0
            else:
                L[turn, -1] = (pd.to_datetime(time_array[turn]) - pd.to_datetime(time_array[turn -1])).days    
            
        turn += 1
    ID += 1
    
    # 4 make cutoff
    C = len(cols)
    percent = 1
    
    mask = np.sum(np.isnan(L), axis=1) < C*percent
    new_time_array, new_L = time_array[mask], L[mask]
    input_list.append(new_L)
    time_list.append(new_time_array)
    
    process += 1
    
    if process == 1000:
        print('1000...')
    
    if process == 3000:
        print('3000...')
        
    if process == 5000:
        print('5000...')

    if process == 7000:
        print('7000...')
    
    if process == 9000:
        print('9000...')
        
    if process == 11000:
        print('11000...')

# #5 drop id
# drop_save = []
# var_mask = np.sum(np.isnan(new_L), axis=0) / new_L.shape[0]
# if np.sum(var_mask < 0.5):
#     drop_save.append(해당id) 

np.save('./input_list_first_output_allcomps_new', input_list)
np.save('./time_list_allcomps_new', time_list)
np.save('./all_comp_time_new', time_list)
# np.save('./t2_time_list_allcomps', t2_time_list)
np.save('./id_list_allcomps_new', id_array)

def add_t2(input_list): # col3 = hba1c, col -3 = diag diabetes
    for i in range(len(input_list)):
        signal = 0
        hba1c_sample = input_list[i][:, 3]
        for j in range(input_list[i].shape[0]):
            if math.isnan(hba1c_sample[j]):
                continue
            
            if hba1c_sample[j] < 6.5:
                signal = 0
                continue
            
            if (hba1c_sample[j] > 6.5) & (signal == 0):
                signal = 1
            elif (hba1c_sample[j] > 6.5) & (signal == 1): # diag diabetes
#                 diag_point = j
                input_list[i][:, -3][j:] = 1
                break
    return input_list
    
def add_diag_time(input_list):
    for i in range(len(input_list)):
        temp_diag = input_list[i][:, -3] # diabetes col
        temp_time = input_list[i][:, -1] # time delta col
        empty_duration = np.zeros((len(input_list[i]))).reshape(len(input_list[i]), -1)
        zero=True
        for j in range(len(input_list[i])):
            if zero:
                if not temp_diag[j]:
                    continue
                else:
                    zero=False
                    continue
                    
            empty_duration[j] = empty_duration[j-1] + temp_time[j]
        input_list[i] = np.concatenate((input_list[i], empty_duration), axis=1)
    return input_list
             
def clip_comp_get_xy(input_list):
    y = np.zeros((len(input_list)))
    for i in range(len(input_list)):
        temp_comp = input_list[i][:, -3] # target complication col
        for j in range(len(input_list[i])):
            if not temp_comp[j]:
                continue
            else:
                input_list[i] = np.delete(input_list[i][:j+1], (-3), axis=1)
                y[i] = 1
                break
        
        if y[i] == 0:
            input_list[i] = np.delete(input_list[i][:j+1], (-3), axis=1)
        
    return input_list, y


input_add_t2 = add_t2(input_list.copy())
input_add_t2_diag = add_diag_time(input_add_t2.copy())
x, y = clip_comp_get_xy(input_add_t2_diag.copy())
np.save('./only_eye_x_all_comps_new', x)
np.save('./only_eye_y_all_comps_new', y)

x_time_list = []
for i in range(len(x)):
    x_time_list.append(time_list[i][:x[i].shape[0]])

np.save('./all_comps_time_list_new', x_time_list)
