import numpy as np
import pandas as pd
import os

meta = pd.read_csv('metadata.tsv', sep='\t')
meta['Donor(s)'].unique()

#### timecourse data ####
df = meta.loc[meta['Donor(s)'] == '/mouse-donors/ENCDO509HIY/'].copy(deep=True)
df['timecourse'] = True
print(len(df.index))

# tissue
short_tissue_map = {'gastrocnemius': 'gastroc',
                    'heart': 'heart',
                    'adrenal gland': 'adrenal',
                    'left cerebral cortex': 'cortex',
                    'layer of hippocampus': 'hippocampus'}
df['tissue'] = df['Biosample term name'].map(short_tissue_map)

exp_df = pd.read_csv('experiment_report.tsv', sep='\t', skiprows=[0])
exp_df = exp_df[['Accession', 'Description', 'Biosample age', 'Biosample summary']]
df = df.merge(exp_df, how='left', left_on='Experiment accession', right_on='Accession')

# sex
df['sex'] = 'm'
df.loc[df['Biosample summary'].str.contains('female'), 'sex'] = 'f'

# time pt
short_time_map = {'4 days': '4d',
                  '10 days': '10d',
                  '14 days': '14d',
                  '25 days': '25d',
                  '36 days': '36d',
                  '2 months': '2mo',
                  '18-20 months': '18-20mo'}
df['age'] = df['Biosample age'].map(short_time_map)

# fix for bug
if np.nan in df.age.unique().tolist():
    print('hellow')
    df.loc[df.age.isnull(), 'age'] = '18-20mo'

# rep
df['rep'] = df['Biological replicate(s)']

df['hr'] = df['tissue']+'_'+df['age']+'_'+df['sex']+'_'+df['rep'].astype(str)
df['hr'] = df['hr'].str.replace('_', '-')

# save samples file
temp = df[['File accession', 'hr']]
samples = temp.copy(deep=True)

# also save file to convert ENCFF --> human readable
file_hr = df[['File accession', 'hr']].copy(deep=True)

#### other data (non timecourse) ####
df = meta.loc[meta['Donor(s)'] != '/mouse-donors/ENCDO509HIY/'].copy(deep=True)
df['timecourse'] = False

df = df[['File accession', 'Experiment accession', 'Biosample term name', 'Biosample type', 'Technical replicate(s)', 'Biological replicate(s)']]
m = pd.read_csv('biosamp_term_name_map.tsv', sep='\t', header=None, usecols=[0,1,2])
m = m.drop_duplicates()
m.columns = ['Experiment accession', 'Biosample term name', 'new_biosamp_name']
m.drop_duplicates(inplace=True)
df = df.merge(m, how='left', on=['Experiment accession', 'Biosample term name'])
df['new_biosamp_name'] = df.new_biosamp_name.str.lower()
df['new_biosamp_name'] = df.new_biosamp_name.str.replace(' ', '_')
temp = df[['Experiment accession', 'new_biosamp_name', 'File accession']].groupby(['Experiment accession', 'new_biosamp_name']).count().reset_index()
temp['biorep'] = temp.groupby('new_biosamp_name').cumcount()+1

temp = temp[['Experiment accession', 'biorep']]
temp.biorep = temp.biorep.astype(str)
df = df.merge(temp, on='Experiment accession')
df['techrep'] = df.groupby('Experiment accession').cumcount()+1
df['hr'] = df.new_biosamp_name+'_'+df.biorep+'_'+df.techrep.astype(str)
df['hr'] = df['hr'].str.replace('_', '-')
file_hr = pd.concat([file_hr, df[['File accession', 'hr']]])

# save info about whether each dataset is a cell line or tissue
temp = df[['hr', 'Biosample type']].copy(deep=True)
d = {'cell line': 'cell_line',
     'in vitro differentiated cells': 'cell_line',
     'primary cell': 'cell_line',
     'tissue': 'tissue'}
temp['biosample_type'] = df['Biosample type'].map(d)
temp = temp[['hr', 'biosample_type']]

samples = pd.concat([df[['File accession', 'hr']],samples])

try:
    os.rename('file_to_hr.tsv', 'file_to_hr_back.tsv')
except:
    pass

# make the samples file
try:
    old_samples = pd.read_csv('file_to_hr_back.tsv', header=None, sep='\t')
    old_samples.columns = ['file_id', 'hr']

    # only include old files
    print(len(samples.index))
    samples = samples.loc[~samples['File accession'].isin(old_samples.file_id.tolist())]
    print(len(samples.index))

except:
    pass

samples = samples[['hr']]
samples.to_csv('samples.txt', header=None, index=None)
file_hr.to_csv('file_to_hr.tsv', index=False, header=None, sep='\t')
