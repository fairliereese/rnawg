import pandas as pd
import os

# table of experiment <-> tissue
df = pd.read_csv('metadata.tsv', sep='\t')
df = df[['Experiment accession', 'Biosample term name']]
df.drop_duplicates(inplace=True)
df['formatted_name'] = '| '+df['Biosample term name']+' | '+df['Experiment accession']+' |'
df = df['formatted_name']
df.to_csv('biosample_exp_table.txt', index=False, header=False)

# commands to make the file accession <-> sample ID
df = pd.read_csv('metadata.tsv', sep='\t')
df = df[['File accession', 'Experiment accession', 'Biosample term name', 'Biosample type', 'Technical replicate(s)', 'Biological replicate(s)']]
map = pd.read_csv('biosamp_term_name_map.tsv', sep='\t', header=None, usecols=[0,1,4])
map.columns = ['Experiment accession', 'Biosample term name', 'new_biosamp_name']
map.drop_duplicates(inplace=True)
df = df.merge(map, how='left', on=['Experiment accession', 'Biosample term name'])
df['new_biosamp_name'] = df.new_biosamp_name.str.lower()
df['new_biosamp_name'] = df.new_biosamp_name.str.replace(' ', '_')
# temp = df[['Experiment accession', 'biosamp_name', 'File accession']].groupby(['Experiment accession', 'biosamp_name']).count().reset_index()
# temp['biorep'] = temp.groupby('biosamp_name').cumcount()+1
temp = df[['Experiment accession', 'new_biosamp_name', 'File accession']].groupby(['Experiment accession', 'new_biosamp_name']).count().reset_index()
temp['biorep'] = temp.groupby('new_biosamp_name').cumcount()+1

temp = temp[['Experiment accession', 'biorep']]
temp.biorep = temp.biorep.astype(str)
df = df.merge(temp, on='Experiment accession')
df['techrep'] = df.groupby('Experiment accession').cumcount()+1
df['hr'] = df.new_biosamp_name+'_'+df.biorep+'_'+df.techrep.astype(str)
temp = df[['File accession', 'hr']]

try:
    os.rename('file_to_hr.tsv', 'file_to_hr_back.tsv')
except:
    pass

temp.to_csv('file_to_hr.tsv', sep='\t', header=None, index=False)

# save info about whether each dataset is a cell line or tissue
temp = df[['hr', 'Biosample type']].copy(deep=True)
d = {'cell line': 'cell_line',
     'in vitro differentiated cells': 'cell_line',
     'primary cell': 'cell_line',
     'tissue': 'tissue'}
temp['biosample_type'] = df['Biosample type'].map(d)
temp = temp[['hr', 'biosample_type']]
temp.to_csv('hr_to_biosample_type.tsv', sep='\t', index=False)

# df[['Experiment accession', 'biosamp_name', 'biorep', 'techrep', 'hr']].sort_values('Experiment accession')
# len(df.hr.unique())
# len(df['File accession'].unique())

# make the samples file
try:
    old_samples = pd.read_csv('file_to_hr_back.tsv', header=None, sep='\t')
    old_samples.columns = ['file_id', 'hr']

    # # only include old files
    # df = df.loc[~df.hr.isin(old_samples.hr.tolist())]
    # only include new files
    df = df.loc[~df['File accession'].isin(old_samples.file_id.tolist())]
except:
    pass


df = df['hr']
df.to_csv('samples.txt', header=None, index=None)
