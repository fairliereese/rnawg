import pandas as pd

# table of experiment <-> tissue
df = pd.read_csv('metadata.tsv', sep='\t')
df = df[['Experiment accession', 'Biosample term name']]
df.drop_duplicates(inplace=True)
df['formatted_name'] = '| '+df['Biosample term name']+' | '+df['Experiment accession']+' |'
df = df['formatted_name']
df.to_csv('biosample_exp_table.txt', index=False, header=False)

# commands to make the file accession <-> sample ID
df = pd.read_csv('metadata.tsv', sep='\t')
df = df[['File accession', 'Experiment accession', 'Biosample term name', 'Technical replicate(s)', 'Biological replicate(s)']]
map = pd.read_csv('biosamp_term_name_map.tsv', sep='\t', header=None, usecols=[1,2])
map.columns = ['biosamp_name', 'new_biosamp_name']
map.drop_duplicates(inplace=True)
df = df.merge(map, how='left', left_on='Biosample term name', right_on='biosamp_name')
df['new_biosamp_name'] = df.new_biosamp_name.str.lower()
df['new_biosamp_name'] = df.new_biosamp_name.str.replace(' ', '_')
temp = df[['Experiment accession', 'biosamp_name', 'File accession']].groupby(['Experiment accession', 'biosamp_name']).count().reset_index()
temp['biorep'] = temp.groupby('biosamp_name').cumcount()+1
temp = temp[['Experiment accession', 'biorep']]
temp.biorep = temp.biorep.astype(str)
df = df.merge(temp, on='Experiment accession')
df['hr'] = df.new_biosamp_name+'_'+df.biorep+'_'+df['Biological replicate(s)'].astype(str)
df = df[['File accession', 'hr']]
df.to_csv('file_to_hr.tsv', sep='\t', header=None, index=False)

# make the samples file
df = df['biosamp_name']
df.to_csv('samples.txt', header=None, index=None)
