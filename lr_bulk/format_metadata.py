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
df = df[['File accession', 'Experiment accession', 'Biosample term name', 'Technical replicate(s)']]
map = pd.read_csv('biosamp_term_name_map.tsv', sep='\t', header=None, usecols=[1,2])
map.columns = ['biosamp_name', 'new_biosamp_name']
map.drop_duplicates(inplace=True)
df = df.merge(map, left_on='Biosample term name', right_on='biosamp_name')
df['new_biosamp_name'] = df.new_biosamp_name.str.lower()
df = df[['File accession', 'new_biosamp_name', 'Technical replicate(s)']]
df['biosamp_name'] = df.new_biosamp_name+'_'+df['Technical replicate(s)']
df = df[['File accession', 'biosamp_name']]
df.to_csv('file_to_hr.tsv', sep='\t', header=None, index=False)
