import pandas as pd
df = pd.read_csv('metadata.tsv', sep='\t')
df = df[['Experiment accession', 'Biosample term name',
         'Biosample treatments', 'Biosample treatments duration']]
df = df.drop_duplicates()

# merge with preexising shorthand names
temp = pd.read_csv('biosamp_term_name_map.tsv', sep='\t', header=None,
                   names=['Experiment accession', 'Biosample term name',
                          'shorthand'])
temp.drop('Biosample term name', axis=1, inplace=True)

df = df.merge(temp, on='Experiment accession', how='left')
df = df.sort_values(by='Biosample term name')
df.to_csv('temp_biosamp_term_name_map.tsv', sep='\t', index=False, header=None)
