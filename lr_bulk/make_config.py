import pandas as pd

df = pd.read_csv('samples.txt', header=None, names=['name'])
df['sample'] = df['name']
df['platform'] = 'PacBio'
df['fname'] = 'processing/'+df.name+'_labeled.sam'

fname = 'talon/talon_config.csv'
df.to_csv(fname, header=False, index=False)

# more manageable chunks
inds = [i for i in range(63) if i % 14 == 0]
if len(df.index) not in inds:
    inds += [len(df.index)]

i = 0
for start, end in zip(inds[:-1], inds[1:]):
    temp = df.iloc[start:end]
    i += 1
    fname = 'talon/talon_config_{}.csv'.format(i)
    temp.to_csv(fname, header=False, index=False)
