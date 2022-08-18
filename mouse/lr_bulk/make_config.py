# TODO - have this sort by things that are already in the db so that
# I can make config files for everything that's not

import pandas as pd
import glob
import os

df = pd.read_csv('samples.txt', header=None, names=['name'])
df['fname'] = 'processing/'+df.name+'_labeled.sam'
df['name'] = df['name'].str.replace('-', '_')
df['sample'] = df['name']
df['platform'] = 'PacBio'
df = df[['name', 'sample', 'platform', 'fname']]
n_samples = len(df.index)

fname = 'talon/talon_config.csv'
# append if we already have a config file
if os.path.exists(fname):
    df.to_csv(fname, header=False, index=False, mode='a')
else:
    df.to_csv(fname, header=False, index=False)

# decide what index is going to be next
inds = []
for f in glob.glob('talon/talon_config_*.csv'):
    ind = f.rsplit('_')[-1]
    ind = ind.split('.')[0]
    if ind != 'back':
        inds.append(int(ind))

if len(inds) == 0:
    i = 0
else:
    i = max(inds)

# more manageable chunks
inds = [i for i in range(n_samples) if i % 14 == 0]
if len(df.index) not in inds:
    inds += [len(df.index)]

for start, end in zip(inds[:-1], inds[1:]):
    temp = df.iloc[start:end]
    i += 1
    fname = 'talon/talon_config_{}.csv'.format(i)
    temp.to_csv(fname, header=False, index=False)
