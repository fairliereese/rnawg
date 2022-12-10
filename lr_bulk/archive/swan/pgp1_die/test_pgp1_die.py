import swan_vis as swan

sg = swan.read('swan.p')

# run die tests between each set of pgp1 things
pgp1_conds = sg.adata.obs.loc[sg.adata.obs['sample'].str.contains('pgp1'), 'sample'].unique().tolist()
tested = []
for c1 in pgp1_conds:
    for c2 in pgp1_conds:
        test1 = (c1, c2)
        test2 = (c2, c1)
        if test1 in tested or test2 in tested:
            continue
        elif c1 == c2:
            continue
        else:
            print()
            print(c1)
            print(c2)
            tested.append(test1)
            tested.append(test2)
            fname = 'iso_die_{}_{}.tsv'.format(c1, c2)
            die_table, test_results = sg.die_gene_test(kind='iso',
                                                       obs_col='sample',
                                                       obs_conditions=[c1, c2],
                                                       verbose=True)
            die_table.to_csv(fname, sep='\t')
