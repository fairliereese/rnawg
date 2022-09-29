def plot_dorito(counts,
                top='splicing_ratio',
                subset=None,
                gene=None,
                hue=None,
                cmap='magma',
                mmap=None,
                density=False,
                density_scale=1,
                density_cmap='viridis',
                density_vmax=None,
                sectors=False,
                sect_alpha=0.5,
                sect_beta=0.5,
                sect_gamma=0.5,
                log_density=False,
                scatter=True,
                size=None,
                legend=True,
                log_size=False,
                jitter=False,
                alpha=1,
                scale=True,
                title=None,
                opref='figures/'):
    """
    Plot a dorito from counts with the given subset in a given 
    color
    
    Parameters:
        counts (pandas DataFrame): DF of the counts per gene 
            of ic, tss, tes from get_ic_tss_tes or
            compute_triplets (or both!!!)
        top (str): Column name to plot as apex of dorito.
            Choose from 'ic' or 'splicing_ratio'
        subset (dict of lists): List mapping counts column names
            to values in said columns to include in the data
        hue (str): Column from counts to color by
        cmap (str or dict of str): Either a dictionary mapping
            categorical column values from hue or a valid 
            matplotlib continuous named color palette
        mmap (str or dict of str): Dictionary mapping categorical
            column values from hue to marker styles
        scale (bool): Whether to scale values b/w 1 and 0. 
        alpha (float): Alpha value of points
        title (str): Title to give plot
        opref (str): Output file prefix to save fig
    """
    
    #### subset dataset and transform numbers as needed ####
    temp = counts.copy(deep=True)

    # if we have a gene name, limit to those entries
    if gene:

        temp = temp.loc[temp.gname == gene]
        if temp.empty:
            temp = temp.loc[temp.gid == gene]
            gene = temp.gname.values[0]

    # if we have a list of allowed sources, limit to those entries
    if subset:
        for col, val in subset.items():
            if type(val) != list:
                val = [val]
            temp = temp.loc[temp[col].isin(val)]
            
    # scale and assign which columns to use 
    c = dict()
    if scale:
        if top == 'splicing_ratio':
            temp['total'] = temp.n_tss+temp.n_tes+temp.splicing_ratio
        elif top == 'n_ic':
            temp['total'] = temp.n_tss+temp.n_tes+temp.n_ic
        temp['tss_ratio'] = temp.n_tss/temp.total
        temp['tes_ratio'] = temp.n_tes/temp.total
        temp['top_ratio'] = temp[top]/temp.total

        c['a'] = 'tss_ratio'
        c['b'] = 'top_ratio'
        c['c'] = 'tes_ratio'
    else:
        c['a'] = 'n_tss'
        c['b'] = top
        c['c'] = 'n_tes'
    
    if scale == True:
        scale = 1
        mult = 0.2
    else: 
        scale = max_pts(temp, c)
        
    # density
    if density:
        if hue:
            if counts[hue].dtype.name == 'object':
                pad = 0.1
            else:
                pad = 0.0
        else:
            pad = 0.1
        figure, tax, temp = density_dorito(temp, c, 
                                 density_scale, 
                                 density_cmap, 
                                 density_vmax,
                                 log_density,
                                 pad=pad)
        scale = density_scale
        figure.set_size_inches(13,10)
        
    # if we're jittering, adjust the points for each thing
    if jitter:
        temp, c = jitter_dorito(temp, c, density_scale) 

    # figure layout parameters
    fontsize = 18
    offset = 0.1
    mult = scale/5

    # if we don't already have a fig and axis from density,
    # make one
    if not density:
        figure, tax = ternary.figure(scale=scale, permutation='210')
        figure.set_facecolor('white')

    # plot gridlines below the scatterplot
    tax.gridlines(linewidth=3, multiple=mult,
                  color='white', zorder=1, linestyle=None)
    
    # scatter
    if scatter:

        figure, tax = scatter_dorito(temp, c, hue,
                                    size, log_size,
                                    cmap, mmap, alpha,
                                    density, legend,
                                    figure, tax)
        
    # sectors
    if sectors:
        line_dorito(sect_alpha, sect_beta, sect_gamma, 
                    scale, tax, figure)

    # title handler
    if not title:
        if gene:
            title = '$\it{}$\n'.format(gene)
        else:
            title = ''
    else:
        if gene:
            title = '{} $\it{}$\n'.format(title, gene)
        else:
            title = '{}\n'.format(title)

    tax.set_title(title, fontsize=20)
    tax.boundary(linewidth=2, c='#e5ecf6')
    labels = ['{:.1f}'.format(n) for n in np.arange(0, 1.2, .2)]
    tax.ticks(ticks=labels,
              axis='lbr', linewidth=1, multiple=mult,
              tick_formats="%.1f", offset=0.014,
              fontsize=14)
    # tax.ticks(axis='lbr', linewidth=1, multiple=mult,
    #           tick_formats="%.1f", offset=0.014,
    #           fontsize=14)
    
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    tax.set_background_color('#e5ecf6')

    if top == 'splicing_ratio':
        top_label = 'Splicing ratio $\\beta$'
    elif top == 'intron_chain':
        top_label = 'Intron chains $\\delta$'
    # tax.left_corner_label('# TSSs $\\alpha$', fontsize=fontsize)
    # tax.top_corner_label(top_label, fontsize=fontsize)
    # tax.right_corner_label('# TESs $\\gamma$', fontsize=fontsize)
    tax.left_axis_label('TSS $\\alpha$', fontsize=fontsize, offset=0.12)
    tax.right_axis_label(top_label, fontsize=fontsize, offset=0.12)
    tax.bottom_axis_label('TES $\\gamma$', fontsize=fontsize, offset=0.00)
    
    figure.set_facecolor('white')
        
    # tax.show()
    
    # save figure
    fname = opref
    if gene:
        fname += '_{}'.format(gene)
    if density:
        fname += '_density'
    if scatter:
        fname += '_scatter'
    if hue:
        fname += '_{}'.format(hue)
    fname += '.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight')               
    
    return temp

def scatter_dorito(counts,
                   c, 
                   hue,
                   size,
                   log_size,
                   cmap, 
                   mmap,
                   alpha,
                   density,
                   legend,
                   figure, 
                   tax):
    """
    Parameters 
        counts (pandas DataFrame): subset the thing
        c (dict of str): Dictionary of column names to plot as a, b, c 
            indexed by 'a', 'b', 'c'
    """
    
    def scale_col(points, counts, col, log=False, how='color'):
            if log:
                log_col = '{}_log'.format(col)
                counts[log_col] = np.log10(counts[col])
                col = log_col
            vals = counts[[col]]
            max_val = vals[col].max()
            min_val = vals[col].min()
            min_max_scaler = preprocessing.MinMaxScaler(feature_range=(10, 300))
            vals = min_max_scaler.fit_transform(vals)
            max_scaled = max(vals)
            min_scaled = min(vals)
            
            # replace nans w/ 100
            vals = [100 if np.isnan(v) else v for v in vals]
            
            return vals, min_val, max_val, min_scaled, max_scaled
        
    # defaults
    points = [(x[0], x[1], x[2]) for x in zip_pts(counts, c)]
    labels = ['' for i in range(len(points))]
    hue_type = None
    figsize = (10,10)
    colors = '#e78424'
    if len(points) < 60:
        sizes = [100 for i in range(len(points))]
    else:
        sizes =  [20 for i in range(len(points))]
    markers = 'o'
    vmin = 0
    vmax = 1
    plotted = False
        
    # get color
    if hue:
            
        # categorical
        if counts[hue].dtype.name == 'object':
            hue_type = 'cat'
            colors = counts[hue].map(cmap).tolist()
            labels = counts[hue].tolist()

        # continuous
        else:
            hue_type = 'cont'
            colors, abs_min, abs_max, vmin, vmax = scale_col(points, counts, hue)
    
    # get sizes
    if size:
        sizes,_,_,_,_ = scale_col(points, counts, size, log_size)
        print(sizes[:5])
        
    # marker style
    if mmap:
        markers = [mmap[val] if val in mmap.keys() else 'o' for val in counts[hue].unique()]
    elif hue_type == 'cat':
        markers = ['o' for val in counts[hue].unique()]
        
    # figure size handling
    if hue_type == 'cat' and density: figsize = (13,10)
    elif hue_type == 'cat' and not density: figsize = (10,10)
    elif hue_type == 'cont' and density: figsize = (16,10)
    elif hue_type == 'cont' and not density: figsize = (13,10)
    elif density: figsize = (13,10)
    figure.set_size_inches(figsize[0], figsize[1])
    
    # actual scatter call
    if hue_type == 'cat':
        pdb.set_trace()
        for point, color, size, label, marker in zip(points, colors, sizes, labels, markers):
            tax.scatter([point], vmin=vmin, vmax=vmax,
                    s=size, c=color, cmap=cmap,
                    marker=marker,label=label,
                    alpha=alpha, zorder=3)
    else:   
        tax.scatter(points, vmin=vmin, vmax=vmax,
                    s=sizes, c=colors, cmap=cmap, marker=markers,
                    alpha=alpha, zorder=3)
    
    # legend handling
    if hue_type == 'cat' and legend:
        if density: x = 1.6
        else: x = 1.4
        tax.legend(bbox_to_anchor=(x, 1.05),
                   loc='upper right', prop={'size': 14})
        
        # fix marker size
        ax = tax.get_axes()
        lgnd = ax.get_legend()
        for handle in lgnd.legendHandles:
            handle._sizes = [100]
    
    # colorbar handling
    if hue_type == 'cont':
        ax = tax.get_axes()
        norm = plt.Normalize(vmin=abs_min, vmax=abs_max)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        cb = plt.colorbar(sm, ax=ax, pad=0.1)
        for t in cb.ax.get_yticklabels():
             t.set_fontsize(16)
        if hue == 'tss' or hue == 'tes':
            cb.set_label('# {}s'.format(hue.upper()), size=16)
        elif hue == 'intron_chain':
            cb.set_label('# {}s'.format(hue), size=16)
        
    return figure, tax  