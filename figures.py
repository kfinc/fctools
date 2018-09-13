import seaborn as sns

def matrix_networks_plot(M, network_colors, dpi = 300, colorbar = False, group = None, ses = None, suffix = None, out_dir = None):
    """Creates and saves matrixplot with networks color labels. """
    g = sns.clustermap(M,
                   cmap="RdBu_r",
                   row_cluster=False,
                   col_cluster=False,
                   row_colors=network_colors,
                   #col_colors=network_colors,
                   linewidths=0,
                   yticklabels=False,
                   xticklabels=False)

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([.97, .2, .03, .45])
    g.ax_heatmap.set_title(f'{group}: {ses}', size = 15)
    g.cax.set_visible(colorbar)
    
    if out_dir == None:
        "Figure not saved"
    else:
    
        if suffix != None:
            g.savefig(f'{out_dir}{group}_{ses}_{suffux}.eps', dpi=dpi)
        else:
            g.savefig(f'{out_dir}{group}_{ses}.eps', dpi=dpi)
