import matplotlib.pyplot as plt

# This code has functions to read in and plot data
# Or add new data to existing data (ax should not be provided in this case)
# Parameters:
#   - ax: if ax is provided, then plot additional data points onto the ax
#   - pc_index_1/pc_index_2: PCs to plot
#   - skip_1st_line: if input file has header, then need to skip the header
# Return:
#   - fig, ax for further modification or saving purpose
def plot_pca(input_fn, ax=None, pc_index_1=1, pc_index_2=2, skpi_1st_line = True,
             figure_title = None, subset_label = None, subset_color = None):
    # Plot pc 1 and 2 by default, can change this later
    pc_index_1 = pc_index_1
    pc_index_2 = pc_index_2
    # Empty list to store pc values for plotting
    pc1_lst = []
    pc2_lst = []
    with open(input_fn, 'r') as fh:
        if skpi_1st_line:
            line = fh.readline().strip()    # Skip the first header line
        line = fh.readline().strip()
        while line != '':
            # Index needs to + 1 in .eigenve file
            pc1 = line.split()[pc_index_1 + 1]
            pc2 = line.split()[pc_index_2 + 1]
            pc1_lst.append(float(pc1))
            pc2_lst.append(float(pc2))
            line = fh.readline().strip()

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
        ax.plot(pc1_lst, pc2_lst, 'bo', markersize=2, alpha=0.3,
                markeredgecolor='k', markeredgewidth=1)
        if figure_title is None:
            figure_title = 'PC' + str(pc_index_1) + ' and ' + 'PC' + str(pc_index_2)
        else:
            figure_title = figure_title
        ax.set_title(figure_title)
        ax.set_xlabel('PC' + str(pc_index_1))
        ax.set_ylabel('PC' + str(pc_index_2))
        return fig, ax
    else:  # Add data to existing ax
        if subset_label is None:
            subset_label = 'Subset'
        else:
            subset_label = subset_label

        if subset_color is None:
            subset_color = 'r'
        else:
            subset_color = subset_color
        ax.plot(pc1_lst, pc2_lst, color=subset_color, marker='o', linestyle=None,
                markersize=2, alpha=0.5, markeredgecolor='r', markeredgewidth=1,
                label=subset_label)
        ax.legend()