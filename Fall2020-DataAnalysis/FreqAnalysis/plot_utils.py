import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def add_stats(data,sig_matrix,ax,type='boxplot'):
    """
    Add asterisks to the plot to indicate significance; makes it so the bars and
     *s don't overlap
    Inputs:
        data - depending on the type, this is used to ensure the significance lines
            do not overlap with the data. For boxplot, format it in a list of numpy arrays.
            For bar, it is a list of the max value on the plot
        sig_matrix - numpy array with a row for each significant pair. The columns
            are as follows [factor1(0-N) factor2(0-N) p-value]
        ax - the plot axes to place the significance lines
        type - indicated what type of plot you are adding to-- this is used to make
            sure the lines do not overlap with the plot
    Outputs: n/a
    """

    # Checks to make sure the plot exists
    if sig_matrix.shape[0]==0:
        return

    # Matrix for storing where lines have been placed above every column
    sig_fill = np.zeros((sig_matrix.shape[0],len(data)))

    # Set variables for spacing lines
    [ymin,ymax]=ax.get_ylim()
    min_dist = (ymax-ymin)/50.0 # indicates the minimum distance from the data
    tol = min_dist # indicates the minimum spacing between two lines

    # Iterates through each significant p-value (row in sig_matrix)
    for i in range(sig_matrix.shape[0]):

        # obtain values from sig_matrix
        fac1 = int(sig_matrix[i,0])
        fac2 = int(sig_matrix[i,1])
        pval = sig_matrix[i,2]

        # determine whether to step up or down when iterating
        if fac1>fac2:
            temp = fac2
            fac2 = fac1
            fac1 = temp
        step=1

        # Find the minimum height of the line so that is doesn't overlap with the data
        min_liney = 0
        for ii in range(fac1,fac2+1,step):
            if type=='boxplot':
                ii_max = np.max(data[ii])
            elif type=='bar':
                ii_max = data[ii]
            if ii_max>min_liney:
                min_liney = ii_max
        min_liney += min_dist + tol

        # increase the location of y_topline until one is suitable
        notfound = True
        y_topline = min_liney
        if i==0:
            notfound=False
        while notfound:
            notfound=False
            for ii in range(fac1,fac2+1): #iterate though overlapping trial conditions
                for j in range(i): # iterate through all of the line locations
                    if abs(y_topline-sig_fill[j,ii])<tol: #This line location is not suitable
                        y_topline += tol
                        notfound = True
                        break


        # fill array with chosen locations
        for ii in range(fac1,fac2+1):
            sig_fill[i,ii] = y_topline

        # Plot the line
        # get points for line notches
        if type=='boxplot':
            fac1_max = np.max(data[fac1])
            fac2_max = np.max(data[fac2])
        elif type=='bar':
            fac1_max = data[fac1]
            fac2_max = data[fac2]

        ax.plot([fac1,fac1,fac2,fac2],
                [fac1_max+min_dist,y_topline,y_topline,fac2_max+min_dist],
                 '-k', markersize=5)


        xloc = fac1+((fac2-fac1)/2.0)
        yloc = y_topline - (tol/3)

        if pval<0.001:
            ax.text(xloc,yloc,'***',horizontalalignment='center',fontsize=8,fontweight='bold')
        elif pval<0.01:
            ax.text(xloc,yloc,'**',horizontalalignment='center',fontsize=8,fontweight='bold')
        elif pval<0.05:
            ax.text(xloc,yloc,'*',horizontalalignment='center',fontsize=8,fontweight='bold')

        # if i==3:
        #     return
    return

def add_labels(ax,x1,x2,y,name,text_buffer):
    """
    Adds arrows and labels to the bottom of plot by providing the locations in
    plot coordinates
    Inputs:
        ax - indicated plot axis
        x1 - list of x-values for the starting point of the arrows
        x2 - list of x-values for the ending point of the arrows
        y - float for the y-location of the arrow
        name - list of names for the arrows
        buffer - the size of the buffer between the arrow and text
    Outputs: n/a
    """

    text_height = y-text_buffer
    for i in range(len(x1)):
        text_x = x1[i] + (x2[i]-x1[i])/2
        ax.annotate('', xy=(x1[i],y),xytext=(x2[i],y),                     #draws an arrow from one set of coordinates to the other
                    arrowprops=dict(arrowstyle='<|-|>',facecolor='black'),   #sets style of arrow and colour
                    annotation_clip=False)                               #This enables the arrow to be outside of the plot
        ax.text(text_x, text_height,name[i],horizontalalignment='center', fontname="sans-serif", fontsize=10)

    return

def make_stackedbar(data1,data2,title,xlabel,ylabel,labels,colors,alphas,figure_size):
    """
    Creates and formats a stacked bar plot
    Inputs:
        data1 - The original data for the lower bars formatted it in a list of numpy arrays
        data2 - The original data for the upper bars formatted it in a list of numpy arrays
        title/xlabel/ylabel - strings for each label
        labels - list of strings corresponding to each dataset
        colors - list of strings corresponding to the desired color for each dataset
        alphas - list of floats from 0 to 1 corresponding to the desired
            tranparency for each list item in both datasets
        figure_size - tuple for figure size
    Outputs: the fig and ax labels for the figure and the upper_data_bound for
        each condition on the plot for adding significant asterisks
    """
    fig, ax = plt.subplots(figsize=figure_size,dpi=300)

    n = len(data1)
    n2 = len(data2)
    if n!=n:
        print('Data not of the same length')

    # Defines x-axis values
    ind = np.arange(n)

    width = 0.5 # width of each plot

    # plot data bars sequentially
    upper_data_bound = [] # store the upper bar for stat testing
    for i in range(n):
        data1_mean=np.mean(data1[i])
        data1_std=np.std(data1[i])/np.sqrt(data1[i].shape[0])

        data2_mean=np.mean(data2[i])
        data2_std=np.std(data2[i])/np.sqrt(data2[i].shape[0])

        upper_data_bound.append(data1_mean+data2_mean+data2_std)

        if alphas[i]==1:
            p1 = ax.bar(ind[i],data1_mean,width,yerr=data1_std,
                ecolor='black',capsize=5,color=colors[0],alpha=alphas[i])
            p2 = ax.bar(ind[i],data2_mean,width,bottom=data1_mean,yerr=data2_std,
                ecolor='black',capsize=5,color=colors[1],alpha=alphas[i])
        else:
            ax.bar(ind[i],data1_mean,width,yerr=data1_std,
                ecolor='black',capsize=5,color=colors[0],alpha=alphas[i])
            ax.bar(ind[i],data2_mean,width,bottom=data1_mean,yerr=data2_std,
                ecolor='black',capsize=5,color=colors[1],alpha=alphas[i])

    # place grid in back
    ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                   alpha=0.5)
    ax.set_axisbelow(True)

    # Add titles and labels
    plt.xlabel(xlabel,fontname="sans-serif", fontsize=11)
    plt.ylabel(ylabel,fontname="sans-serif", fontsize=11)
    plt.title(title,fontname="sans-serif", fontsize=11,fontweight='bold')
    for label in (ax.get_yticklabels()):
        label.set_fontsize(8)

    #figure legend
    L = fig.legend([p1[0], p2[0]], ['Lives', 'Treasures'],ncol=2, fontsize=10,loc='upper center',
        bbox_to_anchor=(0.46, .79, 0.1, 0.1))
    plt.setp(L.texts, family='Arial')

    # x-ticks x-axis
    plt.xticks(ind, labels, fontname="sans-serif", fontsize=10)
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)

    return [fig,ax,upper_data_bound]

def make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size):

    """
    Creates and formats boxplot
    Inputs:
        data - The original data formatted it in a list of numpy arrays
        title/xlabel/ylabel - strings for each label
        labels - list of strings corresponding to each list item in data
        box_colors - list of strings corresponding to the desired color for each
            list item in data
        box_alpha - list of floats from 0 to 1 corresponding to the desired
            tranparency for each list item in data
        figure_size - tuple for figure size
    Outputs: the fig and ax labels for the figure
    """

    fig, ax = plt.subplots(figsize=figure_size,dpi=300)
    # medianprops = dict(linewidth=2.5, color='black')
    medianprops = dict(linewidth=1.5, color='black')
    bp = ax.boxplot(data, notch=0, medianprops=medianprops, labels=labels)#, patch_artist=True,boxprops=dict(facecolor=color_combine, color=c))
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    ax.set_axisbelow(True) # Hide these grid behind plot objects
    ax.set_title(title, fontsize=10, fontweight='bold',fontname="sans-serif")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    #     label.set_fontsize(8)
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)
    fig.subplots_adjust(bottom=0.2)

    numboxes = len(data)

    medians = np.empty(numboxes)
    for i in range(numboxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax.add_patch(Polygon(box_coords, facecolor=box_colors[i], alpha=box_alpha[i],zorder=2))

    # for i in range(numboxes):
    #     x_coordinate = i+1
    #     y_mean = np.mean(data[i])
    #     y_std = np.std(data[i])/np.sqrt(data[i].shape[0])
    #     ax.plot([x_coordinate,x_coordinate],
    #             [y_mean-y_std,y_mean+y_std]
    #             ,'black',markersize=7,zorder=3)
    #     ax.plot(x_coordinate, y_mean, 'o',
    #                  color='w', marker='o', markersize=7, markeredgecolor='black',zorder=4)#, linewidth=0)

    return [fig,ax]
