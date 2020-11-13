import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def add_stats(ax,data,sig_matrix):

    # takes
        # data (formatted for boxplot)-- this is used to find the maximum value so that the significance lines are placed above
        # sig_matrix -- numpy array with row=each significant pair
                                        #col: factor1(0-N) factor2(0-N) p-value
    print(len(sig_matrix.shape))
    if sig_matrix.shape[0]==0:
        # print('exited')
        return
    elif len(sig_matrix.shape)==1: # The numpy array is 1D
    # add dummy row to make numpy matrix 2D
        sig_matrix = np.append(sig_matrix,[0,0,5])
        sig_matrix = np.reshape(sig_matrix,(2,3))

    [ymin,ymax]=ax.get_ylim()
    min_dist = (ymax-ymin)/50.0
    sig_fill = np.zeros(10)
    for i in range(sig_matrix.shape[0]):
        fac1 = int(sig_matrix[i,0])
        fac2 = int(sig_matrix[i,1])
        pval = sig_matrix[i,2]
        if pval<0.05:
            fac1_max = np.max(data[fac1])
            fac2_max = np.max(data[fac2])
            if fac1_max>fac2_max:
                y_topline = fac1_max+min_dist*2
            else:
                y_topline = fac2_max+min_dist*2

            # decides how far up to place the sig line so that there is no overlap
            level1 = sig_fill[fac1]
            sig_fill[fac1] += 1
            level2 = sig_fill[fac2]
            sig_fill[fac2] += 1
            level = np.max([level1,level2])
            # found1 = False
            # found2 = False
            # while found1==False and found2==False:
            #
            #     if found1==False and sig_fill[fac1]!=level1:
            #         level1+=1
            #     else:
            #         found1 = True
            #
            #      and sig_fill[level,fac2]==0:
            #         sig_fill[level,fac1] = 1
            #         sig_fill[level,fac2] = 1
            #         found = True
            #     level+=1
            y_topline += min_dist * level

            # ax.plot([fac2+1,fac2+1],
            #         [y_topline,fac2_max+min_dist],
            #          '-k', markersize=5)
            # ax.plot([fac1+1,fac1+1,fac2+1,fac2+1],
            #         [fac1_max+min_dist,y_topline,y_topline,fac2_max+min_dist],
            #          '-k', markersize=5)
            xloc = fac1+1.0+((fac2-fac1)/2.0)
            yloc = y_topline - (min_dist/3)
            if pval<0.001:
                ax.text(xloc,yloc,'***',horizontalalignment='center',fontsize=8,fontweight='bold')
                # ax.plot([fac1+1,fac1+1,fac2+1,fac2+1],
                #         [fac1_max+min_dist,y_topline,y_topline,fac2_max+min_dist],
                #          '-k', markersize=5)
            elif pval<0.01:
                ax.text(xloc,yloc,'**',horizontalalignment='center',fontsize=8,fontweight='bold')
                # ax.plot([fac1+1,fac1+1,fac2+1,fac2+1],
                #         [fac1_max+min_dist,y_topline,y_topline,fac2_max+min_dist],
                #          '-k', markersize=5)
            elif pval<0.05:
                ax.text(xloc,yloc,'*',horizontalalignment='center',fontsize=8,fontweight='bold')
            ax.plot([fac1+1,fac1+1,fac2+1,fac2+1],
                    [fac1_max+min_dist,y_topline,y_topline,fac2_max+min_dist],
                     '-k', markersize=5)

    return

def make_boxplot(fig,ax,data,labels,box_colors,box_alpha):
    # data is a python list

    medianprops = dict(linewidth=2.5, color='black')
    bp = ax.boxplot(data, notch=0, medianprops=medianprops, labels=labels, widths = 0.6)#, patch_artist=True,boxprops=dict(facecolor=color_combine, color=c))
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    ax.set_axisbelow(True) # Hide these grid behind plot objects
    # for tick in ax.get_xticklabels():
    #     tick.set_rotation(45)
    # fig.subplots_adjust(bottom=0.2)

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
        ax.add_patch(Polygon(box_coords, facecolor=box_colors[i], alpha=box_alpha[i]))

    x_coordinates = np.arange(1,numboxes+1)
    for i in range(numboxes):
        y_coordinate = np.mean(data[i])
        std_i = np.std(data[i])/np.sqrt(len(data[i]))
        ax.plot([x_coordinates[i],x_coordinates[i]],
                [y_coordinate-std_i,y_coordinate+std_i]
                ,color='w',markersize=15)
        ax.plot(x_coordinates[i], y_coordinate, 'o',
                 color='w', marker='o', markersize=7, markeredgecolor='black')#, linewidth=0)
    return
