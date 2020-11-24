import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def xy_spectrum(w,xdata,ydata,w_resonance,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend):

    """
    w - frequency values for the x-axis of plots
    xdata - data for left subplot in the form of a list of numpy arrays
    ydata - data for right subplot in the form of a list of numpy arrays
    w_resonance - list of resonant frequencies of the ball for vertical lines
    title,xlabel,ylabel - strings for labeling the plot
    legend - list of string labels for the legend corresponding to linestyles,
            colors, and xdata/ydata (should be the same length)
    linestyles - list of string labels for line style
    colors - list of string hex colors
    ymin/ymax - plot y min and max
    ymax_pend - max y-value for pendulum line
    """
    # create plot
    figure_size = (8,3.55) # inches
    plt.figure(1, figsize=figure_size, dpi=150)
    fig, ax =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.05)
    legend_lines = []

    # plot lines in order of legend
    if len(w_resonance)>0:
        l = ax[0].plot([w_resonance[0], w_resonance[0]],[ymin, ymax_pend],linestyle=':',color='k')
        legend_lines.append(l)
    if len(xdata)>0:
        for i in range(0,len(xdata)):
            l = ax[0].plot(w,xdata[i],linestyle=linestyles[i],color=colors[i])
            legend_lines.append(l)
    fig.legend(legend_lines,'AutoUpdate','off',labels=legend,loc="center right", fontsize=9)

    # plot dashed lines at resonance
    if len(w_resonance)>0:
        for plot_num in range(0,2):
            for i in range(0,len(w_resonance)):
                ax[plot_num].plot([w_resonance[i], w_resonance[i]],[ymin, ymax_pend],linestyle=':',color='k')

    # plot data
    if len(xdata)>0:
        for i in range(0,len(xdata)):
            ax[0].plot(w,xdata[i],linestyle=linestyles[i],color=colors[i])
        for i in range(0,len(ydata)):
            ax[1].plot(w,ydata[i],linestyle=linestyles[i],color=colors[i])

    # set plot parameters
    for plot_num in range(0,2):
        ax[plot_num].set_xscale('log')
        ax[plot_num].set_yscale('log')
        ax[plot_num].set_xlim((w[1],w[len(w)-1]))
        ax[plot_num].set_ylim(ymin,ymax)
        ax[plot_num].grid(True, color="#E0E0E0")
        for label in (ax[plot_num].get_xticklabels() + ax[plot_num].get_yticklabels()):
            label.set_fontsize(8)

    # create titles
    ax[0].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax[0].transAxes, fontsize=10)
    ax[0].text(.05,.93,'A.',horizontalalignment='left',transform=ax[0].transAxes, fontsize=10, fontweight='bold')
    ax[1].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax[1].transAxes, fontsize=10)
    ax[1].text(.05,.93,'B.',horizontalalignment='left',transform=ax[1].transAxes, fontsize=10, fontweight='bold')
    fig.text(0.5, 0.91,title, ha='center', fontsize=10, fontweight='bold')
    fig.text(0.5, 0.01,xlabel, ha='center', fontsize=10)
    fig.text(0.065, 0.5,ylabel, va='center', rotation='vertical', fontsize=10)
    fig.subplots_adjust(right=0.78)
    return [fig,ax]

def mag_spectrum(w,data,w_resonance,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend):

    """
    w - frequency values for the x-axis of plots
    data - data for plot in the form of a list of numpy arrays
    w_resonance - list of resonant frequencies of the ball for vertical lines
    title,xlabel,ylabel - strings for labeling the plot
    legend - list of string labels for the legend corresponding to linestyles,
            colors, and xdata/ydata (should be the same length)
    linestyles - list of string labels for line style
    colors - list of string hex colors
    ymin/ymax - plot y min and max
    ymax_pend - max y-value for pendulum line
    """
    # create plot
    figure_size = (6,3.55) # inches
    plt.figure(1, figsize=figure_size, dpi=150)
    fig, ax =plt.subplots(nrows=1, ncols=1, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.05)
    legend_lines = []

    # plot lines in order of legend
    if len(w_resonance)>0:
        l = ax.plot([w_resonance[0], w_resonance[0]],[ymin, ymax_pend],linestyle=':',color='k')
        legend_lines.append(l)
    if len(data)>0:
        for i in range(0,len(data)):
            l = ax.plot(w,data[i],linestyle=linestyles[i],color=colors[i])
            legend_lines.append(l)

    # plot dashed lines at resonance
    if len(w_resonance)>0:
        for i in range(0,len(w_resonance)):
            ax.plot([w_resonance[i], w_resonance[i]],[ymin, ymax_pend],linestyle=':',color='k')

    # plot data
    if len(data)>0:
        for i in range(0,len(data)):
            ax.plot(w,data[i],linestyle=linestyles[i],color=colors[i])

    # set plot parameters
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlim((10^-1,w[len(w)-1]))
    # ax.set_ylim(ymin,ymax)
    ax.set_xlim((w[1],w[len(w)-1]))
    ax.set_ylim(ymin,ymax)
    ax.grid(True, color="#E0E0E0")
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(8)

    # create titles
    fig.text(0.5, 0.91,title, ha='center', fontsize=10, fontweight='bold')
    fig.text(0.5, 0.01,xlabel, ha='center', fontsize=10)
    fig.text(0.065, 0.5,ylabel, va='center', rotation='vertical', fontsize=10)
    fig.legend(legend_lines,labels=legend,loc="center right", fontsize=9)
    fig.subplots_adjust(right=0.75)
    return [fig,ax]
