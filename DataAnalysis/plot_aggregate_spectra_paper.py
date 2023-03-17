import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils.plot_freq_utils import *
from utils.perform_transform import *


file_spectrum = "stroke-freq-spectrum.csv"
df_header = pd.read_csv(file_spectrum,header=None).to_numpy()[0]
df = pd.read_csv(file_spectrum)
w = df_header[5:]

mean_e_at_resonance = np.zeros((3,3))
std_e_at_resonance = np.zeros((3,3))

legend = ['Task Frequency','Nonparetic Arm\n(N=10)','Mild Stroke\n(N=3)','Moderate-Severe Stroke\n(N=7)']
title = 'Stroke Impairs High-Frequency Motion Bandwidth'
xlabel = 'Frequency (Hz)'
ylabel = 'Normalized Frequency Amplitude'
colors = ['#601a4aff', '#ee442fff','#63acbeff']
alphas = [.5,.5,.5]

cutoff = 5  # desired cutoff frequency of the filter, Hz
DT = 0.05
Fs = 1/DT
window = .3000000000001
FMA_cutoff = 35

# create plot
figure_size = (6,2.55) # inches
fig, ax =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=300)
fig.subplots_adjust(hspace=0.05)
fig.subplots_adjust(wspace=0.05)
legend_lines = []


frequencylabels = ['1Hz','1.5Hz','2.5Hz']
freq_pendulum = [1,1.5,2.5]
ymin = 0
ymax_pend = 1.35
ymax = 1.5
for freq in range(len(frequencylabels)):
    df_freq = df[df['BallFreq']==frequencylabels[freq]]
    ax[freq].plot([freq_pendulum[freq], freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')

    # line for nonparetic arm
    df_nonparetic = df_freq[df_freq['Arm']=='nonparetic']
    list = []
    list_ee = []
    subjects = df_nonparetic.Subject.unique()
    for sub in subjects:
        df_sub = df_nonparetic[df_nonparetic['Subject']==sub].to_numpy()[:,5:]
        A_mag = np.mean(df_sub,axis=0)
        list.append(A_mag)
        list_ee.append(find_energy_at_resonance(w,A_mag,freq_pendulum[freq],window))
    list = np.array(list,dtype=float)
    mean_e_at_resonance[0,freq] = np.mean(list_ee)
    std_e_at_resonance[0,freq] = np.std(list_ee)/ len(list_ee)
    y = np.mean(list,axis=0)
    error = np.std(list,axis=0)/np.sqrt(list.shape[0])
    # y = butter_lowpass_filter(y,cutoff,Fs)
    # error = butter_lowpass_filter(error,cutoff,Fs)
    ax[freq].plot(w, y, color=colors[0])
    ax[freq].fill_between(w.tolist(), y-error, y+error, color=colors[0],alpha=alphas[0],lw=0)


    # line for mild participants
    df_mild = df_freq[(df_freq['Arm']=='paretic') & (df_freq['FMA']>=FMA_cutoff)]
    list = []
    subjects = df_mild.Subject.unique()
    for sub in subjects:
        df_sub = df_mild[df_mild['Subject']==sub].to_numpy()[:,5:]
        A_mag = np.mean(df_sub,axis=0)
        list.append(A_mag)
        list_ee.append(find_energy_at_resonance(w,A_mag,freq_pendulum[freq],window))
    list = np.array(list,dtype=float)
    mean_e_at_resonance[1,freq] = np.mean(list_ee)
    std_e_at_resonance[1,freq] = np.std(list_ee)/ len(list_ee)
    y = np.mean(list,axis=0)
    error = np.std(list,axis=0)/np.sqrt(list.shape[0])
    # y = butter_lowpass_filter(y,cutoff,Fs)
    # error = butter_lowpass_filter(error,cutoff,Fs)
    ax[freq].plot(w, y, color=colors[1])
    ax[freq].fill_between(w.tolist(), y-error, y+error, color=colors[1],alpha=alphas[1],lw=0)

    # line for modsevere participants
    df_modsevere = df_freq[(df_freq['Arm']=='paretic') & (df_freq['FMA']<FMA_cutoff)]
    list = []
    subjects = df_modsevere.Subject.unique()
    for sub in subjects:
        df_sub = df_modsevere[df_modsevere['Subject']==sub].to_numpy()[:,5:]
        A_mag = np.mean(df_sub,axis=0)
        list.append(A_mag)
        list_ee.append(find_energy_at_resonance(w,A_mag,freq_pendulum[freq],window))
    list = np.array(list,dtype=float)
    mean_e_at_resonance[2,freq] = np.mean(list_ee)
    std_e_at_resonance[2,freq] = np.std(list_ee)/ len(list_ee)
    y = np.mean(list,axis=0)
    error = np.std(list,axis=0)/np.sqrt(list.shape[0])
    # y = butter_lowpass_filter(y,cutoff,Fs)
    # error = butter_lowpass_filter(error,cutoff,Fs)
    ax[freq].plot(w, y, color=colors[2])
    ax[freq].fill_between(w.tolist(), y-error, y+error, color=colors[2],alpha=alphas[2],lw=0)

    ax[freq].set_xlim((w[1],w[len(w)-1]))
    # ax.set_xlim((10^-1,w[len(w)-1]))
    ax[freq].set_ylim(ymin,ymax)
    # ax.set_xlim((w[1],w[len(w)-1]))
    # ax.set_ylim(ymin,ymax)

    # ax[freq].grid(True, color="#E0E0E0")
    for label in (ax[freq].get_xticklabels() + ax[freq].get_yticklabels()):
        label.set_fontsize(8)


# create titles
fig.text(0.5, 0.91,title, ha='center', fontsize=9, fontweight='bold')
fig.text(0.5, 0.01,xlabel, ha='center', fontsize=9)
fig.text(0.06, 0.5,ylabel, va='center', rotation='vertical', fontsize=9)
if len(legend)>0:
    fig.legend(legend_lines,labels=legend,loc="center right", fontsize=9)
    fig.subplots_adjust(right=0.69)

fig.savefig('Plots_stroke/paper_spectrum.pdf')
fig.savefig('Plots_stroke/paper_spectrum.png')


# Line plots for raw data
figure_size = (2,2.55)
fig_scat, ax_scat = plt.subplots(figsize=figure_size,dpi=300)
xlabel = 'Resonant Frequency of Task'
ylabel = 'Fraction of Total Energy'
title = ''
conditions = ['Nonparetic Arm\n(N=10)','Mild Stroke\n(N=3)','Moderate-Severe Stroke\n(N=7)']

# place grid in back
ax_scat.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
               alpha=0.5)
ax_scat.set_axisbelow(True)

# Add titles and labels
plt.xlabel(xlabel,fontname="sans-serif", fontsize=9)
plt.ylabel(ylabel,fontname="sans-serif", fontsize=9)
plt.title(title,fontname="sans-serif", fontsize=9,fontweight='bold')
for label in (ax_scat.get_yticklabels()):
    label.set_fontsize(8)

# x-ticks x-axis
plt.xticks(freq_pendulum, frequencylabels, fontname="sans-serif", fontsize=9)
for tick in ax_scat.get_xticklabels():
    tick.set_rotation(0)

p = np.zeros(3, dtype=object)
for group in range(3):
    # for freq in range(0,3):
        # i = arm *2 + SL
        # print(arm,SL,data_mean)
        # p[i] = ax_scat.errorbar(ind,data_mean,yerr=data_std,color=colors[arm],ls=linestyles[SL],marker="o",ecolor=colors[arm],capsize=5)
    ax_scat.errorbar(freq_pendulum,mean_e_at_resonance[group,:],yerr=std_e_at_resonance[group,:],
                        color=colors[group],marker="o",ecolor=colors[group],capsize=3,ms=3)

# L = ax_scat.legend(p, conditions, fontsize=9,loc='upper right')
# L = fig.legend(p, conditions, fontsize=9,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
# plt.setp(L.texts, family='sans-serif')
fig_scat.subplots_adjust(bottom=0.15,left=.2)
fig_scat.savefig('Plots_stroke/paper_e_at_res.pdf')
fig_scat.savefig('Plots_stroke/paper_e_at_res.png')
# plt.close(fig_scat)




plt.show()
