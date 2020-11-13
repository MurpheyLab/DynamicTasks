import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
from preform_transform import calculate_amplitude
from make_boxplot import make_boxplot, add_stats
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import csv

# Edit these variables before running
DIR = "Z:" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
# subject_num = 4
DT = 0.01
Fs = 1/DT
freq_pendulum = (1/(2*np.pi))*np.sqrt(9.81/0.07) #resonant frequency of pendulum
individual_plots = 0

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.1 # set the resolution for the frequcny bins
dt_measure = .5 # seconds for each point
num_i_measure = int(dt_measure*Fs)

# Label factors as strings for ezANOVA analysis
supportlevels = ['0%Max','20%Max','50%Max']
tasks = ['Task1','Task2','Task3','Task4','Task5']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
color_array = ['#fc0303','#fc9403','#d2d904','#04d92b','#04d9d2','#040fd9','#8b04d9']
group_name = ['Controls','Stroke']
group_name = ['Stroke']
linestyles = ['--','-']
markerstyles = ['o','D']

workspace_limits = np.zeros((17,5))
# myfile << armWeight << ',' << ymin << ',' << ymax << ',' << xmin << ',' << xmax << '\n';
workspace_limits[1,:] = [24.9512,-0.0447694,0.219882,-0.0522197,-0.033374]
workspace_limits[2,:] = [28.7598,-0.111891,0.206901,-0.0173099,0.00216103]
workspace_limits[3,:] = [46.3867,-0.0784368,0.19,-0.04,-0.0103121]
workspace_limits[4,:] = [42.749,-0.0319042,0.143313,-0.0852239,-0.0443593]
workspace_limits[5,:] = [22.6807,-0.0104683,0.232982,0.00131819,0.0204978]
workspace_limits[6,:] = [22.3145,-0.0601319,0.0241804,-0.04,-0.01]
workspace_limits[11,:] = [24.7803,-0.157845,0.202717,-0.0471689,0.068145]
workspace_limits[12,:] = [23.7061,-0.0708732,0.219661,-0.0760239,0.0782709]
workspace_limits[13,:] = [22.9492,-0.0959921,0.208944,-0.0747804,0.0768332]
workspace_limits[14,:] = [42.4072,-0.150891,0.205025,-0.0481668,0.119187]
workspace_limits[15,:] = [21.4355,-0.113227,0.174447,-0.118506,0.0248991]
workspace_limits[16,:] = [25.8057,-0.110519,0.135236,-0.073724,0.013615]

# Initialize arrays for calculating effect size
# row = subject
# col = in order [stroke0%close; stroke0%far; control0%close; control0%far]
data_closevsfar = np.zeros((5,4))

subs_test = [1,2,5]
minsub = 1
if individual_plots==1:
    numsub = 5
else:
    numsub = minsub
for group in range(len(group_name)):
    for subject_num in range(minsub,numsub+1):
        # Sets the correct subfile depending on which group the subject is in
        if (group_name[group] == 'Stroke'):
            subname = "S0" + str(subject_num)
        else:
            subject_num += 10
            subname = "S" + str(subject_num)
        subfile = "/OutputData_S" + str(subject_num)

        # Set up subject parameters
        mx = workspace_limits[subject_num,4] - workspace_limits[subject_num,3]
        my = workspace_limits[subject_num,2] - workspace_limits[subject_num,1]

        cutoff = .7
        if individual_plots==1:
            figure_size = (6,3.55) # inches
            fig, ax = plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
            figure_size = (5,3.55) # inches
            # plt.figure(subject_num,figsize=figure_size, dpi=150)
            fig_all, ax_all = plt.subplots(nrows=1, ncols=1, sharey='row', squeeze=True, figsize=figure_size, dpi=150)

        for j in range(0,3): #iterate through support levels
            points_all = np.zeros((1,2))
            points_close = []
            points_far = []
            for i in range(1,6): #iterate through tasks
                for k in range(1,46): #iterate through trials
                    try:
                        # Open the trial files
                        trialfile = DIR+subname+subfile+"_Task"+str(i)+"_SL"+str(j)+"_Trial"+ str(k)+".csv"
                        data = genfromtxt(trialfile,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names
                        print(trialfile)
                        # Determine the indecies for the start and end of trial
                        delay = data[-1,0]-20.0 # Gets the final time to calculate the delay in data collection
                        delayInd = int(delay/DT)
                        score = data[-1,7]
                        # gets end of trial by moving backwards until the score is less than 20
                        trial_complete_index = int(data.shape[0]-1)
                        if score==20:
                            while data[trial_complete_index,7]==20:
                                trial_complete_index-=1

                        # Get points for workspace plots
                        i_start = delayInd
                        i_end = i_start+num_i_measure
                        while i_end<trial_complete_index:
                            # Use signal directly 1-x position 2- y position 8- x force 9- y force
                            # y_Fx = data[i_start:(i_end-1),8]
                            # 1/2/3/4-far 5/-close
                            y_Fx = np.sqrt(np.square(data[i_start:(i_end-1),8])+np.square(data[i_start:(i_end-1),9]))
                            [A_Fx,frq] = calculate_amplitude(y_Fx,Fs)
                            # print(frq)
                            resonant_peak_Fx = 0
                            low_freq_energy_Fx = 0
                            high_freq_energy_Fx = 0
                            freq_list = []
                            for w_i in range(0,len(frq)):
                                # Criteria for high frequency 2.75<frq<10Hz -- mulitply by 2 for freq magnitude
                                # if (frq[w_i]>2.75*2) and (frq[w_i]<10*2):
                                # if (frq[w_i]>freq_pendulum*2-1) and (frq[w_i]<freq_pendulum*2+1):
                                    # print(frq[w_i])
                                # if (frq[w_i]>freq_pendulum+1) and (frq[w_i]<10):
                                # if (frq[w_i]>1*2):# and (frq[w_i]<10*2):
                                if (frq[w_i]>freq_pendulum*2+1) and (frq[w_i]<10*2):
                                    high_freq_energy_Fx += A_Fx[w_i]
                                    freq_list.append(A_Fx[w_i])

                            dw = frq[1]-frq[0]
                            # high_freq_energy_Fx = np.log(np.sum(np.square(freq_list))*dw)
                            high_freq_energy_Fx = np.sum(np.square(freq_list))*dw

                            x_loc_mean = np.mean(data[i_start:(i_end-1),1])
                            y_loc_mean = np.mean(data[i_start:(i_end-1),2])
                            y_loc = 1-((x_loc_mean-workspace_limits[subject_num,3])/mx)
                            x_loc = (y_loc_mean-workspace_limits[subject_num,1])/my

                            if high_freq_energy_Fx == 0:
                                print('no high freq energy')

                            if y_loc>0.275 and y_loc<1:
                            # if x_loc>0 and x_loc<1:
                                if y_loc<cutoff:
                                # if x_loc<cutoff:
                                    points_close.append(high_freq_energy_Fx)
                                else:
                                    points_far.append(high_freq_energy_Fx)
                                points_all=np.append(points_all,[[y_loc,high_freq_energy_Fx]],axis=0)
                                # points_all=np.append(points_all,[[x_loc,high_freq_energy_Fx]],axis=0)
                            i_start = i_end # i_start+int(num_i_measure*0.25)
                            i_end = i_start+num_i_measure

                    except:
                        pass

            labels = ['Close','Far']
            box_colors = ['#63ACBE','#EE442F']
            box_alpha = [1,1]
            # indent for each subject figure
            if j==0:
                if subject_num==minsub:
                    p_close_0 = points_close
                    p_far_0 = points_far
                # elif subject_num==subs_test[1]:# or subject_num==subs_test[2]:
                # p_close_0 = np.append(p_close_0,points_close)
                # p_far_0 = np.append(p_far_0,points_far)
                if individual_plots==1:
                    axnum = 0
                    ax[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
                    ax[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
                    data = [points_close,points_far]
                    make_boxplot(fig,ax[axnum],data,labels,box_colors,box_alpha)
                    [tstat,pval] = stats.ttest_ind(data[0],data[1], equal_var = False)
                    sig_matrix = np.array([0,1,pval])
                    add_stats(ax[axnum],data,sig_matrix)
            elif j==1:
                if subject_num==minsub:
                    p_close_20 = points_close
                    p_far_20 = points_far
                # elif subject_num==subs_test[1]:# or subject_num==subs_test[2]:
                # p_close_20 = np.append(p_close_20,points_close)
                # p_far_20 = np.append(p_far_20,points_far)

                if individual_plots==1:
                    axnum = 1
                    ax[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
                    ax[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
                    data = [points_close,points_far]
                    make_boxplot(fig,ax[axnum],data,labels,box_colors,box_alpha)
                    [tstat,pval] = stats.ttest_ind(data[0],data[1], equal_var = False)
                    sig_matrix = np.array([0,1,pval])
                    add_stats(ax[axnum],data,sig_matrix)
            elif j==2:
                if subject_num==minsub:
                    p_close_50 = points_close
                    p_far_50 = points_far
                # elif subject_num==subs_test[1]:# or subject_num==subs_test[2]:
                # p_close_50 = np.append(p_close_50,points_close)
                # p_far_50 = np.append(p_far_50,points_far)
                if individual_plots==1:
                    axnum = 2
                    ax[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
                    ax[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
                    data = [points_close,points_far]
                    make_boxplot(fig,ax[axnum],data,labels,box_colors,box_alpha)
                    [tstat,pval] = stats.ttest_ind(data[0],data[1], equal_var = False)
                    sig_matrix = np.array([0,1,pval])
                    add_stats(ax[axnum],data,sig_matrix)

            if individual_plots==1:
                # Evalutate for a range of cutoff values
                cutoffvals = np.linspace(0.3, 0.95, 30)
                numpoints = points_all.shape[0]
                pvals = []
                for i in range(len(cutoffvals)):
                    p_close = []
                    p_far = []
                    for row in range(numpoints):
                        if points_all[row,0]<cutoffvals[i]:
                            p_close.append(points_all[row,1])
                        else:
                            p_far.append(points_all[row,1])
                    [tstat,pval] = stats.ttest_ind(p_close,p_far, equal_var = False)
                    if pval<.00015:
                        pval = .00015
                    pvals.append(pval)
                # plt.figure(subject_num)
                ax_all.plot(cutoffvals,pvals)

        sub_index = subject_num-1
        close_far_index = 0
        if (group_name[group] == 'Controls'):
            sub_index -= 10
            close_far_index = 2
        data_closevsfar[sub_index,close_far_index] = np.mean(p_close_0)
        data_closevsfar[sub_index,close_far_index+1] = np.mean(p_far_0)

        if individual_plots==1:
            fig.subplots_adjust(hspace=0.05)
            fig.subplots_adjust(wspace=0.05)
            for axnum in range(3):
                # ax[axnum].grid(True)
                # ax[axnum].set_xlim([0,1])
                [ymin,ymax]=ax[axnum].get_ylim()
                min_dist = (ymax-ymin)/50.0
                ax[axnum].set_ylim([ymin,ymax+min_dist*2])
                for label in (ax[axnum].get_yticklabels()):
                    label.set_fontsize(9)
                for label in (ax[axnum].get_xticklabels()):
                    label.set_fontsize(10)
            ax_all.plot([.3,.95],[.05,.05],'-k')
            fig.text(0.5, 0.91, 'High Frequency Content vs. Workspace Location For Subject '+str(subject_num), ha='center', fontsize=10, fontweight='bold')
            fig.text(0.5, 0.01, 'Position of Hand Relative to the Body', ha='center', fontsize=10)
            fig.text(0.06, 0.5, 'High Frequency Content', va='center', rotation='vertical', fontsize=10)
            fig.savefig('sub'+str(subject_num)+'workspaceloc.png')


            ax_all.grid(True)
            ax_all.set_yscale('log')
            ax_all.set_ylim(.0001,1.)
            for label in (ax_all.get_yticklabels()):
                label.set_fontsize(9)
            for label in (ax_all.get_xticklabels()):
                label.set_fontsize(9)
            fig_all.text(0.5, 0.91, 'p-values for a range of cutoff choices for subject '+str(subject_num), ha='center', fontsize=10, fontweight='bold')
            # fig_all.title('p-values for a range of cutoff choices for subject '+str(subject_num), ha='center', fontsize=10, fontweight='bold')
            fig_all.text(0.5, .01, 'Cut-off value', ha='center', fontsize=10)
            # fig_all.ylabel('P-value', ha='center', fontsize=10)
            fig_all.text(0.01, 0.5, 'P-value', va='center', rotation='vertical', fontsize=10)
            ax_all.legend(supportlevels,loc="upper center", fontsize=9)
            fig_all.savefig('sub'+str(subject_num)+'pvals.png')

print(data_closevsfar)

# for repeated measures
mean1 = np.mean(data_closevsfar[0])
mean2 = np.mean(data_closevsfar[1])
std1 = np.std(data_closevsfar[0])
std2 = np.std(data_closevsfar[1])
std_pooled = np.sqrt((std1**2+std2**2)/2)
print(mean1,mean2,std1,std2,std_pooled)
d = abs((mean1-mean2)/std_pooled)
print('d for repeated measures',d)

# for stroke vs control

stroke_all = np.array(data_closevsfar[1]-data_closevsfar[0])
controls_all = np.array(data_closevsfar[3]-data_closevsfar[2])
mean1 = np.mean(stroke_all.flatten())
mean2 = np.mean(controls_all.flatten())
std1 = np.std(stroke_all.flatten())
std2 = np.std(controls_all.flatten())
std_pooled = np.sqrt((std1**2+std2**2)/2)
d = abs((mean1-mean2)/std_pooled)
# print(mean1,mean2,std1,std2)
print('d for stroke vs. controls',d)


figure_size = (5,3.25) # sets the size of the figure in inches
fig, ax1 = plt.subplots(figsize=figure_size, dpi=150)

ymin = -0.03
ymax = 0.2 #35 # 1.1
ax1.set_ylim((ymin,ymax))
ax1.add_patch(Rectangle((.5,ymin), 2, ymax-ymin, facecolor='#601A4A', alpha=.15))
ax1.add_patch(Rectangle((2.5,ymin), 2, ymax-ymin, facecolor='#EE442F', alpha=.15))
ax1.add_patch(Rectangle((4.5,ymin), 2, ymax-ymin, facecolor='#63ACBE', alpha=.15))

data = [p_close_0,p_far_0,p_close_20,p_far_20,p_close_50,p_far_50]
# labels = ['0%-close','0%-far','20%-close','20%-far','50%-close','50%-far']
labels = ['Proximal','Distal','Proximal','Distal','Proximal','Distal']
box_colors = ['#601A4A', '#601A4A','#EE442F','#EE442F', '#63ACBE', '#63ACBE']
box_colors = ['#909090', '#E8E8E8','#909090','#E8E8E8', '#909090', '#E8E8E8']
box_alpha = [None,0.3,None,0.3,None,0.3]
box_alpha = [None,None,None,None,None,None]
# box_alpha = [0.4,0.1,0.4,0.1,0.4,0.1]
make_boxplot(fig,ax1,data,labels,box_colors,box_alpha)

combos = [[0,2],[2,4],[0,4],[1,3],[3,5],[1,5],[0,1],[2,3],[4,5]]
combos = [[0,1],[2,3],[4,5]]
sig_matrix = np.zeros((len(combos),3))
for i in range(len(combos)):
    [tstat,pval] = stats.ttest_ind(data[combos[i][0]],data[combos[i][1]], equal_var = False)
    # [tstat,pval] = stats.ttest_rel(data[combos[i][0]],data[combos[i][1]])
    pval /= 2
    sig_matrix[i,:] = [combos[i][0],combos[i][1],pval]
print(sig_matrix)
add_stats(ax1,data,sig_matrix)

ax1.set_title('High-Frequency Content Over Workspace for Subject 1', fontweight='bold',fontname="Arial", fontsize=11, color='black')#, fontweight='bold')
# ax1.set_title('Frequency Content Around Resonance for Subject '+str(minsub), fontweight='bold',fontname="Arial", fontsize=11, color='black')#, fontweight='bold')
# ax1.set_title('Frequency Content Above Resonance for Subject '+str(minsub), fontweight='bold',fontname="Arial", fontsize=11, color='black')#, fontweight='bold')
# ax1.set_xlabel('Loading Level (% Max) - Distance From Body',fontname="Arial", fontsize=11)
fig.subplots_adjust(bottom=0.2)
# ax1.set_ylabel('High-Frequency Content',fontname="Arial", fontsize=11, color='black')
ax1.set_ylabel('Frequency Content',fontname="Arial", fontsize=11, color='black')
for label in ax1.get_xticklabels():
    label.set_fontsize(10)
    label.set_color('black')
for label in ax1.get_yticklabels():
    label.set_fontsize(9)
    label.set_color('black')
box_coords = []

arrow_height = -.055 #-.16
text_height = -.18
ax1.annotate('', xy=(.5,arrow_height),xytext=(2.5,arrow_height),                     #draws an arrow from one set of coordinates to the other
            arrowprops=dict(arrowstyle='<|-|>',facecolor='black'),   #sets style of arrow and colour
            annotation_clip=False)                               #This enables the arrow to be outside of the plot
ax1.text(.166666, text_height,'0% Added',horizontalalignment='center',transform=ax1.transAxes, fontname="Arial", fontsize=10)
ax1.annotate('', xy=(2.5,arrow_height),xytext=(4.5,arrow_height),                     #draws an arrow from one set of coordinates to the other
            arrowprops=dict(arrowstyle='<|-|>',facecolor='black'),   #sets style of arrow and colour
            annotation_clip=False)                               #This enables the arrow to be outside of the plot
ax1.text(.5, text_height,'20% Added',horizontalalignment='center',transform=ax1.transAxes, fontname="Arial", fontsize=10)
ax1.annotate('', xy=(4.5,arrow_height),xytext=(6.5,arrow_height),                     #draws an arrow from one set of coordinates to the other
            arrowprops=dict(arrowstyle='<|-|>',facecolor='black'),   #sets style of arrow and colour
            annotation_clip=False)                               #This enables the arrow to be outside of the plot
ax1.text(.8333333, text_height,'50% Added',horizontalalignment='center',transform=ax1.transAxes, fontname="Arial", fontsize=10)
ax1.text(0.0, text_height,'Loading:',horizontalalignment='right',transform=ax1.transAxes, fontname="Arial", fontsize=10, fontweight='bold')
ax1.text(0.0, -.09,'Location:',horizontalalignment='right',transform=ax1.transAxes, fontname="Arial", fontsize=10, fontweight='bold')
fig.subplots_adjust(left=0.2)

fig.savefig('sub'+str(minsub)+'_workspace.png')

plt.show()
