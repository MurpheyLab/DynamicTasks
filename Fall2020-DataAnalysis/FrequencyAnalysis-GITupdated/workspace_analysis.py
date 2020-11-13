import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
from preform_transform import calculate_amplitude
import matplotlib.pyplot as plt
import csv

# Edit these variables before running
DIR = "Z:" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
number_of_subjects = 6 # in each group, controls and stroke
DT = 0.01
Fs = 1/DT
freq_pendulum = (1/(2*np.pi))*np.sqrt(9.81/0.07) #resonant frequency of pendulum

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.1 # set the resolution for the frequcny bins
dt_measure = 0.6 # seconds for each point
num_i_measure = int(dt_measure*Fs)

# Label factors as strings for ezANOVA analysis
supportlevels = ['0%Max','20%Max','50%Max']
tasks = ['Task1','Task2','Task3','Task4','Task5']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
color_array = ['#fc0303','#fc9403','#d2d904','#04d92b','#04d9d2','#040fd9','#8b04d9']
group_name = ['Controls','Stroke']
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


fignum = 1

for group in range(0,2):#len(group_name)):
    # pointsall = np.zeros((1,3))

    figure_size = (8,3.55) # inches
    fig, ax =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig_time, ax_time =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig_timecorr, ax_timecorr =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig_timecorry, ax_timecorry =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    fig_yposcorr, ax_yposcorr =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
    for j in range(0,3): #iterate through support levels
        points = np.zeros((1,3))
        points_time = np.zeros((1,3))
        points_timecorr = np.zeros((1,3))
        points_yposcorr = np.zeros((1,3))
        for subject_num in range(1,1+number_of_subjects): #iterate though subjects

            # Sets the correct subfile depending on which group the subject is in
            if (group_name[group] == 'Stroke'):
                subname = "S0" + str(subject_num)
            else:
                subject_num += 10
                subname = "S" + str(subject_num)
            subfile = "/OutputData_S" + str(subject_num)

            mx = workspace_limits[subject_num,4] - workspace_limits[subject_num,3]
            my = workspace_limits[subject_num,2] - workspace_limits[subject_num,1]

        # figure_size = (8,3.55) # inches
        # fignum += 1
        # fig, ax =plt.subplots(nrows=1, ncols=3, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        # for j in range(0,3): #iterate through support levels
        #     points = np.zeros((1,3))
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

                        # Get points for time plots
                        for iii in range(delayInd,trial_complete_index,30):
                            index = 1. * (iii - delayInd)/(trial_complete_index - delayInd)
                            y_loc = 1-(data[iii,1]-workspace_limits[subject_num,3])/mx
                            x_loc = (data[iii,2]-workspace_limits[subject_num,1])/my
                            if np.max(points_time)>0:
                                if y_loc<1:
                                    points_time = np.append(points_time, [[x_loc,y_loc,index]], axis=0)
                            else:
                                if y_loc<1:
                                    points_time = [[x_loc,y_loc,index]]

                        # Get points for workspace plots
                        i_start = delayInd
                        i_end = i_start+num_i_measure
                        while i_end<trial_complete_index:
                            # Use signal directly 1-x position 2- y position 8- x force 9- y force
                            # y_Fx = data[i_start:(i_end-1),8]
                            y_Fx = np.sqrt(np.square(data[i_start:(i_end-1),8])+np.square(data[i_start:(i_end-1),9]))
                            [A_Fx,frq] = calculate_amplitude(y_Fx,Fs)
                            resonant_peak_Fx = 0
                            low_freq_energy_Fx = 0
                            high_freq_energy_Fx = 0
                            for w_i in range(0,len(frq)):
                                if (frq[w_i]<2):
                                    low_freq_energy_Fx += A_Fx[w_i]
                                    last_under_1_Hz = w_i
                                else:
                                    high_freq_energy_Fx += A_Fx[w_i]
                            ratio_energy_Fx = high_freq_energy_Fx/low_freq_energy_Fx

                            # print(frq)
                            # print(last_under_1_Hz)
                            dw = frq[1]-frq[0]
                            high_freq_energy_Fx = np.sum(np.square(A_Fx[0:last_under_1_Hz]))*dw
                            # print(high_freq_energy_Fx)

                            if True: #data[i_end-1,1]<data[i_start,1]: # if the person is reaching out
                                x_loc_mean = np.mean(data[i_start:(i_end-1),1])
                                y_loc_mean = np.mean(data[i_start:(i_end-1),2])

                                y_loc = 1-(x_loc_mean-workspace_limits[subject_num,3])/mx
                                x_loc = (y_loc_mean-workspace_limits[subject_num,1])/my
                                if high_freq_energy_Fx == 0:
                                    print('no high freq energy')
                                if np.max(points)>0:
                                    if y_loc<1:
                                        points = np.append(points, [[x_loc,y_loc,high_freq_energy_Fx]], axis=0)#np.append(points,zeros((1,3)),axis=0)
                                else:
                                    if y_loc<1:
                                        points = [[x_loc,y_loc,high_freq_energy_Fx]]

                                if np.max(points_yposcorr)>0:
                                    if y_loc>0.3 and y_loc<1:
                                        points_yposcorr = np.append(points_yposcorr, [[x_loc,y_loc,high_freq_energy_Fx]], axis=0)#np.append(points,zeros((1,3)),axis=0)
                                else:
                                    if y_loc>0.3 and y_loc<1:
                                        points_yposcorr = [[x_loc,y_loc,high_freq_energy_Fx]]

                                index = 1. * (i_start+((i_end-i_start)/2) - delayInd)/(trial_complete_index - delayInd)
                                if np.max(points_timecorr)>0:
                                    if y_loc<1:
                                        points_timecorr = np.append(points_timecorr, [[index,high_freq_energy_Fx,y_loc]], axis=0)#np.append(points,zeros((1,3)),axis=0)
                                else:
                                    if y_loc<1:
                                        points_timecorr = [[index,high_freq_energy_Fx,y_loc]]

                            i_start = i_start+int(num_i_measure*0.25)
                            i_end = i_start+num_i_measure

                    except:
                        pass

        # indent for each subject figure
        if j==0:
            axnum = 0
            # High Freq Content over workspace
            vmin_val=np.min(points[:,2])
            vmax_val=np.max(points[:,2])
            ax[axnum].scatter(points[:,0],points[:,1],s=0.5,c=points[:,2],cmap='jet',vmin=vmin_val,vmax=vmax_val)
            ax[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
            ax[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
            #Time over workspace
            vmin_val_time=np.min(points_time[:,2])
            vmax_val_time=np.max(points_time[:,2])
            ax_time[axnum].scatter(points_time[:,0],points_time[:,1],s=0.5,c=points_time[:,2],cmap='jet',vmin=vmin_val_time,vmax=vmax_val_time)
            ax_time[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10)
            ax_time[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10, fontweight='bold')

            #Time correlation with freq. content over workspace
            ax_timecorr[axnum].scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
            ax_timecorr[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorr[axnum].plot(xp,yp,color='#000000')
            ax_timecorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)

            #Time correlation with y-position over workspace
            ax_timecorry[axnum].scatter(points_timecorr[:,0],points_timecorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,2])
            ax_timecorry[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorry[axnum].plot(xp,yp,color='#000000')
            ax_timecorry[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)

            #Ypos correlation with high freq. content over workspace
            ax_yposcorr[axnum].scatter(points_yposcorr[:,1],points_yposcorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_yposcorr[:,1],points_yposcorr[:,2])
            ax_yposcorr[axnum].text(.15,.93,'0% Max',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10)
            ax_yposcorr[axnum].text(.05,.93,'A.',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(.3, 1, 10)
            yp = xp*slope + intercept
            ax_yposcorr[axnum].plot(xp,yp,color='#000000')
            ax_yposcorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)
            ax_yposcorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)

        elif j==1:
            axnum = 1
            # High Freq Content over workspace
            ax[axnum].scatter(points[:,0],points[:,1],s=0.5,c=points[:,2],cmap='jet',vmin=vmin_val,vmax=vmax_val)
            ax[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
            ax[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
            #Time over workspace
            vmin_val_time=np.min(points_time[:,2])
            vmax_val_time=np.max(points_time[:,2])
            ax_time[axnum].scatter(points_time[:,0],points_time[:,1],s=0.5,c=points_time[:,2],cmap='jet',vmin=vmin_val_time,vmax=vmax_val_time)
            ax_time[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10)
            ax_time[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10, fontweight='bold')

            #Time correlation with freq. content over workspace
            ax_timecorr[axnum].scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
            ax_timecorr[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorr[axnum].plot(xp,yp,color='#000000')
            ax_timecorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)

            #Time correlation with y-position over workspace
            ax_timecorry[axnum].scatter(points_timecorr[:,0],points_timecorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,2])
            ax_timecorry[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorry[axnum].plot(xp,yp,color='#000000')
            ax_timecorry[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)

            #Ypos correlation with high freq. content over workspace
            ax_yposcorr[axnum].scatter(points_yposcorr[:,1],points_yposcorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_yposcorr[:,1],points_yposcorr[:,2])
            ax_yposcorr[axnum].text(.15,.93,'20% Max',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10)
            ax_yposcorr[axnum].text(.05,.93,'B.',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(.3, 1, 10)
            yp = xp*slope + intercept
            ax_yposcorr[axnum].plot(xp,yp,color='#000000')
            ax_yposcorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)
            ax_yposcorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)

        elif j==2:
            axnum = 2
            # High Freq Content over workspace
            cb3 = ax[axnum].scatter(points[:,0],points[:,1],s=0.5,c=points[:,2],cmap='jet',vmin=vmin_val,vmax=vmax_val)
            ax[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10)
            ax[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax[axnum].transAxes, fontsize=10, fontweight='bold')
            #Time over workspace
            vmin_val_time=np.min(points_time[:,2])
            vmax_val_time=np.max(points_time[:,2])
            cb3_time = ax_time[axnum].scatter(points_time[:,0],points_time[:,1],s=0.5,c=points_time[:,2],cmap='jet',vmin=vmin_val_time,vmax=vmax_val_time)
            ax_time[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10)
            ax_time[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax_time[axnum].transAxes, fontsize=10, fontweight='bold')
            #Time correlation with freq. content over workspace
            ax_timecorr[axnum].scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
            ax_timecorr[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorr[axnum].plot(xp,yp,color='#000000')
            ax_timecorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)
            ax_timecorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorr[axnum].transAxes, fontsize=10)

            #Time correlation with y-position over workspace
            ax_timecorry[axnum].scatter(points_timecorr[:,0],points_timecorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,2])
            ax_timecorry[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(0, 1, 10)
            yp = xp*slope + intercept
            ax_timecorry[axnum].plot(xp,yp,color='#000000')
            ax_timecorry[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            ax_timecorry[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left',transform=ax_timecorry[axnum].transAxes, fontsize=10)
            #Ypos correlation with high freq. content over workspace
            ax_yposcorr[axnum].scatter(points_yposcorr[:,1],points_yposcorr[:,2],s=0.5,c='#1f77b4')
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_yposcorr[:,1],points_yposcorr[:,2])
            ax_yposcorr[axnum].text(.15,.93,'50% Max',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10)
            ax_yposcorr[axnum].text(.05,.93,'C.',horizontalalignment='left',transform=ax_yposcorr[axnum].transAxes, fontsize=10, fontweight='bold')
            xp = np.linspace(.3, 1, 10)
            yp = xp*slope + intercept
            ax_yposcorr[axnum].plot(xp,yp,color='#000000')
            ax_yposcorr[axnum].text(.05,.88,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)
            ax_yposcorr[axnum].text(.05,.83,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),transform=ax_yposcorr[axnum].transAxes,horizontalalignment='left', fontsize=10)


    # High Freq Content over workspace
    for axnum in range(3):
        # ax[axnum].grid(True)
        ax[axnum].set_xlim([0,1])
        ax[axnum].set_ylim([0,1.1])
        for label in (ax[axnum].get_xticklabels() + ax[axnum].get_yticklabels()):
            label.set_fontsize(8)
    # fig.text(0.5, 0.91, 'High Frequency Content in Force Magnitude For Sub '+str(subject_num), ha='center', fontsize=10, fontweight='bold')
    fig.text(0.5, 0.91, 'High Frequency Content in Force Magnitude For ' + str(group_name[group]), ha='center', fontsize=10, fontweight='bold')
    fig.text(0.5, 0.01, 'X-position (Side-to-side)', ha='center', fontsize=10)
    fig.text(0.065, 0.5, 'Y-position (Reaching Forward)', va='center', rotation='vertical', fontsize=10)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cb3, cax=cbar_ax)
    # fig.savefig('sub'+str(subject_num)+'_freqmag_workspace'+'.png')
    fig.savefig(str(group_name[group])+'_all_freqmag_workspace'+'.png')

    #Time over workspace
    for axnum_time in range(3):
        # ax[axnum].grid(True)
        ax_time[axnum_time].set_xlim([0,1])
        ax_time[axnum_time].set_ylim([0,1.1])
        for label in (ax_time[axnum_time].get_xticklabels() + ax_time[axnum_time].get_yticklabels()):
            label.set_fontsize(8)
    fig_time.text(0.5, 0.91, 'Position in Workspace Over For ' + str(group_name[group]), ha='center', fontsize=10, fontweight='bold')
    fig_time.text(0.5, 0.01, 'X-position (Side-to-side)', ha='center', fontsize=10)
    fig_time.text(0.065, 0.5, 'Y-position (Reaching Forward)', va='center', rotation='vertical', fontsize=10)
    fig_time.subplots_adjust(right=0.8)
    cbar_ax_time = fig_time.add_axes([0.85, 0.15, 0.05, 0.7])
    fig_time.colorbar(cb3_time, cax=cbar_ax_time)
    fig_time.savefig(str(group_name[group])+'_all_time_workspace'+'.png')

    #Time correlation with freq. content over workspace
    for axnum_timecorr in range(3):
        # ax[axnum].grid(True)
        for label in (ax_timecorr[axnum_timecorr].get_xticklabels() + ax_timecorr[axnum_timecorr].get_yticklabels()):
            label.set_fontsize(8)
    fig_timecorr.text(0.5, 0.91, 'Time vs. High Frequency Content For ' + str(group_name[group]), ha='center', fontsize=10, fontweight='bold')
    fig_timecorr.text(0.5, 0.01, 'Time (0-start 1-end)', ha='center', fontsize=10)
    fig_timecorr.text(0.065, 0.5, 'High Frequency Content', va='center', rotation='vertical', fontsize=10)
    fig_timecorr.savefig(str(group_name[group])+'_all_timecorr_workspace'+'.png')

    #Time correlation with y-position over workspace
    for axnum_timecorry in range(3):
        # ax[axnum].grid(True)
        ax_timecorry[axnum_timecorry].set_xlim([0,1])
        ax_timecorry[axnum_timecorry].set_ylim([0,1.2])
        for label in (ax_timecorry[axnum_timecorry].get_xticklabels() + ax_timecorry[axnum_timecorry].get_yticklabels()):
            label.set_fontsize(8)
    fig_timecorry.text(0.5, 0.91, 'Time vs. Y-position For ' + str(group_name[group]), ha='center', fontsize=10, fontweight='bold')
    fig_timecorry.text(0.5, 0.01, 'Time (0-start 1-end)', ha='center', fontsize=10)
    fig_timecorry.text(0.065, 0.5, 'Y-position (Reaching Forward)', va='center', rotation='vertical', fontsize=10)
    fig_timecorry.savefig(str(group_name[group])+'_all_timecorry_workspace'+'.png')

    #Ypos correlation with high freq. content over workspace
    for axnum_yposcorr in range(3):
        # ax[axnum].grid(True)
        for label in (ax_yposcorr[axnum_yposcorr].get_xticklabels() + ax_yposcorr[axnum_yposcorr].get_yticklabels()):
            label.set_fontsize(8)
    fig_yposcorr.text(0.5, 0.91, 'Y-position vs. High Frequency Content in Workspace For ' + str(group_name[group]), ha='center', fontsize=10, fontweight='bold')
    fig_yposcorr.text(0.5, 0.01, 'Y-position (Reaching Forward)', ha='center', fontsize=10)
    fig_yposcorr.text(0.065, 0.5, 'High Frequency Content', va='center', rotation='vertical', fontsize=10)
    fig_yposcorr.savefig(str(group_name[group])+'_all_yposfreqmagcorr_workspace'+'.png')


    # plt.figure(50+group*10)
    # plt.title('High Frequency Content in Force Magnitude vs. Y-position For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    # plt.scatter(pointsall[:,1],pointsall[:,2],s=0.5,c='#1f77b4')
    # slope, intercept, r_value, p_value, std_err = sp.stats.linregress(pointsall[:,1],pointsall[:,2])
    # print(slope, intercept, r_value, p_value, std_err)
    # xp = np.linspace(0.3, 1, 10)
    # yp = xp*slope + intercept
    # plt.plot(xp,yp,color='#000000')
    # plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    # plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    # # z = np.polyfit(pointsall[:,1],pointsall[:,2], 1)
    # # xp = np.linspace(0, 1, 10)
    # # yp = xp*z[0] + z[1]
    # # plt.plot(xp,yp)
    # plt.savefig(group_name[group]+'_all_correlation_workspace'+'.png')
    #
    # plt.figure(80)
    # plt.title('Index Frequency Content Correlation For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    # plt.scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
    # slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
    # print(slope, intercept, r_value, p_value, std_err)
    # xp = np.linspace(0, 1, 10)
    # yp = xp*slope + intercept
    # plt.plot(xp,yp,color='#000000')
    # plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    # plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    # plt.savefig(group_name[group]+'_timecorrelation_workspace'+'.png')

    # if group == 1:
    #     plt.figure(61)
    #     plt.title('High Frequency Content in Force Magnitude For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    #     plt.scatter(pointsall[:,1],pointsall[:,2],s=0.5,c='#1f77b4')
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(pointsall[:,1],pointsall[:,2])
    #     print(slope, intercept, r_value, p_value, std_err)
    #     xp = np.linspace(0.3, 1, 10)
    #     yp = xp*slope + intercept
    #     plt.plot(xp,yp,color='#000000')
    #     plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     # z = np.polyfit(pointsall[:,1],pointsall[:,2], 1)
    #     # xp = np.linspace(0, 1, 10)
    #     # yp = xp*z[0] + z[1]
    #     # plt.plot(xp,yp)
    #     plt.savefig(group_name[group]+'_correlation_workspace'+'.png')
    #
    #     plt.figure(80)
    #     plt.title('Index Frequency Content Correlation For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    #     plt.scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
    #     print(slope, intercept, r_value, p_value, std_err)
    #     xp = np.linspace(0, 1, 10)
    #     yp = xp*slope + intercept
    #     plt.plot(xp,yp,color='#000000')
    #     plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.savefig(group_name[group]+'_timecorrelation_workspace'+'.png')
    # else:
    #     plt.figure(50)
    #     plt.title('High Frequency Content in Force Magnitude For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    #     plt.scatter(pointsall[:,1],pointsall[:,2],s=0.5,c='#1f77b4')
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(pointsall[:,1],pointsall[:,2])
    #     print(slope, intercept, r_value, p_value, std_err)
    #     xp = np.linspace(0.3, 1, 10)
    #     yp = xp*slope + intercept
    #     plt.plot(xp,yp,color='#000000')
    #     plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.savefig(group_name[group]+'_correlation_workspace'+'.png')
    #
    #     plt.figure(70)
    #     plt.title('Index Frequency Content Correlation For '+group_name[group])#\n Sub '+str(subject_num)+' '+supportlevels[j])#supportlevels[j])#+' '+tasks[i-1])
    #     plt.scatter(points_timecorr[:,0],points_timecorr[:,1],s=0.5,c='#1f77b4')
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(points_timecorr[:,0],points_timecorr[:,1])
    #     print(slope, intercept, r_value, p_value, std_err)
    #     xp = np.linspace(0, 1, 10)
    #     yp = xp*slope + intercept
    #     plt.plot(xp,yp,color='#000000')
    #     plt.text(.5,1.8,'y='+str(round(slope,2))+'*x+'+str(round(intercept,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.text(.5,1.7,'r^2='+str(round(r_value**2,2))+'  p='+str(round(p_value,2)),horizontalalignment='left', fontsize=10, fontweight='bold')
    #     plt.savefig(group_name[group]+'_timecorrelation_workspace'+'.png')



plt.show()
