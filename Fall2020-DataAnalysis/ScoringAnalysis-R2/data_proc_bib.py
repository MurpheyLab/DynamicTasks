from numpy import genfromtxt
import numpy as np
import csv

################################################################################
# This program reads the NIH prelim data for all subjects, computes various
# metrics for each trial, and saves them to csv files: one for each subject
# within their file and one overall data sheet with all subjects. The csv files
# are formatted to be easily read by R for statistcal analysis,
################################################################################

#set directory where data is mounted
#DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/" # linux
DIR = "E:\\Research\\DowntownData\\Fall2020RoundTwo\\" # windows

number_of_control_subjects = 3
number_of_stroke_subjects = 0
DT = 0.05

list_of_participants = ["S01","S02","S03"] #,"S04","S05","S06","S07"]
list_of_participants_num = ["1","2","3"] #,"4","5","6","7"]
side_arm = ["R","L","R"]
dom_arm = ["D","N","D"]
#impairment = []
# for i in range(0,number_of_stroke_subjects):
#     impairment.append("Stroke")
# for i in range(number_of_stroke_subjects, number_of_stroke_subjects+number_of_control_subjects):
#     impairment.append("Control")

################################################################################
#               Define functions for generating metrics                        #
################################################################################

def calc_metrics(trialfile):
    """
    modified time to completion add 1 sec for every flag not collected
    """
    data = np.genfromtxt(trialfile,delimiter=',') #,dtype=float)
    data = np.delete(data,0,0) # Deletes the first column of column names
    score = data[-1,7] # Gets the last value of the score
    tinbowl = 0.
    for i in range(len(data[:,0])):
        if data[i,3]>-0.14:
            if data[i,6]<0.021:
                tinbowl+=DT
    #tinbowlpercent = tinbowltocompletion / 30.0
    #meanvel,maxvel,meanvellifted = calc_ballvel(data[:,7], tlifted_vec) #4:6 #if just 6, calculates the velocity in the z direction
    #meanjerk,maxjerk = calc_jerk(data[delayInd:,1:3], tlifted_vec) #1:3 #if just 1:2, calculates jerk in the xy plane

    meanenergy1 = np.mean(data[:,16])
    meanenergy2 = calc_ball_energy(data[:,4:7], len(data[:,0]))
    #print([score,tinbowl,meanenergy])
    return [score,tinbowl,meanenergy1,meanenergy2] # meanvel,maxvel,meanvellifted,meanjerk,maxjerk,intenergy,windowenergy]

def calc_ball_energy(pos, trial_len):
    vel = np.gradient(pos,axis=0)
    ballenergy_vec = np.zeros(trial_len)
    for i in range(trial_len):
        ballenergy_vec[i] = 9.81*pos[i,2]+0.5*pow(vel[i,0],2)+0.5*pow(vel[i,1],2)+0.5*pow(vel[i,2],2)
    # ballenergyDER_vec = np.gradient(ballenergy_vec,axis=0)
    # ballenergyDER_lifted = np.multiply(tlifted_vec, ballenergyDER_vec)
    # windowenergy = np.max(windowing(ballenergy_vec,100))
    # intenergy = np.sum(ballenergy_lifted)
    # meanenergy = np.nanmean(np.where(ballenergy_vec!=0,ballenergy_vec,np.nan))
    meanenergy = np.mean(ballenergy_vec)
    # medianenergy = np.nanmedian(np.where(ballenergy_lifted!=0,ballenergy_lifted,np.nan))
    return meanenergy #,intenergy,windowenergy]

def agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all):
    # Label factors as strings for ezANOVA analysis
    # supportlevels = ['0%Max','20%Max','50%Max']
    freq = ['Freq0','Freq1','Freq2','Freq3']

    # Loop through each file, compute metrics, and save to file/s
    for j in range(0,4): #frequency
        for k in range(1,80): #trials
            for i in range(0,2): #forces
                for b in range(0,2): #ball moving
                    if (i!=0 or b!=0):
                        try:
                            #print(DIR+sub+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL0_F"+str(i)+"_B"+str(b)+".csv")
                            row = calc_metrics(DIR+sub+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL0_F"+str(i)+"_B"+str(b)+".csv")
                            #print(row)
                            # inserts the trial information before metrics
                            row.insert(0,j) # frequency
                            row.insert(0,b) # ball moving
                            row.insert(0,i) # forces
                            if b==0 and i==1:
                                condition = 1
                            elif b==1 and i==1:
                                condition = 0
                            else:
                                condition = 2
                            row.insert(0,condition)
                            row.insert(0,k) # trial num
                            row.insert(0,list_of_participants[subject_num])
                            row.insert(0,dom_arm[subject_num])
                            row.insert(0,side_arm[subject_num])
                            with open (file_all,'a', newline="") as csvfile_all:
                                testwriter_all = csv.writer(csvfile_all,delimiter=',')
                                testwriter_all.writerow(row)
                        except:
                            pass


################################################################################
#         Generate an "all-metrics-comparison" file for comparison analysis    #
################################################################################

file_all = DIR+"all-metrics-bib.csv"
columns_all = ["Side","Dominance","Subject","TrialNum","Condition","Forces","BallMov","Freq","Score","TIB","MeanEnergy1","MeanEnergy2"]
with open(file_all,'w', newline="") as csvfile_all:
    testwriter_all = csv.writer(csvfile_all,delimiter=',')
    testwriter_all.writerow(columns_all)

################################################################################
#                 Loop through all subjects (stroke & control)                 #
################################################################################
for subject_num in range(0,len(list_of_participants)):
    sub = list_of_participants[subject_num]
    subfile = "\OutputData_S" + list_of_participants_num[subject_num]
    print(list_of_participants[subject_num])
    agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all)
