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
DIR = "E:\\Research\\DowntownWork\\DowntownData\\Stroke2021\\" # windows

number_of_control_subjects = 0
number_of_stroke_subjects = 3
DT = 1./60. #0.05

list_of_participants = ["S202","S203","S208","S209","S211","S212"] #["S201","S202","S203","S204","S205","S206","S207"]
list_of_participants_num = ["202","203","208","209","211","212"] #["1","2","3","4","5","6","7"]
side_arm = ["L","R","L","R","L","R"]
dom_arm = ["D","D","D","D","D","D"]
#impairment = []
# for i in range(0,number_of_stroke_subjects):
#     impairment.append("Stroke")
# for i in range(number_of_stroke_subjects, number_of_stroke_subjects+number_of_control_subjects):
#     impairment.append("Control")

################################################################################
#               Define functions for generating metrics                        #
################################################################################

def calc_metrics(trialfile):
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
    scoretotal = 0; scorecount = 0
    scorematrix = np.zeros((2,2,4,160))
    # Loop through each file, compute metrics, and save to file/s
    for i in range(0,2): #arm
        for supportlevel in range(1,3):
            for j in range(0,4): #frequency
                scoretotal = 0; scorecount = 0
                for k in range(1,160): #trials
                    try:
                        row = calc_metrics(DIR+sub+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL"+str(supportlevel)+"_A"+str(i)+".csv")
                        print(DIR+sub+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL"+str(supportlevel)+"_A"+str(i)+".csv")
                        scorematrix[i,supportlevel-1,j,k] = row[0]
                        scoretotal = scoretotal + row[0]
                        scorecount = scorecount + 1
                        # inserts the trial information before metrics
                        row.insert(0,j) # frequency
                        row.insert(0,i) # arm
                        row.insert(0,supportlevel)
                        row.insert(0,k) # trial num
                        row.insert(0,list_of_participants[subject_num])
                        row.insert(0,dom_arm[subject_num])
                        row.insert(0,side_arm[subject_num])
                        with open (file_all,'a', newline="") as csvfile_all:
                            testwriter_all = csv.writer(csvfile_all,delimiter=',')
                            testwriter_all.writerow(row)
                    except:
                        pass
                print(scoretotal, scorecount)
                scoreave = scoretotal / scorecount

                # find average of best 3 trials
                # bestscores = [0,0,0]
                # for nn in range (0,3):
                #      bestscores[nn] = np.amax(scorematrix[i,supportlevel-1,j,:])
                #      index = np.argmax(scorematrix[i,supportlevel-1,j,:])
                #      scorematrix[i,supportlevel-1,j,index] = 0
                #      scoreavebest = np.average(bestscores)

                scoreavebest = np.amax(scorematrix[i,supportlevel-1,j,:])

                with open (file_corr,'a', newline="") as csvfile_corr:
                    testwriter_corr = csv.writer(csvfile_corr,delimiter=',')
                    testwriter_corr.writerow([list_of_participants[subject_num],'BiB'+str(j)+'Hz', 'Arm'+str(i), 'SL'+str(supportlevel), scoreavebest])


################################################################################
#         Generate an "all-metrics-comparison" file for comparison analysis    #
################################################################################

file_all = DIR+"all-metrics-bib-stroke.csv"
columns_all = ["Side","Dominance","Subject","TrialNum","SupportLevel","Arm","Freq","Score","TIB","MeanEnergy1","MeanEnergy2"]
with open(file_all,'w', newline="") as csvfile_all:
    testwriter_all = csv.writer(csvfile_all,delimiter=',')
    testwriter_all.writerow(columns_all)

##########################################################################
#         Generate an "subj-ave-scores" file for correlation analysis    #
##########################################################################

file_corr = DIR+"subj-ave-scores.csv"

################################################################################
#                 Loop through all subjects (stroke & control)                 #
################################################################################
for subject_num in range(0,len(list_of_participants)):
    sub = list_of_participants[subject_num]
    subfile = "\OutputData_S" + list_of_participants_num[subject_num]
    print(list_of_participants[subject_num])
    agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all)
