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

number_of_control_subjects = 6
number_of_stroke_subjects = 0
DT = 1./60.

list_of_participants = ["S01","S02","S03","S04","S05","S06","S07"]
list_of_participants_num = ["1","2","3","4","5","6","7"]
side_arm = ["R","L","R","L","L","R","R"]
dom_arm = ["D","N","D","N","N","D","D"]
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
    score = data[-1,12] # Gets the last value of the score

    [TBmedian,TBave] = calculateTBmetrics(data)

    return [score,TBmedian,TBave] #TBmode

def calculateTBmetrics(data):
    temp = []
    scorelist = []
    temptime = 0
    currentscore = 0
    flag = 0
    for i in range(1,len(data[:,1])):
        #### calculate TB based on subsections
        # if int(data[i,12])!=int(currentscore):
        #     temp[int(currentscore)] = data[i,0] - temptime
        #     currentscore = currentscore+1
        #     temptime = data[i,0]
        #### calculate TB based on movement
        if int(data[i,7])>0 and int(data[i-1,7])>0 and flag==0:
            temp.append(data[i,0] - temptime)
            scorelist.append(data[i,12])
            temptime = data[i,0]
            flag = 1
        elif int(data[i,7])==0 and int(data[i-1,7])==0 and flag==1:
            flag = 0

    #print("temp: ", temp, "\n")
    #print("score list: ", scorelist, "\n")

    tempindex = 0
    flag = 0
    for i in range(len(temp)):
        if (scorelist[i] % 5)==0 and flag==0:
            temp = np.delete(temp,tempindex, 0)
            tempindex = tempindex -1
            flag = 1
        if (scorelist[i] % 5)!=0:
            flag = 0
        tempindex = tempindex + 1

    TBmedianVAL = np.median(temp)
    # TBmodeVAL = np.mode(temp)
    TBaveVAL = np.mean(temp)

    return [TBmedianVAL,TBaveVAL] #TBmodeVAL,

def agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all):
    # Label factors as strings for ezANOVA analysis
    # supportlevels = ['0%Max','20%Max','50%Max']

    scoretotal = 0; scorecount = 0
    # Loop through each file, compute metrics, and save to file/s
    for j in range(0,2): #R/L
        for k in range(1,14): #trials
            category = 2
            try:
                #print(DIR+sub+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL0"+"_F"+str(i)+".csv")
                row = calc_metrics(DIR+sub+subfile+"_Trial"+str(k)+"_Task"+str(j)+"_SL0.csv")
                scoretotal = scoretotal + row[0]
                scorecount = scorecount + 1
                #print(row)
                # inserts the trial information before metrics
                if ((side_arm[subject_num]=="R" and j==0) or (side_arm[subject_num]=="L" and j==1)):
                    category = 0
                else:
                    category = 1
                row.insert(0,category)
                row.insert(0,k) # trial num
                row.insert(0,list_of_participants[subject_num])
                row.insert(0,dom_arm[subject_num])
                row.insert(0,side_arm[subject_num])
                with open (file_all,'a', newline="") as csvfile_all:
                    testwriter_all = csv.writer(csvfile_all,delimiter=',')
                    testwriter_all.writerow(row)
            except:
                pass

    scoreave = scoretotal / scorecount
    with open (file_corr,'a', newline="") as csvfile_corr:
        testwriter_corr = csv.writer(csvfile_corr,delimiter=',')
        testwriter_corr.writerow([list_of_participants[subject_num],'NaH',scoreave])


################################################################################
#         Generate an "all-metrics-comparison" file for comparison analysis    #
################################################################################

file_all = DIR+"all-metrics-nah.csv"
columns_all = ["Side","Dominance","Subject","TrialNum","category","Score","TBmedian","TBave"]
with open(file_all,'w', newline="") as csvfile_all:
    testwriter_all = csv.writer(csvfile_all,delimiter=',')
    testwriter_all.writerow(columns_all)

##########################################################################
#         Generate an "subj-ave-scores" file for correlation analysis    #
##########################################################################

file_corr = DIR+"subj-ave-scores.csv"
columns_corr = ["Subject","Game","Score"]
with open(file_corr,'a', newline="") as csvfile_corr:
    testwriter_corr = csv.writer(csvfile_corr,delimiter=',')
    testwriter_corr.writerow(columns_corr)

################################################################################
#                 Loop through all subjects                                    #
################################################################################
for subject_num in range(0,len(list_of_participants)):
    sub = list_of_participants[subject_num]
    subfile = "\\NailHammerData_S" + list_of_participants_num[subject_num]
    print(list_of_participants[subject_num])
    agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all)
