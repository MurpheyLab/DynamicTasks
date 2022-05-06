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

number_of_control_subjects = 6
number_of_stroke_subjects = 0
DT = 1./60.

list_of_participants = ["S202","S203","S207","S208","S209","S212"] #,"S211","S212"]
list_of_participants_num = ["202","203","207","208","209","S212"] #,"211","212"]
side_arm = ["L","R","R","L","R","L"]
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
        if (scorelist[i] % 7)==0 and flag==0:
            temp = np.delete(temp,tempindex, 0)
            tempindex = tempindex -1
            flag = 1
        if (scorelist[i] % 7)!=0:
            flag = 0
        tempindex = tempindex + 1

    TBmedianVAL = np.median(temp)
    # TBmodeVAL = np.mode(temp)
    TBaveVAL = np.mean(temp)

    return [TBmedianVAL,TBaveVAL] #TBmodeVAL,

def agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all):
    # Label factors as strings for ezANOVA analysis
    # supportlevels = ['0%Max','20%Max','50%Max']
    subname = list_of_participants[subject_num]

    scoretotal = 0; scorecount = 0
    # Loop through each file, compute metrics, and save to file/s
    for side in range(0,2): #iterate through paretic/nonparetic
        for suplev in range(0,2): #iterate through supportlevels
            for k in range(0,25): #iterate through trials
                #print(DIR+subname+subfile+"_Trial"+str(k)+"_Task2_SL"+str(suplev+1)+"_A"+str(side)+".csv")
                try:
                    trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Task2_SL"+str(suplev+1)+"_A"+str(side)+".csv"
                    #data = genfromtxt(trialfile,delimiter=',',dtype=float)
                    #data = np.delete(data,0,0) # Deletes the first column of column names
                    row = calc_metrics(trialfile)
                    print(row)
                    scoretotal = scoretotal + row[0]
                    print(DIR+subname+subfile+"_Trial"+str(k)+"_Task2_SL"+str(suplev+1)+"_A"+str(side)+".csv")
                    scorecount = scorecount + 1
                    #print(row)
                    # inserts the trial information before metrics
                    row.insert(0,suplev)
                    row.insert(0,side)
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
                testwriter_corr.writerow([list_of_participants[subject_num],'NaH',side, suplev,scoreave])


################################################################################
#         Generate an "all-metrics-comparison" file for comparison analysis    #
################################################################################

file_all = DIR+"all-metrics-nah.csv"
columns_all = ["Side","Dominance","Subject","TrialNum","Arm","Suplev","Score","TBmedian","TBave"]
with open(file_all,'w', newline="") as csvfile_all:
    testwriter_all = csv.writer(csvfile_all,delimiter=',')
    testwriter_all.writerow(columns_all)

##########################################################################
#         Generate an "subj-ave-scores" file for correlation analysis    #
##########################################################################

file_corr = DIR+"subj-ave-scores.csv"
columns_corr = ["Subject","Game","Arm","Suplev","Score"]
with open(file_corr,'a', newline="") as csvfile_corr:
    testwriter_corr = csv.writer(csvfile_corr,delimiter=',')
    testwriter_corr.writerow(columns_corr)

################################################################################
#                 Loop through all subjects                                    #
################################################################################
for subject_num in range(0,len(list_of_participants)):
    sub = list_of_participants[subject_num]
    subfile = "\\NailHammerData_" + list_of_participants[subject_num]
    print(list_of_participants[subject_num])
    agg_data(subject_num, sub, subfile, file_all, csvfile_all, testwriter_all)
