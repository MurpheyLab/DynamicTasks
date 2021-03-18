import csv
import numpy as np
import random

minsub = 1
maxsub = 8

for sub in range(minsub,maxsub+1):
    filename = "S0"+str(sub)+"_testorder_stroke.csv"
    with open(filename, 'wb') as csvfile:
        testwriter = csv.writer(csvfile,delimiter = ',')
        testwriter.writerow(['DAY 1---WORKSPACE SETUP'])

        # Define experiment
        num_support = 2
        num_freq = 4
        num_conditions = num_support*num_freq #number of experimental conditions with
        num_repetitions = 5
        num_arms = 2
        arms_list = ['Nonparetic','Paretic']
        trialnum = 1

        # create list of experimental trial options
        SL_options = np.zeros(num_conditions)
        F_options = np.zeros(num_conditions)
        for i in range(num_conditions):
            SL_options[i] = i // 4
            F_options[i] = i % num_freq

        arm_first = random.randint(0,num_arms-1)
        for day in range(2):
            for arm in range(num_arms):
                if arm==0:
                    if day==1:
                        arm_first = not(arm_first)
                    arm_i = arm_first
                    if day==0:
                        testwriter.writerow([arms_list[not(arm_i)],'30%max (SL1)'])
                        testwriter.writerow([arms_list[not(arm_i)],'0%max (SL0)'])
                    else:
                        testwriter.writerow([])
                        testwriter.writerow([])
                        testwriter.writerow(['DAY 2'])
                    testwriter.writerow([arms_list[not(arm_i)],'practice ball-in-bowl'])
                    if day==0:
                        testwriter.writerow([arms_list[arm_i],'30%max (SL1)'])
                        testwriter.writerow([arms_list[arm_i],'0%max (SL0)'])
                    testwriter.writerow([arms_list[arm_i],'practice ball-in-bowl'])
                    testwriter.writerow(['DAY '+str(day+1)+'---BALL-IN-BOWL'])
                    columns = ['Arm','Supportlevel','Frequency','trial']
                    testwriter.writerow(columns)
                elif arm==1:
                    arm_i = not(arm_first)
                    testwriter.writerow(['Mandatory 5min break + switch arms'])

                # set experimental order
                trial_order = np.arange(0,num_conditions,1)
                random.shuffle(trial_order)

                for trial in range(num_conditions):
                    if trial == 4:
                        testwriter.writerow(['Mandatory 5min break'])
                    trial_int = trial_order[trial]
                    for trial_count in range(0,num_repetitions):
                        row=[arms_list[arm_i],SL_options[trial_int],F_options[trial_int],trialnum]
                        testwriter.writerow(row)
                        trialnum+=1

        # NAIL AND HAMMER
        trialnum = 1
        for arm in range(num_arms):

            if arm==0:
                arm_i = not(arm_i)
                testwriter.writerow(['Mandatory 5min break'])
                testwriter.writerow([])
                testwriter.writerow(['NAIL-AND-HAMMER'])
                testwriter.writerow(['Arm','SupportLevel'])
                testwriter.writerow([arms_list[not(arm_i)],'training'])
                testwriter.writerow([arms_list[arm_i],'training'])
            elif arm==1:
                arm_i = not(arm_i)
                testwriter.writerow(['Mandatory 5min break + switch arms'])

            SL_int = random.randint(0,num_support-1)
            for trial_count in range(0,num_repetitions):
                row=[arms_list[arm_i],SL_int,trialnum]
                testwriter.writerow(row)
                trialnum+=1

            for trial_count in range(0,num_repetitions):
                row=[arms_list[arm_i],int(not(SL_int)),trialnum]
                testwriter.writerow(row)
                trialnum+=1