import csv
import numpy as np
import random

minsub = 1
maxsub = 2

for sub in range(minsub,maxsub+1):
    filename = "S0"+str(sub)+"_testorder.csv"
    columns = ['frequency','haptic_forces','trial']
    with open(filename, 'wb') as csvfile:
        testwriter = csv.writer(csvfile,delimiter = ',')
        testwriter.writerow(['ball-in-bowl'])
        testwriter.writerow(columns)

    # Define experiment
    num_hapticforcelevels = 3
    num_frequencies = 4
    num_repetitions = 3
    num_sets = 3

    # Arrange haptic force order
    HF_order = np.zeros((num_frequencies,num_hapticforcelevels*2))
    order_repeat = True
    while order_repeat:
        for freq in range(num_frequencies):
            HF_check = np.zeros(num_hapticforcelevels)
            i = 0
            while i < num_hapticforcelevels:
                HF = random.randint(0,num_hapticforcelevels-1)
                if HF_check[HF]==0:
                    HF_order[freq,i] = HF
                    HF_check[HF] += 1
                    i += 1
            while i < num_hapticforcelevels*2:
                HF_order[freq,i] = HF_order[freq,i-3]
                i += 1

        # check to see if there are >3 repeats in any column
        order_repeat = False
        for i in range(3):
            count = np.zeros(num_hapticforcelevels)
            for j in range(num_frequencies):
                count[int(HF_order[j,i])] += 1
            if np.max(count) > 2:
                order_repeat = True
                print('Try again to find a suitable haptic force order')
                continue

    print(HF_order)

    # Arrange frequency order for each set
    freq_order = np.zeros((num_sets,num_frequencies))
    order_repeat = True
    while order_repeat:
        for set in range(num_sets):
            freq_check = np.zeros(num_frequencies)
            i = 0
            while i < num_frequencies:
                freq = random.randint(0,num_frequencies-1)
                if freq_check[freq]==0:
                    freq_check[freq] = 1
                    freq_order[set,i] = freq
                    i += 1

        order_repeat = False
        # check to see if there are >1 repeats in any column
        for i in range(4):
            count = np.zeros(num_frequencies)
            for j in range(num_sets):
                count[int(freq_order[j,i])] += 1
            if np.max(count) > 1:
                order_repeat = True
                # print('Try again to find a suitable frequency order')
                continue
        if freq_order[0,3]==freq_order[1,0] or freq_order[1,3]==freq_order[2,0]:
            order_repeat = True
            # print('Try again to find a suitable frequency order')
            continue

    print(freq_order)

    trialnum = 1
    for set in range(num_sets):
        if set == 0:
            HF_index = 0
        elif set==1:
            HF_index = 2
        else:
            HF_index = 4
        for freq in range(num_frequencies):
            freq_level = int(freq_order[set,freq])

            with open(filename,'ab') as csvfile:
                testwriter = csv.writer(csvfile,delimiter=',')

                # First haptic force set
                HF_level = HF_order[freq_level,HF_index]
                for trial_count in range(0,num_repetitions):
                    row=[freq_level,HF_level,trialnum]
                    testwriter.writerow(row)
                    trialnum+=1

                # Second haptic force set
                HF_level = HF_order[freq_level,HF_index+1]
                for trial_count in range(0,num_repetitions):
                    row=[freq_level,HF_level,trialnum]
                    testwriter.writerow(row)
                    trialnum+=1

        # stationary ball
        if set < 2:
            with open(filename,'ab') as csvfile:
                testwriter = csv.writer(csvfile,delimiter=',')
                for trial_count in range(0,num_repetitions):
                    row=['stationary ball','',trialnum]
                    testwriter.writerow(row)
                    trialnum+=1
        with open(filename,'ab') as csvfile:
            testwriter = csv.writer(csvfile,delimiter=',')
            row=['Mandatory 5min break']
            testwriter.writerow(row)

    columns = ['left-v-right','trial']
    with open(filename, 'ab') as csvfile:
        testwriter = csv.writer(csvfile,delimiter = ',')
        testwriter.writerow(['nail-and-hammer'])
        testwriter.writerow(columns)

    # Define experiment
    num_LR = 2
    num_repetitions = 3
    num_sets = 2

    # For each set, randomly chose from all experimental combos
    trialnum = 1
    for i in range(0,num_sets):
        LR_check = np.zeros(num_LR)
        while np.sum(LR_check) < num_LR:
            LR = random.randint(0,num_LR-1)
            if LR_check[LR]==0:
                with open(filename,'ab') as csvfile:
                    testwriter = csv.writer(csvfile,delimiter=',')
                    for trial_count in range(0,num_repetitions):
                        row=[LR,trialnum]
                        testwriter.writerow(row)
                        trialnum+=1
                LR_check[LR]=1
