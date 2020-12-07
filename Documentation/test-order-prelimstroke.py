import csv
import numpy as np
import random

minsub = 2
maxsub = 2

for sub in range(minsub,maxsub+1):
    filename = "S0"+str(sub)+"_testorder_stroke.csv"
    columns = ['Frequency','trial']
    with open(filename, 'wb') as csvfile:
        testwriter = csv.writer(csvfile,delimiter = ',')
        testwriter.writerow(['ball-in-bowl'])
        testwriter.writerow(columns)

    # Define experiment
    num_support = 4
    num_repetitions = 3
    num_sets = 2
    num_arms = 2
    trialnum = 1

    for arm in range(num_arms):
        with open(filename,'ab') as csvfile:
            testwriter = csv.writer(csvfile,delimiter=',')
            if arm == 0:
                row=['Non-paretic']
            else:
                row=['Paretic']
            testwriter.writerow(row)

        # Arrange support order for each set
        SL_order = np.zeros((num_sets,num_support))
        order_repeat = True
        while order_repeat:
            for set in range(num_sets):
                SL_check = np.zeros(num_support)
                i = 0
                while i < num_support:
                    SL = random.randint(0,num_support-1)
                    if SL_check[SL]==0:
                        SL_check[SL] = 1
                        SL_order[set,i] = SL
                        i += 1

            order_repeat = False
            # check to see if there are >1 repeats in any column
            for i in range(4):
                count = np.zeros(num_support)
                for j in range(num_sets):
                    count[int(SL_order[j,i])] += 1
                if np.max(count) > 1:
                    order_repeat = True
                    # print('Try again to find a suitable SLuency order')
                    continue
            if SL_order[0,3]==SL_order[1,0]:
                order_repeat = True
                # print('Try again to find a suitable SLuency order')
                continue

        print(SL_order)

        with open(filename,'ab') as csvfile:
            testwriter = csv.writer(csvfile,delimiter=',')
            for set in range(num_sets):
                row=['Set '+str(set+1)]
                testwriter.writerow(row)
                for SL in range(num_support):
                    SL_level = int(SL_order[set,SL])
                    for trial_count in range(0,num_repetitions):
                        row=[SL_level,trialnum]
                        testwriter.writerow(row)
                        trialnum+=1


    columns = ['left-v-right','trial']
    with open(filename, 'ab') as csvfile:
        testwriter = csv.writer(csvfile,delimiter = ',')
        testwriter.writerow(['nail-and-hammer'])
        testwriter.writerow(columns)

    # Define experiment
    num_support = 4
    num_LR = 2
    num_repetitions = 3
    num_sets = 2
    num_arms = 2
    trialnum = 1

    for arm in range(num_arms):
        with open(filename,'ab') as csvfile:
            testwriter = csv.writer(csvfile,delimiter=',')
            if arm == 0:
                row=['Non-paretic']
            else:
                row=['Paretic']
            testwriter.writerow(row)

        # # Arrange support order for each set
        # SL_order = np.zeros((num_sets,num_support))
        # order_repeat = True
        # while order_repeat:
        #     for set in range(num_sets):
        #         SL_check = np.zeros(num_support)
        #         i = 0
        #         while i < num_support:
        #             SL = random.randint(0,num_support-1)
        #             if SL_check[SL]==0:
        #                 SL_check[SL] = 1
        #                 SL_order[set,i] = SL
        #                 i += 1
        #
        #     order_repeat = False
        #     # check to see if there are >1 repeats in any column
        #     for i in range(4):
        #         count = np.zeros(num_support)
        #         for j in range(num_sets):
        #             count[int(SL_order[j,i])] += 1
        #         if np.max(count) > 1:
        #             order_repeat = True
        #             # print('Try again to find a suitable SLuency order')
        #             continue
        #     if SL_order[0,3]==SL_order[1,0]:
        #         order_repeat = True
        #         # print('Try again to find a suitable SLuency order')
        #         continue
        # print(SL_order)

        LR_order = np.zeros((num_sets,num_support))
        order_repeat = True
        while order_repeat:
            # for set in range(num_sets):
            for i in range(num_support):
                LR = random.randint(0,num_LR-1)
                LR_order[0,i] = LR
                LR_order[1,i] = not LR

            if np.sum(LR_order[0,:]) == 2:
                order_repeat = False
        print(LR_order)

        with open(filename,'ab') as csvfile:
            testwriter = csv.writer(csvfile,delimiter=',')
            for set in range(num_sets):
                row=['Set '+str(set+1)]
                testwriter.writerow(row)
                # for i in range(num_support):
                #     SL_level = int(SL_order[set,i])
                LR_level = int(LR_order[set,i])
                for trial_count in range(0,num_repetitions):
                    row=[LR_level,trialnum]
                    testwriter.writerow(row)
                    trialnum+=1
