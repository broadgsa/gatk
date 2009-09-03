#!/util/bin/python

import math, string, sys

for input_file in sys.argv[1:]:
    
    FILENAME = input_file

    # Read the depth of coverage distribution from file
    depth = []; count = []
    for line in open(FILENAME):
        try:
            int(line[0])
            split = string.split(line,' ')
            if int(split[1]) != 0:
                depth.append(int(split[0]))
                count.append(int(split[1]))
        except ValueError: pass

    # Calculate the mode
    mode = depth[count.index(max(count))]

    # Find the index of maximum extent of 'good' data 
    idx = depth.index(mode); dist = 10**9
    while abs(count[idx] - count[0]) < dist:
        dist = abs(count[idx] - count[0])
        idx += 1
    model_count = count[:idx + 1]

    # Calculate mean & variance of 'good' data
    model_tot = sum(model_count); adder = 0
    for i in range(len(model_count)): adder += depth[i]*model_count[i]
    mean = adder/float(model_tot)

    adder = 0
    for i in range(len(model_count)): adder += model_count[i]*(mean - depth[i])**2
    var = adder/float(model_tot)
    stdev = var**0.5

    # Output Thresholds
    print FILENAME
    print 'Mean Depth of Coverage: ' + str(mean)
    print 'Plus One Std Dev: ' + str(mean + stdev)
    print 'Plus Two Std Dev: ' + str(mean + 2*stdev)
    print 'Plus Three Std Dev: ' + str(mean + 3*stdev)
    print 'Plus Four Std Dev: ' + str(mean + 4*stdev)
    print 'Plus Five Std Dev: ' + str(mean + 5*stdev)
