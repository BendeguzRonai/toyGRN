#### IMPORTS #######################################################################################################
import numpy as np
import pandas as pd

#### QUICK EMPTY DATAFRAME #########################################################################################


def empty ():
    
    df = pd.DataFrame(
        {
            "target" : [],
            "simtime" : [],
            "gridsize" : [],
            "mutations" : [],
            "seed" : [],
            "sel_str" : [],
            "mixed" : [],
            "time" : [],
            "meanF" : [],
            "maxF" : [],
            "mean_distance" : [],
            "least_distance" : [],
            "size" : []
        }
    )

    return df

#### INPUT PARSER FUNCTION #########################################################################################

def parser (path, verbose = False):

    ret_times = []
    ret_meanFs = []
    ret_maxFs = []
    ret_meanDs = []
    ret_leastDs = []
    ret_sizes = []
    current_exp = None

    try:
        with open (path, "r") as f:
            
            lines = [line.split() for line in f]
            for l in lines:
                try:
                    
                    if "Experiment" in l[0]:
                        if current_exp != None:
                            ret_times.append(np.array(times))
                            ret_meanFs.append(np.array(meanFs))
                            ret_maxFs.append(np.array(maxFs))
                            ret_meanDs.append(np.array(meanDs))
                            ret_leastDs.append(np.array(leastDs))
                            ret_sizes.append(np.array(sizes))
                        current_exp = int(l[1])
                        times = []
                        meanFs = []
                        maxFs = []
                        meanDs = []
                        leastDs = []
                        sizes = []
                        
                    if current_exp != None and "Time" in l[0]:
                        times.append(int(l[1]))
                        meanFs.append(float(l[4]))
                        maxFs.append(float(l[7]))
                        meanDs.append(float(l[10]))
                        leastDs.append(float(l[13]))
                        sizes.append(int(l[15]))
                except:
                    continue
                    
            if current_exp != None:
                        ret_times.append(np.array(times))
                        ret_meanFs.append(np.array(meanFs))
                        ret_maxFs.append(np.array(maxFs))
                        ret_meanDs.append(np.array(meanDs))
                        ret_leastDs.append(np.array(leastDs))
                        ret_sizes.append(np.array(sizes))
        if verbose:
            print(len(ret_times), len(ret_meanFs), len(ret_maxFs), len(ret_meanDs), len(ret_leastDs), len(ret_sizes))
            
    except Exception as e:
        print ("Error occured, check file at " + path)
        print(e)

    return ret_times, ret_meanFs, ret_maxFs, ret_meanDs, ret_leastDs, ret_sizes


#### PADDER FUNCTION ###############################################################################################

def found_padder (array, simtime, intervals, istime = False):

    desired_length = int(simtime / intervals + 1)

    if not istime:
        return np.pad(array, (0, desired_length - len(array)), mode="constant", constant_values = array[-1])
    
    elif istime:
        lasttime = array [-1]
        if lasttime > simtime - 10000:
            return array    
        elif lasttime % 10000 == 0:
            next10k = lasttime + 10000
        else:
            next10k = lasttime - (lasttime % 10000) + 20000
        return np.append(array, [i * 10000 + next10k for i in range(desired_length - len(array))])

        
        



#### PANDAIZER FUNCTION ############################################################################################
    
def pandaizer (path, target,
               simtime = 20000,
               gridsize="50x50",
               mut="point",
               seed = 0,
               sel_str = 1,
               mixed = True,
               num_meas = 20,
               verbose=False):
    '''path, target, simtime = 20000, gridsize="50x50", mut="point", seed = 0, sel_str = 1, mixed = True, num_meas = 20'''

    times, meanFs, maxFs, meanDs, leastDs, sizes = parser(path, verbose=verbose)

    if verbose:
        print(path, target, simtime, mut, "ss", sel_str, "mixed", mixed)
    df = pd.DataFrame(
        {
            "target" : np.full(num_meas, target),
            "simtime" : np.full(num_meas, simtime),
            "gridsize" : np.full(num_meas, gridsize),
            "mutations" : np.full(num_meas, mut),
            "seed" : np.full(num_meas, seed),
            "sel_str" : np.full(num_meas, sel_str),
            "mixed" : np.full(num_meas, "Mixed" if mixed else "Unmixed"),
            "time" : [found_padder(i, simtime, 10000, istime=True) for i in times],
            "meanF" : [found_padder(i, simtime, 10000) for i in meanFs],
            "maxF" : [found_padder(i, simtime, 10000) for i in maxFs],
            "mean_distance" : [found_padder(i, simtime, 10000) for i in meanDs],
            "least_distance" : [found_padder(i, simtime, 10000) for i in leastDs],
            "size" : [found_padder(i, simtime, 10000) for i in sizes]
        }
    )

    return df


#### CONCATENATOR WRAPPER FUNCTION #################################################################################

def concatenator(*dfs):
    return pd.concat(dfs, axis=0, ignore_index=True)
