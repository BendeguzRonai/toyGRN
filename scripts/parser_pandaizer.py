#### IMPORTS ####
import numpy as np
import pandas as pd



#### INPUT PARSER FUNCTION #########################################################################################

def parser (path, verbose = False):
    
    times, meanFs, maxFs = [], [], []
    
    
    with open(path, "r") as file:
        lines = [line.split() for line in file]
        for i, l in enumerate(lines):
            try:
                if l[0] == "Experiment":
                    if lines[i+1][0] == "Time:" and lines[i+1][1] != "0":
                        times.append(int(lines[i+1][1]))
                        meanFs.append(float(lines[i+1][4]))
                        maxFs.append(float(lines[i+1][7]))
                    else:
                        try:
                            if lines[i+2][0] == "Time:":
                                times.append(int(lines[i+2][1]))
                                meanFs.append(float(lines[i+2][4]))
                                maxFs.append(float(lines[i+2][7]))
                            elif lines[i+1][0] == "Time:" and lines[i+1][1] == "0":
                                times.append(int(lines[i+1][1]))
                                meanFs.append(float(lines[i+1][4]))
                                maxFs.append(float(lines[i+1][7]))
                            else:
                                print ("Error: check file at " + path)
                        except:
                            try:
                                if lines[i+1][0] == "Time:" and lines[i+1][1] == "0":
                                    times.append(int(lines[i+1][1]))
                                    meanFs.append(float(lines[i+1][4]))
                                    maxFs.append(float(lines[i+1][7]))
                                else:
                                    print ("Error: check file at " + path)
                            except:    
                                print ("Error: check file at " + path)
                                continue
            except:    
                continue
    if verbose:
        print(len(times), len(meanFs), len(maxFs), path)

    return times, meanFs, maxFs



#### PANDAIZER FUNCTION ############################################################################################

def pandaizer (path, target, simtime = 20000, gridsize="50x50", mut="point", seed = 0, sel_str = 1, mixed = True, num_meas = 20):
    '''path, target, simtime = 20000, gridsize="50x50", mut="point", seed = 0, sel_str = 1, mixed = True, num_meas = 20'''
    times, meanFs, maxFs = parser (path)
    df = pd.DataFrame(
        {
            "target" : np.full(num_meas, target),
            "simtime" : np.full(num_meas, simtime),
            "gridsize" : np.full(num_meas, gridsize),
            "mutations" : np.full(num_meas, mut),
            "seed" : np.full(num_meas, seed),
            "sel_str" : np.full(num_meas, sel_str),
            "mixed" : np.full(num_meas, "Mixed" if mixed else "Unmixed"),
            "time" : times,
            "meanF" : meanFs,
            "maxF" : maxFs
        }
    )

    return df



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
            "maxF" : []
        }
    )

    return df


#### CONCATENATOR WRAPPER FUNCTION #################################################################################
def concatenator(*dfs):
    return pd.concat(dfs, axis=0, ignore_index=True)
