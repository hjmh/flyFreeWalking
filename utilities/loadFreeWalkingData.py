__author__ = 'hannah'

import numpy as np

def loadCtraxData(indat, keyList):
    # data columns to be extracted from ctrax file

    dat = [indat[k] for k in keyList]

    # Reorganise fly position arrays into numpy arrays
    numFrames = len(dat[0])
    time = np.zeros(numFrames)
    xPos = np.zeros(numFrames)
    yPos = np.zeros(numFrames)
    angle = np.zeros(numFrames)

    pointer = 0
    for t in range(numFrames):
        numFlies = dat[3][t].astype('int')[0]

        if not numFlies == 1:
            time[t] = dat[0][t]
            xPos[t] = np.nan
            yPos[t] = np.nan
            angle[t] = np.nan

        else: # everything is okay
            time[t] = dat[0][t]
            xPos[t] = dat[1][pointer]
            yPos[t] = dat[2][pointer]
            angle[t] = dat[5][pointer]

        pointer += numFlies

    return time, xPos, yPos, angle

