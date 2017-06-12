'''

Function to import Ctrax tracking data from mat files.

'''

__author__ = 'hannah'


def importCtraxData(npzFilePath, matFilePath):
    keyList = ['timestamps', 'x_pos', 'y_pos', 'ntargets', 'identity',
               'angle']  # data columns to be extracted from ctrax file

    try:
        # load existing python file
        npzData = np.load(npzFilePath)
        xPos = npzData['xPos']
        yPos = npzData['yPos']
        angle = npzData['angle']
        flyID = npzData['flyID']
        numFrames = len(xPos)
        maxFlyID = max(flyID)

    except:
        # load matlab data and convert
        indat = loadmat(matFilePath)

        dat = [indat[k] for k in keyList]

        # Reorganise fly position arrays into lists (sorted by frame)
        numFrames = len(dat[0])
        xPos = []
        yPos = []
        angle = []
        flyID = []

        pointer = 0
        for t in range(numFrames):
            numFlies = dat[3][t].astype('int')[0]

            xPos.append(dat[1][pointer:pointer + numFlies])
            yPos.append(dat[2][pointer:pointer + numFlies])
            angle.append(dat[5][pointer:pointer + numFlies])
            flyID.append(dat[4][pointer:pointer + numFlies])

            pointer += numFlies

        xPos = np.array(xPos)
        yPos = np.array(yPos)
        angle = np.array(angle)
        flyID = np.array(flyID)
        maxFlyID = max(dat[4])

        # save data for future sessions
        np.savez(baseDir + '/' + genotype + '/' + genotype + '_' + experiment + '/' + fileName +
                 '.npz', xPos=xPos, yPos=yPos, angle=angle, flyID=flyID)

    return xPos, yPos, angle, flyID, numFrames, maxFlyID
