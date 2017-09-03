""" Funtion for analysing a single FlyOver trial (assumes use of FlyOver version 9.2 or 9.3) """

import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import matplotlib.colors as colors
import seaborn as sns

from os import listdir, mkdir
from glob import glob
from os.path import isfile, sep
from sys import path, exit

# Set path to analysis code directory
codeDir = sep.join(os.getcwd().split(sep)[:-2])
path.insert(1, codeDir)

# Import custom plotting functions
from plottingUtilities.basicPlotting import *
from plottingUtilities.velocityDistributionPlots import *

# Import custom data processing functions
from flyVR.utilities.loadSingleTrial import loadSingleVRLogfile
from flyVR.utilities.loadObjectCoords import loadObjectCoords

from trajectoryAnalysis.downsample import *
from trajectoryAnalysis.trajectoryDerivedParams import *
from trajectoryAnalysis.periodicWorldAnalysis import *
from trajectoryAnalysis.downsample import *

sns.set_style('ticks')

def singleTrialAnalysis(fileToAnalyse):

    # ------------------------------------------------------------------------------------------------------------------
    # Extract folder and file name info
    # ------------------------------------------------------------------------------------------------------------------

    print('Data file: \n' + fileToAnalyse + '\n')

    dataDir = sep.join(fileToAnalyse.split(sep)[0:-3]) + sep
    flyID = fileToAnalyse.split(sep)[-2]
    expDir = sep.join(fileToAnalyse.split(sep)[0:-1]) + sep

    if not ('males' in expDir.split(sep)[-4] or 'females'in expDir.split(sep)[-4]):
        print('You selected an invalid data directory.\n' +
              'Expected folder structure of the selected path is some/path/to/experimentName/flyGender/rawData/')
        exit(1)

    FODataFile = fileToAnalyse.split(sep)[-1]
    FODataFiles = [filepath.split(sep)[-1] for filepath in glob(expDir + '*.txt')]
    FODataFiles = sorted(FODataFiles)

    trial = FODataFiles.index(FODataFile) + 1

    # create analysis dir
    analysisDir = dataDir + 'analysis/'
    try:
        mkdir(analysisDir)
    except OSError:
        print('Analysis directory already exists.')

    dataFileParts = FODataFile.split('_')
    titleString = 'fly ' + flyID + ' in ' + dataFileParts[2] + ' of ' + dataFileParts[3] + ', trial' + str(trial)
    print('Analysing experiment...\n')
    print(titleString + '\n')

    if 'Invisible' in FODataFile or 'invisible' in FODataFile:
        objecttype = 'invisible'
    elif 'Empty' in FODataFile or 'empty' in FODataFile:
        objecttype = 'none'
    else:
        objecttype = 'visible'

    # ------------------------------------------------------------------------------------------------------------------
    # Load or read in logfile data
    # ------------------------------------------------------------------------------------------------------------------

    # Read in logfile data, parse header ...............................................................................
    # TODO

    # Read in object coordinates .......................................................................................
    # TODO

    # Compute movement velocities ......................................................................................
    logTime = np.copy(FOData[:, 0])
    time = np.linspace(0, logTime[-1], numFrames)

    angle = convertRawHeadingAngle(FOData[:, 5])

    samplingRate = 20
    N = 5
    vTrans, vRot, vTransFilt, vRotFilt = velocityFromTrajectory(time, angle, FOData[:, 1], FOData[:, 2], N, numFrames)

    # ------------------------------------------------------------------------------------------------------------------
    # Generate basic analysis plots
    # ------------------------------------------------------------------------------------------------------------------
    # Time step plot ...................................................................................................

    tstpfig = plt.figure(figsize=(10, 3))
    gs = gridspec.GridSpec(1, 2, width_ratios=np.hstack((2, 1)))
    tstpfig.suptitle(titleString, fontsize=14)
    histRange = (0, 0.011)

    ax = tstpfig.add_subplot(gs[0])
    ax.plot(FOData[0:-2, 0], (FOData[1:-1, 0]-FOData[0:-2, 0]).astype('float'), '.', alpha=0.1)
    ax.set_ylim(histRange)
    ax.set_xlim((0, time[-1]))
    ax.set_xlabel('time [s]')
    ax.set_ylabel('time step [1/s]')
    myAxisTheme(ax)

    ax = tstpfig.add_subplot(gs[1])
    ax.hist(FOData[1:-1, 0]-FOData[0:-2, 0], 50, histRange)
    ax.set_xlabel('time step [1/s]')
    ax.set_ylabel('count')
    ax.set_title('mean time step = ' + str(round(np.mean((FOData[1:-1, 0]-FOData[0:-2, 0])*1000.0), 2)) + 'ms')
    myAxisTheme(ax)

    tstpfig.tight_layout()

    makeNestedPlotDirectory(analysisDir, 'timeStepPlot/', objecttype + sep)

    tstpfig.savefig(analysisDir + 'timeStepPlot/' + objecttype + sep
                    + FODataFile[0:-4] + '_timeStepPlot_trial' + str(trial) + '.pdf', format='pdf')

    # Plot of walking trace (+ color bar for time) with object locations ................................................
    tStart = 0
    tEnd = len(FOData[:, 1])
    tStep = 72
    frameRange = range(tStart, tEnd, tStep)
    colMap = 'Accent'
    arrowLength = 5

    trajfig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=np.hstack((10, 1)))

    axTraj = trajfig.add_subplot(gs[0])
    axTime = trajfig.add_subplot(gs[1])

    plotPosInRange(axTraj, axTime, frameRange, FOData[:, 0], FOData[:, 1], FOData[:, 2], np.pi/180*FOData[:, 5],
                   colMap, arrowLength, 0.5, 5)
    axTraj.scatter(visibleObjectCoords[:, 0], visibleObjectCoords[:, 1], 50, alpha=0.5,
                   facecolor='black', edgecolors='none')
    axTraj.scatter(invisibleObjectCoords[4:, 0], invisibleObjectCoords[4:, 1], 50, alpha=0.5,
                   facecolors='none', edgecolors='black')
    axTraj.set_xlabel(header[1], fontsize=12)
    axTraj.set_ylabel(header[2], fontsize=12)
    axTraj.set_title('Walking trace of ' + titleString, fontsize=14)
    axTraj.set_xlim([max(-400, min(FOData[:, 1]) - 20), min(400, max(FOData[:, 1]) + 20)])
    axTraj.set_ylim([max(-400, min(FOData[:, 2]) - 20), min(400, max(FOData[:, 2]) + 20)])
    myAxisTheme(axTraj)

    axTime.set_xlabel(header[0], fontsize=12)
    plt.xlim((0, FOData[-1, 0]))
    timeAxisTheme(axTime)

    makeNestedPlotDirectory(analysisDir, 'tracePlot/', objecttype + sep)

    trajfig.savefig(analysisDir + 'tracePlot/' + objecttype + sep
                    + FODataFile[0:-4] + '_traceObjectPlot_trial' + str(trial) + '.pdf', format='pdf')

    # Plot velocity distributions of downs sampled data ................................................................
    rotLim = (-5, 5)
    transLim = (0, 30)
    angleLim = (-np.pi, np.pi)
    summaryVeloFig_ds = velocitySummaryPlot(time_ds, vTrans_ds, vTransFilt_ds, vRot_ds, vRotFilt_ds, angle_ds, rotLim,
                                            transLim, angleLim, 'Down sampled velocity traces, ' + titleString)

    makeNestedPlotDirectory(analysisDir, 'velocityTraces/', objecttype + sep)

    summaryVeloFig_ds.savefig(analysisDir + 'velocityTraces/' + objecttype + sep
                              + FODataFile[0:-4] + '_veloTraces_ds_trial' + str(trial) + '.pdf', format='pdf')

    # ------------------------------------------------------------------------------------------------------------------
    # Compute heading angle relative to closest object (use collapsed coordinates)
    # ------------------------------------------------------------------------------------------------------------------

    # Compute parameters characterising fly's relationship to object
    objLocation = [0, 0]
    objDirection, objDistance, gammaFull, gamma, gammaV\
        = relationToObject(time_ds, xPosMA_ds, yPosMA_ds, angle_ds, objLocation)

    # Change in heading rel. to object .................................................................................
    sf = 0
    ef = len(xPosMA_ds)
    near = 6
    far = arenaRad
    vTransTH = 2

    turnTH = 3*np.std(abs(vRotFilt_ds))
    turnMask = (abs(vRotFilt_ds) > turnTH)

    selectedRangeDistAll = np.logical_and(objDistance > near, objDistance < far)

    selectedRangeDist = np.logical_and(np.logical_and(objDistance > near, objDistance < far), vTrans_ds > vTransTH)
    selectedRangeDistTurn = np.logical_and(np.logical_and(objDistance > near, objDistance < far), turnMask)

    headTurnFig = plt.figure(figsize=(15, 5))

    headTurnFig.suptitle('Relationship between turns and relativel heading angle\n' + flyID + ', trial' + str(trial) +
                         ', percentage turns: ' + str(round(100.0*sum(turnMask)/len(vRotFilt_ds), 2)) +
                         ' (within range ' + str(near) + '-' + str(far) + 'mm, trans. vel. > ' + str(vTransTH) + ')\n',
                         fontsize=13)

    distRange = (near, far)
    angleRange = (0, np.pi)
    vHeadRange = (-5, 5)

    ax0 = headTurnFig.add_subplot(131)
    ax0.plot(xPosMA_ds[sf:ef], yPosMA_ds[sf:ef], '.', color='grey', alpha=0.4)
    ax0.plot(xPosMA_ds[turnMask[sf:ef]], yPosMA_ds[turnMask[sf:ef]], '.', color='red', alpha=0.7)
    ax0.plot(0, 0, marker='o', markersize=15, linestyle='none', alpha=0.5, color='black')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    ax0.set_aspect('equal')
    myAxisTheme(ax0)

    ax1 = headTurnFig.add_subplot(132)
    niceScatterPlot(ax1, objDistance[selectedRangeDist], gamma[selectedRangeDist], distRange, angleRange, 'grey', 0.3)
    niceScatterPlot(ax1, objDistance[selectedRangeDistTurn], gamma[selectedRangeDistTurn], distRange, angleRange,
                    'red', 0.7)
    plt.xlabel('distance from object')
    plt.ylabel('heading angle')

    ax2 = headTurnFig.add_subplot(133)
    niceScatterPlot(ax2, gammaV[selectedRangeDist], gamma[selectedRangeDist], vHeadRange, angleRange, 'grey', 0.3)
    niceScatterPlot(ax2, gammaV[selectedRangeDistTurn], gamma[selectedRangeDistTurn], vHeadRange, angleRange,
                    'red', 0.7)
    plt.xlabel('change in heading angle (pos - increase, neg - decrease)')
    plt.ylabel('heading angle')

    headTurnFig.tight_layout()

    makeNestedPlotDirectory(analysisDir, 'relativeHeading/', objecttype + sep)

    headTurnFig.savefig(analysisDir + 'relativeHeading/' + objecttype + sep
                        + FODataFile[0:-4] + '_headingAndTurns_ds_trial' + str(trial) + '.pdf', format='pdf')

    # Visualise effect of turns ........................................................................................
    turnEffectFig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(3, 3, height_ratios=np.hstack((1, 1.5, 1.5)), width_ratios=np.hstack((1.5, 1, 1)))

    ax0 = turnEffectFig.add_subplot(gs[:, 0])
    ax0.plot(xPosMA_ds[selectedRangeDist], yPosMA_ds[selectedRangeDist], '.', color='grey', alpha=0.4)
    ax0.plot(xPosMA_ds[selectedRangeDistTurn], yPosMA_ds[selectedRangeDistTurn], '.', color='red', alpha=0.7)
    ax0.plot(0, 0, marker='o', markersize=15, linestyle='none', alpha=0.5, color='black')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    ax0.set_xlim((-arenaRad, arenaRad))
    ax0.set_ylim((-arenaRad, arenaRad))
    ax0.set_aspect('equal')
    myAxisTheme(ax0)
    ax0.set_title('Effect of turns in ' + flyID + ', trial' + str(trial) +
                  '\n percentage turns: ' + str(round(100.0*sum(turnMask)/len(vRotFilt_ds), 2)) + '\n', fontsize=13)

    ax = turnEffectFig.add_subplot(gs[0, 1:3])
    gammaFullSelect = gammaFull[selectedRangeDistAll]
    gammaSelect = gamma[selectedRangeDistAll]
    plt.hist(gammaFullSelect[~np.isnan(gammaFullSelect)], bins=50, color='lightblue', alpha=0.8)
    plt.hist(gammaSelect[~np.isnan(gammaSelect)], bins=50, color='grey', alpha=0.5)
    plt.xlabel('relative heading angle [rad], not filtered for vTrans > ' + str(vTransTH))
    plt.ylabel('count')
    ax.set_xlim((-np.pi, np.pi))
    myAxisTheme(ax)

    headingHist = plotVeloHistogram_fancy(gamma[selectedRangeDist], gs[1, 1], (0, np.pi), 'grey', 0.5)
    headingHist.set_ylabel('walking filtered\n count (vTrans > ' + str(vTransTH) + ')')

    headingHistTurn = plotVeloHistogram_fancy(gamma[selectedRangeDistTurn], gs[2, 1], (0, np.pi), 'red', 0.5)
    headingHistTurn.set_ylabel('walking & turn filtered\n count')
    headingHistTurn.set_xlabel('heading angle\n [rad]')

    rotVHist = plotVeloHistogram_fancy(gammaV[selectedRangeDist], gs[1, 2], (-10, 10), 'grey', 0.5)
    rotVHist.set_xlabel('change in heading while walking\n [rad/s]')

    rotVHistFilt = plotVeloHistogram_fancy(gammaV[selectedRangeDistTurn], gs[2, 2], (-10, 10), 'red', 0.5)
    rotVHistFilt.set_xlabel('change in heading during turn\n [rad/s]')

    turnEffectFig.tight_layout()

    makeNestedPlotDirectory(analysisDir, 'effectOfTurn/', objecttype + sep)

    turnEffectFig.savefig(analysisDir + 'effectOfTurn/' + objecttype + sep
                          + FODataFile[0:-4] + '_turnHeadingChange_trial' + str(trial) + '.pdf', format='pdf')

    # Directional modulation of runs (Gomez-Marin and Louis, 2014) .....................................................

    turnModfig = plt.figure(figsize=(12, 7))
    turnMod = turnModfig.add_subplot(111)
    turnModsc = turnMod.scatter(gammaFull[selectedRangeDist], vRotFilt_ds[selectedRangeDist], marker='o', s=15,
                                linewidths=0, c=objDistance[selectedRangeDist], cmap=plt.cm.Spectral, alpha=0.5)
    turnMod.scatter(gammaFull[selectedRangeDistTurn], vRotFilt_ds[selectedRangeDistTurn], marker='o', s=15,
                    linewidths=0.5, c=objDistance[selectedRangeDistTurn], cmap=plt.cm.Spectral, alpha=0.5)
    turnMod.set_xlim(-np.pi, np.pi)
    turnMod.set_ylim(-5, 5)
    turnMod.set_xlabel('relative heading [rad]\n(or local bearing)')
    turnMod.set_ylabel('instantaneous rot. velocity [rad/s]\n(filtered)')

    turnModcb = plt.colorbar(turnModsc)
    turnModcb.set_label('distance from object [mm]')
    myAxisTheme(turnMod)

    turnMod.set_title(flyID + ' in VR (within range ' + str(near) + '-' + str(far) + 'mm)',  fontsize=13)

    makeNestedPlotDirectory(analysisDir, 'headingVsRotation/', objecttype + sep)

    turnModfig.savefig(analysisDir + 'headingVsRotation/' + objecttype + sep
                       + FODataFile[0:-4] + '_headingVsRotation_trial' + str(trial) + '.pdf', format='pdf')

    # ------------------------------------------------------------------------------------------------------------------
    # Save position and velocities for future analysis
    # ------------------------------------------------------------------------------------------------------------------

    toSave = {'time': time_ds,
              'xPos': xPos_ds,
              'yPos': yPos_ds,
              'xPosInMiniarena': xPosMA_ds,
              'yPosInMiniarena': yPosMA_ds,
              'headingAngle': angle_ds,
              'rotVelo': vRot_ds,
              'transVelo': vTrans_ds,
              'objectDistance': objDistance,
              'gamma': gamma}
    # Alternatively use numpy array:
    # toSave=np.vstack((time_ds,xPos_ds, yPos_ds, xPosMA_ds, yPosMA_ds, angle_ds, vRot_ds,vTrans_ds,objDistance,gamma))

    # Save data in this format as *.npy for easy loading..
    np.save(expDir + FODataFile[:-4], toSave)

    # ------------------------------------------------------------------------------------------------------------------
    print('\n \n Analysis ran successfully. \n \n')
    # ------------------------------------------------------------------------------------------------------------------

    return 0

