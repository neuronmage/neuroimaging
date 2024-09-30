#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.81.02), Sat Jan 17 21:43:45 2015
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
from math import ceil

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'wmaze.py'
expInfo = {'participant':'', 'scanner': 1}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + 'wmazebl_data/%s_%s_%s' %(expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
#save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

#Our own log file
datFile=open('wmazebl_data' + os.sep + '%s_wmazebl_%s.txt' %(expInfo['participant'], expInfo['date']),'a')
datFile.write('Run\tTrial\tStim\tStimOnset\tTrialType\tResp\tFBOnset\tCorrect\tIncorrect\tRT\n')

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(size=(1280, 800), fullscr=True, screen=0, allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-1,-1,-1], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    )
# store frame rate of monitor if we can measure it successfully
expInfo['frameRate']=win.getActualFrameRate()
if expInfo['frameRate']!=None:
    frameDur = 1.0/round(expInfo['frameRate'])
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

if expInfo['scanner'] == 1:
    wait_screen_text = 'Waiting for scanner to start.'
else:
    wait_screen_text = 'Press the return key to start the experiment.'

def BaselineTrial(bl_dict):
    # Draw the noise background
    blbkgrnd_image = (np.random.rand(2560,1600) > 0.85) * 1

    target_change = 0.85 - bl_dict['talpha']
    foil_change = 0.85 - bl_dict['falpha']
    
    bltarget_image = (np.random.rand(2560,1600) > target_change) * 1
    blfoil_image = (np.random.rand(2560,1600) > foil_change) * 1
    
    targ_pos = np.random.permutation(2)
    
    if targ_pos[0] == 0: #0 on the left side of screen
        bl_targ_pos = [-0.5, 0]
        bl_foil_pos = [0.5, 0]
    else:
        bl_targ_pos = [0.5, 0]
        bl_foil_pos = [-0.5, 0]
    return bl_targ_pos, bl_foil_pos, blbkgrnd_image, bltarget_image, blfoil_image
    
def QueueBaselineHist(bl_dict, correct):

    if len(bl_dict['bl_hist']) < 10:
        bl_dict['bl_hist'].append(correct)
        return bl_dict
    
    bl_dict['bl_hist'][0:8] = bl_dict['bl_hist'][1:9]
    bl_dict['bl_hist'][9] = correct
    
    if sum(bl_dict['bl_hist']) > 7:
        bl_dict['talpha'] = bl_dict['talpha'] * 0.9
    elif sum(bl_dict['bl_hist']) < 5:
        bl_dict['talpha'] = bl_dict['talpha'] * 1.1
        
    if bl_dict['talpha'] <= bl_dict['falpha']:
        bl_dict['talpha'] = bl_dict['falpha'] + 0.05
    
    return bl_dict

#Initialise components for Routine "waiting"
waiting_textClock=core.Clock()
waiting_text=visual.TextStim(win=win, ori=0, name='waiting_text',
    text=wait_screen_text,
    font='Arial',
    pos=[0, 0], height=0.1,wrapWidth=None,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()
isi_fix = visual.TextStim(win=win, ori=0, name='isi_fix',
    text=u'+',    font=u'Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
    color=u'white', colorSpace='rgb', opacity=1,
    depth=-1.0)
kaleido = visual.ImageStim(win=win, name='kaleido',
    image=None, mask=None,
    ori=0, pos=[0, 0], size=[0.5, 0.5],
    color=[1, 1, 1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
resp_box1 = visual.Rect(win=win, name='resp_box1',
    width=[0.25, 0.25][0], height=[0.25, 0.25][1],
    ori=0, pos=[-0.5, 0],
    lineWidth=1, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=None, fillColorSpace=None,
    opacity=1,interpolate=True)
resp_box2 = visual.Rect(win=win, name='resp_box2',
    width=[0.25, 0.25][0], height=[0.25, 0.25][1],
    ori=0, pos=[0.5, 0],
    lineWidth=1, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=None, fillColorSpace=None,
    opacity=1,interpolate=True)
feedback_box = visual.Rect(win=win, name='feedback_box',
    width=[0.25, 0.25][0], height=[0.25, 0.25][1],
    ori=0, pos=[0,0],
    lineWidth=1, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[1, 1, 1], fillColorSpace='rgb',
    opacity=1, interpolate=True)
feedback_text = visual.TextStim(win=win, ori=0, name='feedback_text',
    text=u'Any text\n\nincluding line breaks',    font=u'Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
    color=u'white', colorSpace='rgb', opacity=1,
    depth=-5.0)
bl_bkgrnd = visual.ImageStim(win=win, name='bl_bkgrnd',
    image=None, mask=None,
    ori=0, pos=[0, 0], size=[2, 2],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
bl_target = visual.ImageStim(win=win, name='bl_target',
    image=None, mask=None,
    ori=0, pos=[-0.5, 0], size=[0.3, 0.3],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
bl_foil = visual.ImageStim(win=win, name='bl_foil',
    image=None, mask=None,
    ori=0, pos=[0.5, 0], size=[0.3, 0.3],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)

all_kaleidos = os.listdir(os.getcwd() + '/' + 'kaleidojpgs')

bl_dict = {'talpha': 0.3, 'falpha': 0.2, 'bl_hist': []}

curr_run_kaleidos = []
for i in range(3):
    curr_run_kaleidos.append(all_kaleidos.pop(int(ceil(random(1) * len(all_kaleidos)))))

kaleido_stim = {'stimA':os.getcwd() + '/' + 'kaleidojpgs' + '/' + curr_run_kaleidos[0],'stimB':os.getcwd() + '/' + 'kaleidojpgs' + '/' + curr_run_kaleidos[1],'stimC':os.getcwd() + '/' + 'kaleidojpgs' + '/' + curr_run_kaleidos[2]}

curr_trls = np.array([0, 1, 2])
for i in range(80):
    if i == 0:
        tot_trials = curr_trls
    else:
        curr_trls = np.random.permutation(curr_trls)
        if tot_trials[-1] == 1 and curr_trls[0] == 1:
            while curr_trls[0] == 1:
                curr_trls = np.random.permutation(curr_trls)
        tot_trials = np.concatenate((tot_trials, curr_trls))
        
corr_resps = []
last_trial = []
for counter, i in enumerate(tot_trials):
    if counter == 0 and i == 0:
        corr_resps.append('2')
        last_trial.append(tot_trials[counter])
    else:
        if i == 0:
            corr_resps.append('2')
            last_trial.append(tot_trials[counter])
        elif i == 2:
            corr_resps.append('1')
            last_trial.append(tot_trials[counter])
        elif i == 1 and last_trial[-1] == 0:
            corr_resps.append('1')
            last_trial.append(tot_trials[counter])
        elif i == 1 and last_trial[-1] == 2:
            corr_resps.append('2')
            last_trial.append(tot_trials[counter])

corr_resps = np.array(corr_resps)
num_trls = len(tot_trials)
trls_firsthalf = tot_trials[:0.5*num_trls]
corr_resps_firsthalf = corr_resps[:0.5*num_trls]
trls_secondhalf = tot_trials[0.5*num_trls:]
corr_resps_secondhalf = corr_resps[0.5*num_trls:]

while len(trls_firsthalf) < 160:
    # Determine random position
    rand_pos = np.random.permutation(len(trls_firsthalf))
    if rand_pos[0] != 0:
        trls_firsthalf = np.insert(trls_firsthalf, rand_pos[0], 3)
        corr_resps_firsthalf = np.insert(corr_resps_firsthalf, rand_pos[0], 3)
while len(trls_secondhalf) < 160:
    # Determine random position
    rand_pos = np.random.permutation(len(trls_secondhalf))
    if rand_pos[0] != 0:
        trls_secondhalf = np.insert(trls_secondhalf, rand_pos[0], 3)
        corr_resps_secondhalf = np.insert(corr_resps_secondhalf, rand_pos[0], 3)
            
# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 
begExpClock=core.Clock() #to track the time since the actual paradigm started (after the scanner screen)

trls_1stand2ndhalf = [trls_firsthalf, trls_secondhalf]
resps_1stand2ndhalf = [corr_resps_firsthalf, corr_resps_secondhalf]

num_runs = [1]
myarray_runs = []
for i in range(1):
    myarray_runs.append({'run': num_runs[i], 'run_trials': trls_1stand2ndhalf[i], 'run_resps': resps_1stand2ndhalf[i] }) #puts data into an array of dictionaries that the TrialHandler function will accept
runs=data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=None,
    trialList=myarray_runs,
    seed=None, name='runs')
thisExp.addLoop(runs)#add the loop to the experiment
thisRun=runs.trialList[0]#so we can initialise stimuli with some values
#abbreviate parameter names if possible (e.g. rgb=thisTrial.rgb)
if thisRun!=None:
    for paramName in thisRun.keys():
        exec(paramName+'=thisRun.'+paramName)
        
for thisRun in runs:
    if thisRun!=None:
        for paramName in thisRun.keys():
            exec(paramName+'=thisRun.'+paramName)
        
    #------Prepare to start Routine "waiting"-------
    t = 0
    waiting_textClock.reset()  # clock 
    frameN = -1
    # update component parameters for each repeat
    sync_pulse = event.BuilderKeyResponse() #create an object of type KeyResponse
    sync_pulse.status=NOT_STARTED
    # keep track of which components have finished
    waiting_textComponents = []
    waiting_textComponents.append(waiting_text)
    for thisComponent in waiting_textComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED

    #-------Start Routine "waiting"-------
    continueRoutine = True
    while continueRoutine:
        # get current time
        t = waiting_textClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *Waiting* updates
        if t >= 0.0 and waiting_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            waiting_text.tStart = t  # underestimates by a little under one frame
            waiting_text.frameNStart = frameN  # exact frame index
            waiting_text.setAutoDraw(True)
        
        #*sync_pulse* updates
        if t>=0.0 and sync_pulse.status==NOT_STARTED:
            #keep track of start time/frame for later
            sync_pulse.tStart=t#underestimates by a little under one frame
            sync_pulse.frameNStart=frameN#exact frame index
            sync_pulse.status=STARTED
            #keyboard checking is just starting
            sync_pulse.clock.reset() # now t=0
            event.clearEvents()
        if sync_pulse.status==STARTED:#only update if being drawn
            theseKeys = event.getKeys(keyList=['return', '5'])
            if len(theseKeys)>0:#at least one key was pressed
                sync_pulse.keys=theseKeys[-1]#just the last key pressed
                sync_pulse.rt = sync_pulse.clock.getTime()
                #abort routine on response
                core.wait(10)
                continueRoutine=False
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineTimer.reset()  # if we abort early the non-slip timer needs reset
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in waiting_textComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
        else:  # this Routine was not non-slip safe so reset non-slip timer
            routineTimer.reset()

    #-------Ending Routine "waiting"-------
    for thisComponent in waiting_textComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
            
    myarray = []
    for i in range(len(run_trials)):
        myarray.append({'trial_type': run_trials[i], 'correct_resp': run_resps[i]}) #puts data into an array of dictionaries that the TrialHandler function will accept
    trials=data.TrialHandler(nReps=1, method='sequential', 
        extraInfo=expInfo, originPath=None,
        trialList=myarray,
        seed=None, name='trials')
    thisExp.addLoop(trials)#add the loop to the experiment
    thisTrial=trials.trialList[0]#so we can initialise stimuli with some values
    #abbreviate parameter names if possible (e.g. rgb=thisTrial.rgb)
    if thisTrial!=None:
        for paramName in thisTrial.keys():
            exec(paramName+'=thisTrial.'+paramName)

    begExpClock.reset()
    for thisTrial in trials:
        currentLoop = trials
        #abbrieviate parameter names if possible (e.g. rgb=thisTrial.rgb)
        if thisTrial!=None:
            for paramName in thisTrial.keys():
                exec(paramName+'=thisTrial.'+paramName)

        #------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock 
        frameN = -1
        routineTimer.add(2.500000)
        # update component parameters for each repeat
        key_resp = event.BuilderKeyResponse()  # create an object of type KeyResponse
        key_resp.status = NOT_STARTED
        # keep track of which components have finished
        trialComponents = []
        trialComponents.append(isi_fix)
        trialComponents.append(kaleido)
        trialComponents.append(resp_box1)
        trialComponents.append(resp_box2)
        trialComponents.append(key_resp)
        trialComponents.append(feedback_box)
        trialComponents.append(feedback_text)
        trialComponents.append(bl_bkgrnd)
        trialComponents.append(bl_target)
        trialComponents.append(bl_foil)
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        if trial_type == 0:
            kaleido.image = kaleido_stim['stimA']
            Trial_Type = 'A'
        elif trial_type == 1:
            kaleido.image = kaleido_stim['stimB']
            Trial_Type = 'B'
        elif trial_type == 2:
            kaleido.image = kaleido_stim['stimC']
            Trial_Type = 'C'
        else:
            bl_targ_pos, bl_foil_pos, blbkgrnd_image, bltarget_image, blfoil_image = BaselineTrial(bl_dict)
            bl_bkgrnd.image = blbkgrnd_image
            bl_target.image = bltarget_image
            bl_target.pos = bl_targ_pos
            bl_foil.image = blfoil_image
            bl_foil.pos = bl_foil_pos
            Trial_Type = 'BL'
        
        #-------Start Routine "trial"-------
        continueRoutine = True
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *isi_fix* updates
            if t >= 0.0 and isi_fix.status == NOT_STARTED:
                # keep track of start time/frame for later
                isi_fix.tStart = t  # underestimates by a little under one frame
                isi_fix.frameNStart = frameN  # exact frame index
                isi_fix.setAutoDraw(True)
            elif isi_fix.status == STARTED and t >= (0.0 + (0.7-win.monitorFramePeriod*0.75)): #most of one frame period left
                isi_fix.setAutoDraw(False)
            
            if trial_type == 3:
                if t >= 0.7 and bl_bkgrnd.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    bl_bkgrnd.tStart = t  # underestimates by a little under one frame
                    bl_bkgrnd.frameNStart = frameN  # exact frame index
                    bl_bkgrnd.setAutoDraw(True)
                    trialTOnset=begExpClock.getTime()
                elif bl_bkgrnd.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    bl_bkgrnd.setAutoDraw(False)
                    
                if t >= 0.7 and bl_target.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    bl_target.tStart = t  # underestimates by a little under one frame
                    bl_target.frameNStart = frameN  # exact frame index
                    bl_target.setAutoDraw(True)
                elif bl_target.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    bl_target.setAutoDraw(False)
                    
                if t >= 0.7 and bl_foil.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    bl_foil.tStart = t  # underestimates by a little under one frame
                    bl_foil.frameNStart = frameN  # exact frame index
                    bl_foil.setAutoDraw(True)
                elif bl_foil.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    bl_foil.setAutoDraw(False)
            
                # *key_resp* updates
                if t >= 0.7 and key_resp.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    key_resp.tStart = t  # underestimates by a little under one frame
                    key_resp.frameNStart = frameN  # exact frame index
                    key_resp.status = STARTED
                    # keyboard checking is just starting
                    key_resp.clock.reset()  # now t=0
                    event.clearEvents(eventType='keyboard')
                elif key_resp.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    key_resp.status = STOPPED
                if key_resp.status == STARTED:
                    theseKeys = event.getKeys(keyList=['1', '2'])
                    
                    # check for quit:
                    if "escape" in theseKeys:
                        endExpNow = True
                    if len(theseKeys) > 0:  # at least one key was pressed
                        key_resp.keys = theseKeys[-1]  # just the last key pressed
                        key_resp.rt = key_resp.clock.getTime()
                        # a response ends the routine
                        #continueRoutine = True
                        
                        if key_resp.keys == '1':
                            feedback_box.pos = [-0.5, 0]
                            if bl_target.pos[0] == -0.5:
                                feedback_text.text = "Yes!"
                                feedback_text.color = u'green'
                                correct, incorrect = [1, 0]
                                bl_dict = QueueBaselineHist(bl_dict, correct)
                            else:
                                feedback_text.text = "No!"
                                feedback_text.color = u'red'
                                correct, incorrect = [0, 1]
                                bl_dict = QueueBaselineHist(bl_dict, correct)
                        elif key_resp.keys == '2':
                            feedback_box.pos = [0.5, 0]
                            if bl_target.pos[0] == 0.5:
                                feedback_text.text = "Yes!"
                                feedback_text.color = u'green'
                                correct, incorrect = [1, 0]
                                bl_dict = QueueBaselineHist(bl_dict, correct)
                            else:
                                feedback_text.text = "No!"
                                feedback_text.color = u'red'
                                correct, incorrect = [0, 1]
                                bl_dict = QueueBaselineHist(bl_dict, correct)
                        
                        # *feedback_box* updates
                        if t >= 0.7 and feedback_box.status == NOT_STARTED:
                            # keep track of start time/frame for later
                            feedback_box.tStart = t
                            feedback_box.frameNStart = frameN
                            feedback_box.setAutoDraw(True)
                        elif feedback_box.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                            feedback_box.setAutoDraw(False)
            else:
                # *kaleido* updates
                if t >= 0.7 and kaleido.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    kaleido.tStart = t  # underestimates by a little under one frame
                    kaleido.frameNStart = frameN  # exact frame index
                    kaleido.setAutoDraw(True)
                    trialTOnset=begExpClock.getTime()
                elif kaleido.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    kaleido.setAutoDraw(False)
                    
                # *resp_box1* updates
                if t >= 0.7 and resp_box1.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    resp_box1.tStart = t  # underestimates by a little under one frame
                    resp_box1.frameNStart = frameN  # exact frame index
                    resp_box1.setAutoDraw(True)
                elif resp_box1.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    resp_box1.setAutoDraw(False)
                
                # *resp_box2* updates
                if t >= 0.7 and resp_box2.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    resp_box2.tStart = t  # underestimates by a little under one frame
                    resp_box2.frameNStart = frameN  # exact frame index
                    resp_box2.setAutoDraw(True)
                elif resp_box2.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    resp_box2.setAutoDraw(False)
                    
                # *key_resp* updates
                if t >= 0.7 and key_resp.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    key_resp.tStart = t  # underestimates by a little under one frame
                    key_resp.frameNStart = frameN  # exact frame index
                    key_resp.status = STARTED
                    # keyboard checking is just starting
                    key_resp.clock.reset()  # now t=0
                    event.clearEvents(eventType='keyboard')
                elif key_resp.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    key_resp.status = STOPPED
                if key_resp.status == STARTED:
                    theseKeys = event.getKeys(keyList=['1', '2'])
                    
                    # check for quit:
                    if "escape" in theseKeys:
                        endExpNow = True
                    if len(theseKeys) > 0:  # at least one key was pressed
                        key_resp.keys = theseKeys[-1]  # just the last key pressed
                        key_resp.rt = key_resp.clock.getTime()
                        # a response ends the routine
                        #continueRoutine = True
                        
                        if key_resp.keys == '1':
                            feedback_box.pos = [-0.5, 0]
                            if correct_resp == '1':
                                feedback_text.text = "Yes!"
                                feedback_text.color = u'green'
                                correct, incorrect = [1, 0]
                            else:
                                feedback_text.text = "No!"
                                feedback_text.color = u'red'
                                correct, incorrect = [0, 1]
                        elif key_resp.keys == '2':
                            feedback_box.pos = [0.5, 0]
                            if correct_resp == '2':
                                feedback_text.text = "Yes!"
                                feedback_text.color = u'green'
                                correct, incorrect = [1, 0]
                            else:
                                feedback_text.text = "No!"
                                feedback_text.color = u'red'
                                correct, incorrect = [0, 1]
                        
                        # *feedback_box* updates
                        if t >= 0.7 and feedback_box.status == NOT_STARTED:
                            # keep track of start time/frame for later
                            feedback_box.tStart = t
                            feedback_box.frameNStart = frameN
                            feedback_box.setAutoDraw(True)
                        elif feedback_box.status == STARTED and t >= (0.7 + (1.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                            feedback_box.setAutoDraw(False)
            
            if key_resp.keys in ['', [], None]:
                feedback_text.text = "?"
                feedback_text.color = u'white'
                correct, incorrect = ['NaN', 'NaN']
                key_resp.rt = 'NaN'
                key_resp.keys = 'NR'
            
            # *feedback_text* updates
            if t >= 1.7 and feedback_text.status == NOT_STARTED:
                # keep track of start time/frame for later
                feedback_text.tStart = t  # underestimates by a little under one frame
                feedback_text.frameNStart = frameN  # exact frame index
                feedback_text.setAutoDraw(True)
                FBTOnset=begExpClock.getTime()
            elif feedback_text.status == STARTED and t >= (1.7 + (0.8-win.monitorFramePeriod*0.75)): #most of one frame period left
                feedback_text.setAutoDraw(False)

            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineTimer.reset()  # if we abort early the non-slip timer needs reset
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # check for quit (the Esc key)
            if endExpNow or event.getKeys(keyList=["escape"]):
                core.quit()
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        #-------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
         
             
        datFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(run,trials.thisTrialN+1,trialTOnset,Trial_Type,key_resp.keys,FBTOnset,correct,incorrect,key_resp.rt))
    
            
win.close()
core.quit()
