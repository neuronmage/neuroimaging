#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.90.3),
    on Mon May  6 15:15:34 2019
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import absolute_import, division
from psychopy import locale_setup, sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'vCAT_localizer'  # from the Builder filename that created this script
expInfo = {u'participant': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=[1280, 800], fullscr=True, screen=0,
    allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-0.5,-0.5,-0.5], colorSpace='rgb',
    blendMode='avg', useFBO=True)
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Initialize components for Routine "instructions"
instructionsClock = core.Clock()
bkg_instructions = visual.ImageStim(
    win=win, name='bkg_instructions',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
text_instructions = visual.TextStim(win=win, name='text_instructions',
    text='During this part of the experiment, you will be asked two kinds of questions. \n\n#1 "Yes" or "No" questions:\nWhen you see a face: Is it female?\nWhen you see a scene: Is there a lake or river?\nPlease press the left button for "yes" and the right button for "no".\n\n#2 "Greater" or "Less" than 5:\nYou will occasionally be presented with a single number: \nPress the left button if the number is less than 5.\nPress the right button if the number is greater than 5.\n\nPlease wait for the scanner to begin.\n',
    font='Ariel',
    pos=(0, 0), height=0.08, wrapWidth=1.5, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-1.0);
expClock = core.Clock()

import os
import numpy as np
from numpy.random import random, randint, choice, permutation
import random
from glob import glob

stim_opacity = 0

##Set Dependencies
#scenes
localizer_scenes = glob(os.getcwd() + '/scene_img/*/*jpg')
#faces 
localizer_faces = glob(os.getcwd() + '/face_img/*/*bmp')

# Initialize components for Routine "delay_init"
delay_initClock = core.Clock()
delay = visual.ImageStim(
    win=win, name='delay',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
text_2 = visual.TextStim(win=win, name='text_2',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-1.0);

# Initialize components for Routine "stim_prep"
stim_prepClock = core.Clock()


# Initialize components for Routine "trial"
trialClock = core.Clock()

isi_bkg_1 = visual.ImageStim(
    win=win, name='isi_bkg_1',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
bkg_1 = visual.ImageStim(
    win=win, name='bkg_1',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-3.0)
fixation_1 = visual.TextStim(win=win, name='fixation_1',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-4.0);
stim_img_1 = visual.ImageStim(
    win=win, name='stim_img_1',
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(1, 1.25),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-5.0)
math_num = visual.TextStim(win=win, name='math_num',
    text='default text',
    font='Arial',
    pos=(0, 0), height=0.8, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1.0,
    depth=-6.0);

# Initialize components for Routine "stim_reset"
stim_resetClock = core.Clock()


# Initialize components for Routine "trial"
trialClock = core.Clock()

isi_bkg_1 = visual.ImageStim(
    win=win, name='isi_bkg_1',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
bkg_1 = visual.ImageStim(
    win=win, name='bkg_1',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-3.0)
fixation_1 = visual.TextStim(win=win, name='fixation_1',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-4.0);
stim_img_1 = visual.ImageStim(
    win=win, name='stim_img_1',
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(1, 1.25),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-5.0)
math_num = visual.TextStim(win=win, name='math_num',
    text='default text',
    font='Arial',
    pos=(0, 0), height=0.8, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1.0,
    depth=-6.0);

# Initialize components for Routine "delay_end"
delay_endClock = core.Clock()
image_delay = visual.ImageStim(
    win=win, name='image_delay',
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
text = visual.TextStim(win=win, name='text',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-1.0);

# Initialize components for Routine "goodbye"
goodbyeClock = core.Clock()
text_goodbye = visual.TextStim(win=win, name='text_goodbye',
    text='Great job!\n\nPlease wait for the next section of the experiment to begin.',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "instructions"-------
t = 0
instructionsClock.reset()  # clock
frameN = -1
continueRoutine = True
# update component parameters for each repeat
start_key = event.BuilderKeyResponse()

# keep track of which components have finished
instructionsComponents = [bkg_instructions, text_instructions, start_key]
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "instructions"-------
while continueRoutine:
    # get current time
    t = instructionsClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *bkg_instructions* updates
    if t >= 0.0 and bkg_instructions.status == NOT_STARTED:
        # keep track of start time/frame for later
        bkg_instructions.tStart = t
        bkg_instructions.frameNStart = frameN  # exact frame index
        bkg_instructions.setAutoDraw(True)
    
    # *text_instructions* updates
    if t >= 0.0 and text_instructions.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_instructions.tStart = t
        text_instructions.frameNStart = frameN  # exact frame index
        text_instructions.setAutoDraw(True)
    
    # *start_key* updates
    if t >= 0.0 and start_key.status == NOT_STARTED:
        # keep track of start time/frame for later
        start_key.tStart = t
        start_key.frameNStart = frameN  # exact frame index
        start_key.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(start_key.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if start_key.status == STARTED:
        theseKeys = event.getKeys(keyList=['space', '5'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            start_key.keys = theseKeys[-1]  # just the last key pressed
            start_key.rt = start_key.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "instructions"-------
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if start_key.keys in ['', [], None]:  # No response was made
    start_key.keys=None
thisExp.addData('start_key.keys',start_key.keys)
if start_key.keys != None:  # we had a response
    thisExp.addData('start_key.rt', start_key.rt)
thisExp.nextEntry()
expClock.reset()
# the Routine "instructions" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# ------Prepare to start Routine "delay_init"-------
t = 0
delay_initClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(5.750000)
# update component parameters for each repeat
# keep track of which components have finished
delay_initComponents = [delay, text_2]
for thisComponent in delay_initComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "delay_init"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = delay_initClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *delay* updates
    if t >= 0.0 and delay.status == NOT_STARTED:
        # keep track of start time/frame for later
        delay.tStart = t
        delay.frameNStart = frameN  # exact frame index
        delay.setAutoDraw(True)
    frameRemains = 0.0 + 5.75- win.monitorFramePeriod * 0.75  # most of one frame period left
    if delay.status == STARTED and t >= frameRemains:
        delay.setAutoDraw(False)
    
    # *text_2* updates
    if t >= 0.0 and text_2.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_2.tStart = t
        text_2.frameNStart = frameN  # exact frame index
        text_2.setAutoDraw(True)
    frameRemains = 0.0 + 5.75- win.monitorFramePeriod * 0.75  # most of one frame period left
    if text_2.status == STARTED and t >= frameRemains:
        text_2.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in delay_initComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "delay_init"-------
for thisComponent in delay_initComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# set up handler to look after randomisation of conditions etc
runs = data.TrialHandler(nReps=7, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='runs')
thisExp.addLoop(runs)  # add the loop to the experiment
thisRun = runs.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisRun.rgb)
if thisRun != None:
    for paramName in thisRun:
        exec('{} = thisRun[paramName]'.format(paramName))

for thisRun in runs:
    currentLoop = runs
    # abbreviate parameter names if possible (e.g. rgb = thisRun.rgb)
    if thisRun != None:
        for paramName in thisRun:
            exec('{} = thisRun[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "stim_prep"-------
    t = 0
    stim_prepClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    # update component parameters for each repeat
    current_run = 'faces'
    math_stim = ''
    # keep track of which components have finished
    stim_prepComponents = []
    for thisComponent in stim_prepComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "stim_prep"-------
    while continueRoutine:
        # get current time
        t = stim_prepClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stim_prepComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "stim_prep"-------
    for thisComponent in stim_prepComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    # the Routine "stim_prep" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    faces = data.TrialHandler(nReps=1, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions('localizer_schedule.xlsx'),
        seed=None, name='faces')
    thisExp.addLoop(faces)  # add the loop to the experiment
    thisFace = faces.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisFace.rgb)
    if thisFace != None:
        for paramName in thisFace:
            exec('{} = thisFace[paramName]'.format(paramName))
    
    for thisFace in faces:
        currentLoop = faces
        # abbreviate parameter names if possible (e.g. rgb = thisFace.rgb)
        if thisFace != None:
            for paramName in thisFace:
                exec('{} = thisFace[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock
        frameN = -1
        continueRoutine = True
        routineTimer.add(1.250000)
        # update component parameters for each repeat
        resp_trial_1 = event.BuilderKeyResponse()
        runs = eval(current_run)
        runs.addData("StimOnset", expClock.getTime() + 0.25) #add 0.25 to account for ISI
        exec_run = False
        
        options = range(1,5) + range(6,10) #all numbers but 5
        i = runs.thisN
        #localizer trial
        if trialtype == 0: 
            math_stim = ''
            mathstim_opacity = 0
            stim_opacity = 1
            if current_run == 'faces':
                stim = np.random.choice(localizer_faces)
            else:
                stim = np.random.choice(localizer_scenes)
            runs.addData('stim', stim.split("/")[-1])
        else:
            stim_opacity = 0
            mathstim_opacity = 1
            math_stim = np.random.choice(options)
            runs.addData('stim', math_stim)
        
        stim_img_1.setOpacity(stim_opacity)
        stim_img_1.setImage(stim)
        math_num.setOpacity(mathstim_opacity)
        math_num.setText(math_stim)
        # keep track of which components have finished
        trialComponents = [resp_trial_1, isi_bkg_1, bkg_1, fixation_1, stim_img_1, math_num]
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        # -------Start Routine "trial"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *resp_trial_1* updates
            if t >= 0.25 and resp_trial_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                resp_trial_1.tStart = t
                resp_trial_1.frameNStart = frameN  # exact frame index
                resp_trial_1.status = STARTED
                # keyboard checking is just starting
                win.callOnFlip(resp_trial_1.clock.reset)  # t=0 on next screen flip
                event.clearEvents(eventType='keyboard')
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if resp_trial_1.status == STARTED and t >= frameRemains:
                resp_trial_1.status = STOPPED
            if resp_trial_1.status == STARTED:
                theseKeys = event.getKeys(keyList=['1', '2'])
                
                # check for quit:
                if "escape" in theseKeys:
                    endExpNow = True
                if len(theseKeys) > 0:  # at least one key was pressed
                    if resp_trial_1.keys == []:  # then this was the first keypress
                        resp_trial_1.keys = theseKeys[0]  # just the first key pressed
                        resp_trial_1.rt = resp_trial_1.clock.getTime()
            
            
            # *isi_bkg_1* updates
            if t >= 0.0 and isi_bkg_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                isi_bkg_1.tStart = t
                isi_bkg_1.frameNStart = frameN  # exact frame index
                isi_bkg_1.setAutoDraw(True)
            frameRemains = 0.0 + 0.3- win.monitorFramePeriod * 0.75  # most of one frame period left
            if isi_bkg_1.status == STARTED and t >= frameRemains:
                isi_bkg_1.setAutoDraw(False)
            
            # *bkg_1* updates
            if t >= 0.25 and bkg_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                bkg_1.tStart = t
                bkg_1.frameNStart = frameN  # exact frame index
                bkg_1.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if bkg_1.status == STARTED and t >= frameRemains:
                bkg_1.setAutoDraw(False)
            
            # *fixation_1* updates
            if t >= 0.0 and fixation_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                fixation_1.tStart = t
                fixation_1.frameNStart = frameN  # exact frame index
                fixation_1.setAutoDraw(True)
            frameRemains = 0.0 + 0.25- win.monitorFramePeriod * 0.75  # most of one frame period left
            if fixation_1.status == STARTED and t >= frameRemains:
                fixation_1.setAutoDraw(False)
            
            # *stim_img_1* updates
            if t >= 0.25 and stim_img_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                stim_img_1.tStart = t
                stim_img_1.frameNStart = frameN  # exact frame index
                stim_img_1.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if stim_img_1.status == STARTED and t >= frameRemains:
                stim_img_1.setAutoDraw(False)
            
            # *math_num* updates
            if t >= 0.25 and math_num.status == NOT_STARTED:
                # keep track of start time/frame for later
                math_num.tStart = t
                math_num.frameNStart = frameN  # exact frame index
                math_num.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if math_num.status == STARTED and t >= frameRemains:
                math_num.setAutoDraw(False)
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
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
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # check responses
        if resp_trial_1.keys in ['', [], None]:  # No response was made
            resp_trial_1.keys=None
        faces.addData('resp_trial_1.keys',resp_trial_1.keys)
        if resp_trial_1.keys != None:  # we had a response
            faces.addData('resp_trial_1.rt', resp_trial_1.rt)
        correct = None
        if mathstim_opacity == 0:
            if stim.split("/")[-2] == 'face_img':
                if (stim.split("/")[-1] == 'female' and resp_trial_1.keys == "1") or (stim.split("/")[-1] == 'male' and resp_trial_1.keys == "2"):
                    correct = 1
                else:
                    print type(resp_trial_1.keys)
                    correct = 0
            else:
                if (stim.split("/")[-1] == 'water' and resp_trial_1.keys == "1") or (stim.split("/")[-1] == 'land' and resp_trial_1.keys == "2"):
                    correct = 1
                else:
                    correct = 0
        else:
            if (math_stim < 5 and resp_trial_1.keys == "1") or (math_stim > 5 and resp_trial_1.keys == "2"):
                correct = 1
            else:
                correct = 0
        
        runs.addData('trial_acc', correct)
        thisExp.nextEntry()
        
    # completed 1 repeats of 'faces'
    
    
    # ------Prepare to start Routine "stim_reset"-------
    t = 0
    stim_resetClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    # update component parameters for each repeat
    current_run = 'scenes'
    
    
    # keep track of which components have finished
    stim_resetComponents = []
    for thisComponent in stim_resetComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "stim_reset"-------
    while continueRoutine:
        # get current time
        t = stim_resetClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stim_resetComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "stim_reset"-------
    for thisComponent in stim_resetComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    # the Routine "stim_reset" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    scenes = data.TrialHandler(nReps=1, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions('localizer_schedule.xlsx'),
        seed=None, name='scenes')
    thisExp.addLoop(scenes)  # add the loop to the experiment
    thisScene = scenes.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisScene.rgb)
    if thisScene != None:
        for paramName in thisScene:
            exec('{} = thisScene[paramName]'.format(paramName))
    
    for thisScene in scenes:
        currentLoop = scenes
        # abbreviate parameter names if possible (e.g. rgb = thisScene.rgb)
        if thisScene != None:
            for paramName in thisScene:
                exec('{} = thisScene[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock
        frameN = -1
        continueRoutine = True
        routineTimer.add(1.250000)
        # update component parameters for each repeat
        resp_trial_1 = event.BuilderKeyResponse()
        runs = eval(current_run)
        runs.addData("StimOnset", expClock.getTime() + 0.25) #add 0.25 to account for ISI
        exec_run = False
        
        options = range(1,5) + range(6,10) #all numbers but 5
        i = runs.thisN
        #localizer trial
        if trialtype == 0: 
            math_stim = ''
            mathstim_opacity = 0
            stim_opacity = 1
            if current_run == 'faces':
                stim = np.random.choice(localizer_faces)
            else:
                stim = np.random.choice(localizer_scenes)
            runs.addData('stim', stim.split("/")[-1])
        else:
            stim_opacity = 0
            mathstim_opacity = 1
            math_stim = np.random.choice(options)
            runs.addData('stim', math_stim)
        
        stim_img_1.setOpacity(stim_opacity)
        stim_img_1.setImage(stim)
        math_num.setOpacity(mathstim_opacity)
        math_num.setText(math_stim)
        # keep track of which components have finished
        trialComponents = [resp_trial_1, isi_bkg_1, bkg_1, fixation_1, stim_img_1, math_num]
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        # -------Start Routine "trial"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *resp_trial_1* updates
            if t >= 0.25 and resp_trial_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                resp_trial_1.tStart = t
                resp_trial_1.frameNStart = frameN  # exact frame index
                resp_trial_1.status = STARTED
                # keyboard checking is just starting
                win.callOnFlip(resp_trial_1.clock.reset)  # t=0 on next screen flip
                event.clearEvents(eventType='keyboard')
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if resp_trial_1.status == STARTED and t >= frameRemains:
                resp_trial_1.status = STOPPED
            if resp_trial_1.status == STARTED:
                theseKeys = event.getKeys(keyList=['1', '2'])
                
                # check for quit:
                if "escape" in theseKeys:
                    endExpNow = True
                if len(theseKeys) > 0:  # at least one key was pressed
                    if resp_trial_1.keys == []:  # then this was the first keypress
                        resp_trial_1.keys = theseKeys[0]  # just the first key pressed
                        resp_trial_1.rt = resp_trial_1.clock.getTime()
            
            
            # *isi_bkg_1* updates
            if t >= 0.0 and isi_bkg_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                isi_bkg_1.tStart = t
                isi_bkg_1.frameNStart = frameN  # exact frame index
                isi_bkg_1.setAutoDraw(True)
            frameRemains = 0.0 + 0.3- win.monitorFramePeriod * 0.75  # most of one frame period left
            if isi_bkg_1.status == STARTED and t >= frameRemains:
                isi_bkg_1.setAutoDraw(False)
            
            # *bkg_1* updates
            if t >= 0.25 and bkg_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                bkg_1.tStart = t
                bkg_1.frameNStart = frameN  # exact frame index
                bkg_1.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if bkg_1.status == STARTED and t >= frameRemains:
                bkg_1.setAutoDraw(False)
            
            # *fixation_1* updates
            if t >= 0.0 and fixation_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                fixation_1.tStart = t
                fixation_1.frameNStart = frameN  # exact frame index
                fixation_1.setAutoDraw(True)
            frameRemains = 0.0 + 0.25- win.monitorFramePeriod * 0.75  # most of one frame period left
            if fixation_1.status == STARTED and t >= frameRemains:
                fixation_1.setAutoDraw(False)
            
            # *stim_img_1* updates
            if t >= 0.25 and stim_img_1.status == NOT_STARTED:
                # keep track of start time/frame for later
                stim_img_1.tStart = t
                stim_img_1.frameNStart = frameN  # exact frame index
                stim_img_1.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if stim_img_1.status == STARTED and t >= frameRemains:
                stim_img_1.setAutoDraw(False)
            
            # *math_num* updates
            if t >= 0.25 and math_num.status == NOT_STARTED:
                # keep track of start time/frame for later
                math_num.tStart = t
                math_num.frameNStart = frameN  # exact frame index
                math_num.setAutoDraw(True)
            frameRemains = 0.25 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
            if math_num.status == STARTED and t >= frameRemains:
                math_num.setAutoDraw(False)
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
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
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # check responses
        if resp_trial_1.keys in ['', [], None]:  # No response was made
            resp_trial_1.keys=None
        scenes.addData('resp_trial_1.keys',resp_trial_1.keys)
        if resp_trial_1.keys != None:  # we had a response
            scenes.addData('resp_trial_1.rt', resp_trial_1.rt)
        correct = None
        if mathstim_opacity == 0:
            if stim.split("/")[-2] == 'face_img':
                if (stim.split("/")[-1] == 'female' and resp_trial_1.keys == "1") or (stim.split("/")[-1] == 'male' and resp_trial_1.keys == "2"):
                    correct = 1
                else:
                    print type(resp_trial_1.keys)
                    correct = 0
            else:
                if (stim.split("/")[-1] == 'water' and resp_trial_1.keys == "1") or (stim.split("/")[-1] == 'land' and resp_trial_1.keys == "2"):
                    correct = 1
                else:
                    correct = 0
        else:
            if (math_stim < 5 and resp_trial_1.keys == "1") or (math_stim > 5 and resp_trial_1.keys == "2"):
                correct = 1
            else:
                correct = 0
        
        runs.addData('trial_acc', correct)
        thisExp.nextEntry()
        
    # completed 1 repeats of 'scenes'
    
    thisExp.nextEntry()
    
# completed 7 repeats of 'runs'


# ------Prepare to start Routine "delay_end"-------
t = 0
delay_endClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(6.000000)
# update component parameters for each repeat
# keep track of which components have finished
delay_endComponents = [image_delay, text]
for thisComponent in delay_endComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "delay_end"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = delay_endClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *image_delay* updates
    if t >= 0.0 and image_delay.status == NOT_STARTED:
        # keep track of start time/frame for later
        image_delay.tStart = t
        image_delay.frameNStart = frameN  # exact frame index
        image_delay.setAutoDraw(True)
    frameRemains = 0.0 + 6.0- win.monitorFramePeriod * 0.75  # most of one frame period left
    if image_delay.status == STARTED and t >= frameRemains:
        image_delay.setAutoDraw(False)
    
    # *text* updates
    if t >= 0.0 and text.status == NOT_STARTED:
        # keep track of start time/frame for later
        text.tStart = t
        text.frameNStart = frameN  # exact frame index
        text.setAutoDraw(True)
    frameRemains = 0.0 + 6- win.monitorFramePeriod * 0.75  # most of one frame period left
    if text.status == STARTED and t >= frameRemains:
        text.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in delay_endComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "delay_end"-------
for thisComponent in delay_endComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# ------Prepare to start Routine "goodbye"-------
t = 0
goodbyeClock.reset()  # clock
frameN = -1
continueRoutine = True
# update component parameters for each repeat
key_resp_2 = event.BuilderKeyResponse()
# keep track of which components have finished
goodbyeComponents = [text_goodbye, key_resp_2]
for thisComponent in goodbyeComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "goodbye"-------
while continueRoutine:
    # get current time
    t = goodbyeClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_goodbye* updates
    if t >= 0.0 and text_goodbye.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_goodbye.tStart = t
        text_goodbye.frameNStart = frameN  # exact frame index
        text_goodbye.setAutoDraw(True)
    
    # *key_resp_2* updates
    if t >= 0.0 and key_resp_2.status == NOT_STARTED:
        # keep track of start time/frame for later
        key_resp_2.tStart = t
        key_resp_2.frameNStart = frameN  # exact frame index
        key_resp_2.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if key_resp_2.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            key_resp_2.keys = theseKeys[-1]  # just the last key pressed
            key_resp_2.rt = key_resp_2.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in goodbyeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "goodbye"-------
for thisComponent in goodbyeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if key_resp_2.keys in ['', [], None]:  # No response was made
    key_resp_2.keys=None
thisExp.addData('key_resp_2.keys',key_resp_2.keys)
if key_resp_2.keys != None:  # we had a response
    thisExp.addData('key_resp_2.rt', key_resp_2.rt)
thisExp.nextEntry()
# the Routine "goodbye" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()





# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
