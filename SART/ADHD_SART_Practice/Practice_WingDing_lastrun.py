#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.83.01), Thu Feb 11 11:47:25 2016
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys # to get file system encoding

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# Store info about the experiment session
expName = u'Practice_SART_RTV'  # from the Builder filename that created this script
expInfo = {u'zzzzzz': u'', u'participant': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' %(expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=u'/Users/amandarenfro/Google Drive/FIU - PhD/MaD Lab/ADHD_RTV/ADHD_SART_Practice/Practice_WingDing.psyexp',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(size=(1280, 800), fullscr=True, screen=0, allowGUI=False, allowStencil=False,
    monitor=u'testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    )
# store frame rate of monitor if we can measure it successfully
expInfo['frameRate']=win.getActualFrameRate()
if expInfo['frameRate']!=None:
    frameDur = 1.0/round(expInfo['frameRate'])
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# Initialize components for Routine "Instructions"
InstructionsClock = core.Clock()
text = visual.TextStim(win=win, ori=0, name='text',
    text='Welcome to the task.\n\nYour job in this experiment is to hit the "1" button for every number that appears other than a "3". \nWhen a "3" appears, you will refrain from hitting the "1" button.\n\nEvery so often, you will be presented with two words. Please select the word from the pair that you prefer by hitting the "1" button for the word on the left, and "2" for the word on the right.\n\nPlease press the space key when you are ready to begin.',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=1.5,
    color='black', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "Data_Loader"
Data_LoaderClock = core.Clock()


# Initialize components for Routine "Focus"
FocusClock = core.Clock()
Initial_fix = visual.TextStim(win=win, ori=0, name='Initial_fix',
    text='+',    font='Arial',
    pos=[0, 0], height=0.3, wrapWidth=None,
    color='black', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()
from random import choice
from numpy import average, std

nonTarget=[1, 2, 4, 5, 6, 7, 8, 9]
target=3
prev_RT=[]
curr_RT=[]
test_var = "practice_schedules.xlsx"


num_stim = visual.TextStim(win=win, ori=0, name='num_stim',
    text='default text',    font='Arial',
    pos=[0, 0], height=0.3, wrapWidth=None,
    color='black', colorSpace='rgb', opacity=1,
    depth=-1.0)

# Initialize components for Routine "Focus"
FocusClock = core.Clock()
Initial_fix = visual.TextStim(win=win, ori=0, name='Initial_fix',
    text='+',    font='Arial',
    pos=[0, 0], height=0.3, wrapWidth=None,
    color='black', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "Probe"
ProbeClock = core.Clock()
import random

probe_stim1 = visual.TextStim(win=win, ori=0, name='probe_stim1',
    text='default text',    font='Arial',
    pos=[-.5, 0], height=0.2, wrapWidth=1.5,
    color='black', colorSpace='rgb', opacity=1,
    depth=-1.0)
probe_stim2 = visual.TextStim(win=win, ori=0, name='probe_stim2',
    text='default text',    font='Arial',
    pos=[.5, 0], height=0.2, wrapWidth=1.5,
    color='black', colorSpace='rgb', opacity=1,
    depth=-3.0)

# Initialize components for Routine "Goodbye"
GoodbyeClock = core.Clock()
goodbye_text = visual.TextStim(win=win, ori=0, name='goodbye_text',
    text='Thank you for your participation!\n\nYou are finished and may now leave.',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=1.5,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

#------Prepare to start Routine "Instructions"-------
t = 0
InstructionsClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
space_key = event.BuilderKeyResponse()  # create an object of type KeyResponse
space_key.status = NOT_STARTED
# keep track of which components have finished
InstructionsComponents = []
InstructionsComponents.append(text)
InstructionsComponents.append(space_key)
for thisComponent in InstructionsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "Instructions"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = InstructionsClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text* updates
    if t >= 0.0 and text.status == NOT_STARTED:
        # keep track of start time/frame for later
        text.tStart = t  # underestimates by a little under one frame
        text.frameNStart = frameN  # exact frame index
        text.setAutoDraw(True)
    
    # *space_key* updates
    if t >= 0.0 and space_key.status == NOT_STARTED:
        # keep track of start time/frame for later
        space_key.tStart = t  # underestimates by a little under one frame
        space_key.frameNStart = frameN  # exact frame index
        space_key.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(space_key.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if space_key.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            space_key.keys = theseKeys[-1]  # just the last key pressed
            space_key.rt = space_key.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in InstructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "Instructions"-------
for thisComponent in InstructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if space_key.keys in ['', [], None]:  # No response was made
   space_key.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('space_key.keys',space_key.keys)
if space_key.keys != None:  # we had a response
    thisExp.addData('space_key.rt', space_key.rt)
thisExp.nextEntry()
# the Routine "Instructions" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
condition_loop = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='condition_loop')
thisExp.addLoop(condition_loop)  # add the loop to the experiment
thisCondition_loop = condition_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb=thisCondition_loop.rgb)
if thisCondition_loop != None:
    for paramName in thisCondition_loop.keys():
        exec(paramName + '= thisCondition_loop.' + paramName)

for thisCondition_loop in condition_loop:
    currentLoop = condition_loop
    # abbreviate parameter names if possible (e.g. rgb = thisCondition_loop.rgb)
    if thisCondition_loop != None:
        for paramName in thisCondition_loop.keys():
            exec(paramName + '= thisCondition_loop.' + paramName)
    
    #------Prepare to start Routine "Data_Loader"-------
    t = 0
    Data_LoaderClock.reset()  # clock 
    frameN = -1
    # update component parameters for each repeat
    from psychopy import data
    data_vals = data.importConditions("practice_Wordlist.xlsx")
    
    
    # keep track of which components have finished
    Data_LoaderComponents = []
    for thisComponent in Data_LoaderComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    #-------Start Routine "Data_Loader"-------
    continueRoutine = True
    while continueRoutine:
        # get current time
        t = Data_LoaderClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in Data_LoaderComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    #-------Ending Routine "Data_Loader"-------
    for thisComponent in Data_LoaderComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    # the Routine "Data_Loader" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=1, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(test_var),
        seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb=thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial.keys():
            exec(paramName + '= thisTrial.' + paramName)
    
    for thisTrial in trials:
        currentLoop = trials
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial.keys():
                exec(paramName + '= thisTrial.' + paramName)
        
        #------Prepare to start Routine "Focus"-------
        t = 0
        FocusClock.reset()  # clock 
        frameN = -1
        routineTimer.add(2.000000)
        # update component parameters for each repeat
        # keep track of which components have finished
        FocusComponents = []
        FocusComponents.append(Initial_fix)
        for thisComponent in FocusComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        #-------Start Routine "Focus"-------
        continueRoutine = True
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = FocusClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *Initial_fix* updates
            if t >= 0.0 and Initial_fix.status == NOT_STARTED:
                # keep track of start time/frame for later
                Initial_fix.tStart = t  # underestimates by a little under one frame
                Initial_fix.frameNStart = frameN  # exact frame index
                Initial_fix.setAutoDraw(True)
            if Initial_fix.status == STARTED and t >= (0.0 + (2.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                Initial_fix.setAutoDraw(False)
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in FocusComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # check for quit (the Esc key)
            if endExpNow or event.getKeys(keyList=["escape"]):
                core.quit()
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        #-------Ending Routine "Focus"-------
        for thisComponent in FocusComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        
        #------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock 
        frameN = -1
        routineTimer.add(2.000000)
        # update component parameters for each repeat
        #In a trial, load current digit. 
        
        i = trials.thisTrialN
        if trial_schedule==2:
            stim = target
        else:
            stim = choice(nonTarget)
        
        num_stim.setText(stim)
        trial_resp = event.BuilderKeyResponse()  # create an object of type KeyResponse
        trial_resp.status = NOT_STARTED
        # keep track of which components have finished
        trialComponents = []
        trialComponents.append(num_stim)
        trialComponents.append(trial_resp)
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        #-------Start Routine "trial"-------
        continueRoutine = True
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            
            # *num_stim* updates
            if t >= 0.0 and num_stim.status == NOT_STARTED:
                # keep track of start time/frame for later
                num_stim.tStart = t  # underestimates by a little under one frame
                num_stim.frameNStart = frameN  # exact frame index
                num_stim.setAutoDraw(True)
            if num_stim.status == STARTED and t >= (0.0 + (2.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                num_stim.setAutoDraw(False)
            
            # *trial_resp* updates
            if t >= 0.0 and trial_resp.status == NOT_STARTED:
                # keep track of start time/frame for later
                trial_resp.tStart = t  # underestimates by a little under one frame
                trial_resp.frameNStart = frameN  # exact frame index
                trial_resp.status = STARTED
                # keyboard checking is just starting
                win.callOnFlip(trial_resp.clock.reset)  # t=0 on next screen flip
                event.clearEvents(eventType='keyboard')
            if trial_resp.status == STARTED and t >= (0.0 + (2-win.monitorFramePeriod*0.75)): #most of one frame period left
                trial_resp.status = STOPPED
            if trial_resp.status == STARTED:
                theseKeys = event.getKeys(keyList=['1'])
                
                # check for quit:
                if "escape" in theseKeys:
                    endExpNow = True
                if len(theseKeys) > 0:  # at least one key was pressed
                    trial_resp.keys = theseKeys[-1]  # just the last key pressed
                    trial_resp.rt = trial_resp.clock.getTime()
                    # a response ends the routine
                    continueRoutine = False
            
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
        
        #-------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        #At the end save whether inhibition trial or not.
        #At the end of each trial, compute current average and SD
        #Set the conditional value that determines if the probe task should initiate
        should_probe_run = bool(probe_schedule)
        if stim==target:
            trials.addData('inhibition_trial', True)
            if type(trial_resp.rt) is not float:
                trial_resp.rt = 0
                trial_resp.keys = None
        
        else:
            trials.addData('inhibition_trial', False)
        
        if stim in nonTarget and trial_resp.keys=='1':
            hit,cr,com,omis=[1,0,0,0]
        elif stim in nonTarget and trial_resp.keys!='1':
            hit,cr,com,omis=[0,0,0,1]
        elif stim==target and trial_resp.keys=='1':
            hit,cr,com,omis=[0,0,1,0]
        elif stim==target and trial_resp.keys!='1':
            hit,cr,com,omis=[0,1,0,0]
        trials.addData('hit', hit)
        trials.addData('cr', cr)
        trials.addData('com', com)
        trials.addData('omis', omis)
        trials.addData('stim', stim)
        # check responses
        if trial_resp.keys in ['', [], None]:  # No response was made
           trial_resp.keys=None
        # store data for trials (TrialHandler)
        trials.addData('trial_resp.keys',trial_resp.keys)
        if trial_resp.keys != None:  # we had a response
            trials.addData('trial_resp.rt', trial_resp.rt)
        
        # set up handler to look after randomisation of conditions etc
        conditional_run = data.TrialHandler(nReps=should_probe_run, method='random', 
            extraInfo=expInfo, originPath=-1,
            trialList=[None],
            seed=None, name='conditional_run')
        thisExp.addLoop(conditional_run)  # add the loop to the experiment
        thisConditional_run = conditional_run.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb=thisConditional_run.rgb)
        if thisConditional_run != None:
            for paramName in thisConditional_run.keys():
                exec(paramName + '= thisConditional_run.' + paramName)
        
        for thisConditional_run in conditional_run:
            currentLoop = conditional_run
            # abbreviate parameter names if possible (e.g. rgb = thisConditional_run.rgb)
            if thisConditional_run != None:
                for paramName in thisConditional_run.keys():
                    exec(paramName + '= thisConditional_run.' + paramName)
            
            #------Prepare to start Routine "Focus"-------
            t = 0
            FocusClock.reset()  # clock 
            frameN = -1
            routineTimer.add(2.000000)
            # update component parameters for each repeat
            # keep track of which components have finished
            FocusComponents = []
            FocusComponents.append(Initial_fix)
            for thisComponent in FocusComponents:
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            
            #-------Start Routine "Focus"-------
            continueRoutine = True
            while continueRoutine and routineTimer.getTime() > 0:
                # get current time
                t = FocusClock.getTime()
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *Initial_fix* updates
                if t >= 0.0 and Initial_fix.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    Initial_fix.tStart = t  # underestimates by a little under one frame
                    Initial_fix.frameNStart = frameN  # exact frame index
                    Initial_fix.setAutoDraw(True)
                if Initial_fix.status == STARTED and t >= (0.0 + (2.0-win.monitorFramePeriod*0.75)): #most of one frame period left
                    Initial_fix.setAutoDraw(False)
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in FocusComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # check for quit (the Esc key)
                if endExpNow or event.getKeys(keyList=["escape"]):
                    core.quit()
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            #-------Ending Routine "Focus"-------
            for thisComponent in FocusComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            
            #------Prepare to start Routine "Probe"-------
            t = 0
            ProbeClock.reset()  # clock 
            frameN = -1
            routineTimer.add(1.700000)
            # update component parameters for each repeat
            #Load word pairs (self and non-self): Word 1 Word 2
            #Randomize side on which each word type is displayed
            #Save which side Self word is displayed
            #Save which words are displayed
            
            i = random.randint(0, len(data_vals)-1)
            non_self = data_vals[i]["non_self"]
            super_self = data_vals[i]["super_self"]
            
            del data_vals[i]
            
            self_side = 1
            if random.random() >= .5:
                word1 = non_self
                word2 = super_self
                self_side = 2
            else:
                word1 = super_self
                word2 = non_self
            
            conditional_run.addData("self_side", self_side)
            conditional_run.addData("super_self", super_self)
            conditional_run.addData("non_self", non_self)
            probe_stim1.setText(word1)
            probe_resp = event.BuilderKeyResponse()  # create an object of type KeyResponse
            probe_resp.status = NOT_STARTED
            probe_stim2.setText(word2)
            # keep track of which components have finished
            ProbeComponents = []
            ProbeComponents.append(probe_stim1)
            ProbeComponents.append(probe_resp)
            ProbeComponents.append(probe_stim2)
            for thisComponent in ProbeComponents:
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            
            #-------Start Routine "Probe"-------
            continueRoutine = True
            while continueRoutine and routineTimer.getTime() > 0:
                # get current time
                t = ProbeClock.getTime()
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                
                # *probe_stim1* updates
                if t >= 0.0 and probe_stim1.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    probe_stim1.tStart = t  # underestimates by a little under one frame
                    probe_stim1.frameNStart = frameN  # exact frame index
                    probe_stim1.setAutoDraw(True)
                if probe_stim1.status == STARTED and t >= (0.0 + (1-win.monitorFramePeriod*0.75)): #most of one frame period left
                    probe_stim1.setAutoDraw(False)
                
                # *probe_resp* updates
                if t >= 0.0 and probe_resp.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    probe_resp.tStart = t  # underestimates by a little under one frame
                    probe_resp.frameNStart = frameN  # exact frame index
                    probe_resp.status = STARTED
                    # keyboard checking is just starting
                    win.callOnFlip(probe_resp.clock.reset)  # t=0 on next screen flip
                    event.clearEvents(eventType='keyboard')
                if probe_resp.status == STARTED and t >= (0.0 + (1.7-win.monitorFramePeriod*0.75)): #most of one frame period left
                    probe_resp.status = STOPPED
                if probe_resp.status == STARTED:
                    theseKeys = event.getKeys(keyList=['1', '2'])
                    
                    # check for quit:
                    if "escape" in theseKeys:
                        endExpNow = True
                    if len(theseKeys) > 0:  # at least one key was pressed
                        probe_resp.keys = theseKeys[-1]  # just the last key pressed
                        probe_resp.rt = probe_resp.clock.getTime()
                        # a response ends the routine
                        continueRoutine = False
                
                # *probe_stim2* updates
                if t >= 0.0 and probe_stim2.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    probe_stim2.tStart = t  # underestimates by a little under one frame
                    probe_stim2.frameNStart = frameN  # exact frame index
                    probe_stim2.setAutoDraw(True)
                if probe_stim2.status == STARTED and t >= (0.0 + (1-win.monitorFramePeriod*0.75)): #most of one frame period left
                    probe_stim2.setAutoDraw(False)
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in ProbeComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # check for quit (the Esc key)
                if endExpNow or event.getKeys(keyList=["escape"]):
                    core.quit()
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            #-------Ending Routine "Probe"-------
            for thisComponent in ProbeComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            
            # check responses
            if probe_resp.keys in ['', [], None]:  # No response was made
               probe_resp.keys=None
            # store data for conditional_run (TrialHandler)
            conditional_run.addData('probe_resp.keys',probe_resp.keys)
            if probe_resp.keys != None:  # we had a response
                conditional_run.addData('probe_resp.rt', probe_resp.rt)
            thisExp.nextEntry()
            
        # completed should_probe_run repeats of 'conditional_run'
        
        # get names of stimulus parameters
        if conditional_run.trialList in ([], [None], None):  params = []
        else:  params = conditional_run.trialList[0].keys()
        # save data for this loop
        conditional_run.saveAsExcel(filename + '.xlsx', sheetName='conditional_run',
            stimOut=params,
            dataOut=['n','all_mean','all_std', 'all_raw'])
        thisExp.nextEntry()
        
    # completed 1 repeats of 'trials'
    
    # get names of stimulus parameters
    if trials.trialList in ([], [None], None):  params = []
    else:  params = trials.trialList[0].keys()
    # save data for this loop
    trials.saveAsExcel(filename + '.xlsx', sheetName='trials',
        stimOut=params,
        dataOut=['n','all_mean','all_std', 'all_raw'])
    thisExp.nextEntry()
    
# completed 1 repeats of 'condition_loop'

# get names of stimulus parameters
if condition_loop.trialList in ([], [None], None):  params = []
else:  params = condition_loop.trialList[0].keys()
# save data for this loop
condition_loop.saveAsExcel(filename + '.xlsx', sheetName='condition_loop',
    stimOut=params,
    dataOut=['n','all_mean','all_std', 'all_raw'])

#------Prepare to start Routine "Goodbye"-------
t = 0
GoodbyeClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
goodbye_key = event.BuilderKeyResponse()  # create an object of type KeyResponse
goodbye_key.status = NOT_STARTED
# keep track of which components have finished
GoodbyeComponents = []
GoodbyeComponents.append(goodbye_text)
GoodbyeComponents.append(goodbye_key)
for thisComponent in GoodbyeComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "Goodbye"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = GoodbyeClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *goodbye_text* updates
    if t >= 0.0 and goodbye_text.status == NOT_STARTED:
        # keep track of start time/frame for later
        goodbye_text.tStart = t  # underestimates by a little under one frame
        goodbye_text.frameNStart = frameN  # exact frame index
        goodbye_text.setAutoDraw(True)
    
    # *goodbye_key* updates
    if t >= 0.0 and goodbye_key.status == NOT_STARTED:
        # keep track of start time/frame for later
        goodbye_key.tStart = t  # underestimates by a little under one frame
        goodbye_key.frameNStart = frameN  # exact frame index
        goodbye_key.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(goodbye_key.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if goodbye_key.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            goodbye_key.keys = theseKeys[-1]  # just the last key pressed
            goodbye_key.rt = goodbye_key.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in GoodbyeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "Goodbye"-------
for thisComponent in GoodbyeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if goodbye_key.keys in ['', [], None]:  # No response was made
   goodbye_key.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('goodbye_key.keys',goodbye_key.keys)
if goodbye_key.keys != None:  # we had a response
    thisExp.addData('goodbye_key.rt', goodbye_key.rt)
thisExp.nextEntry()
# the Routine "Goodbye" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()



win.close()
core.quit()
