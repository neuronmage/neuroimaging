#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.82.01), Wed Sep  2 15:01:38 2015
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

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'WordRatingTask'  # from the Builder filename that created this script
expInfo = {u'session': u'001', u'participant': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + 'data/%s_%s_%s' %(expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=False,
    dataFileName=filename)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(size=(1280, 800), fullscr=True, screen=0, allowGUI=True, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
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
instructions1 = visual.TextStim(win=win, ori=0, name='instructions1',
    text=u'Some words have the ability to remind us about ourselves and our day-to-day life, and these are known as \u201cself-referencing\u201d words. \n\nIn this experiment, you will be asked to rate 800 individual words on a scale ranging from 1 = \u201cNot At All\u201d to 5 = \u201cVery\u201d self-referencing. Use the mouse to select your answer on the sliding response scale then confirm your response by clicking the button below.\n\nPlease hit the space bar when you are ready to begin.\n',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=1.5,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()
word_present = visual.TextStim(win=win, ori=0, name='word_present',
    text='default text',    font='Arial',
    pos=[0, 0], height=.2, wrapWidth=None,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)
rating = visual.RatingScale(win=win, name='rating', marker='triangle', size=1.5, pos=[0.0, -0.4], choices=[u'Not At All', u'2', u'3', u'4', u'Very'], tickHeight=-1)

# Initialize components for Routine "Goodbye"
GoodbyeClock = core.Clock()
goodbye1 = visual.TextStim(win=win, ori=0, name='goodbye1',
    text='Thank you for participating in our experiment.\n\nIf you are completing this with another participant, please sit quietly until they are finished.',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
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
key_space = event.BuilderKeyResponse()  # create an object of type KeyResponse
key_space.status = NOT_STARTED
# keep track of which components have finished
InstructionsComponents = []
InstructionsComponents.append(instructions1)
InstructionsComponents.append(key_space)
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
    
    # *instructions1* updates
    if t >= 0.0 and instructions1.status == NOT_STARTED:
        # keep track of start time/frame for later
        instructions1.tStart = t  # underestimates by a little under one frame
        instructions1.frameNStart = frameN  # exact frame index
        instructions1.setAutoDraw(True)
    
    # *key_space* updates
    if t >= 0.0 and key_space.status == NOT_STARTED:
        # keep track of start time/frame for later
        key_space.tStart = t  # underestimates by a little under one frame
        key_space.frameNStart = frameN  # exact frame index
        key_space.status = STARTED
        # keyboard checking is just starting
        key_space.clock.reset()  # now t=0
        event.clearEvents(eventType='keyboard')
    if key_space.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            key_space.keys = theseKeys[-1]  # just the last key pressed
            key_space.rt = key_space.clock.getTime()
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
if key_space.keys in ['', [], None]:  # No response was made
   key_space.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('key_space.keys',key_space.keys)
if key_space.keys != None:  # we had a response
    thisExp.addData('key_space.rt', key_space.rt)
thisExp.nextEntry()
# the Routine "Instructions" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
word_loop = data.TrialHandler(nReps=1, method='random', 
    extraInfo=expInfo, originPath=None,
    trialList=data.importConditions('Wordlistexcel.xlsx'),
    seed=None, name='word_loop')
thisExp.addLoop(word_loop)  # add the loop to the experiment
thisWord_loop = word_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb=thisWord_loop.rgb)
if thisWord_loop != None:
    for paramName in thisWord_loop.keys():
        exec(paramName + '= thisWord_loop.' + paramName)

for thisWord_loop in word_loop:
    currentLoop = word_loop
    # abbreviate parameter names if possible (e.g. rgb = thisWord_loop.rgb)
    if thisWord_loop != None:
        for paramName in thisWord_loop.keys():
            exec(paramName + '= thisWord_loop.' + paramName)
    
    #------Prepare to start Routine "trial"-------
    t = 0
    trialClock.reset()  # clock 
    frameN = -1
    # update component parameters for each repeat
    word_present.setText(Word)
    rating.reset()
    # keep track of which components have finished
    trialComponents = []
    trialComponents.append(word_present)
    trialComponents.append(rating)
    for thisComponent in trialComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    #-------Start Routine "trial"-------
    continueRoutine = True
    while continueRoutine:
        # get current time
        t = trialClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *word_present* updates
        if t >= 0.0 and word_present.status == NOT_STARTED:
            # keep track of start time/frame for later
            word_present.tStart = t  # underestimates by a little under one frame
            word_present.frameNStart = frameN  # exact frame index
            word_present.setAutoDraw(True)
        # *rating* updates
        if t >= 0.0 and rating.status == NOT_STARTED:
            # keep track of start time/frame for later
            rating.tStart = t  # underestimates by a little under one frame
            rating.frameNStart = frameN  # exact frame index
            rating.setAutoDraw(True)
        continueRoutine &= rating.noResponse  # a response ends the trial
        
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
    # store data for word_loop (TrialHandler)
    word_loop.addData('rating.response', rating.getRating())
    word_loop.addData('rating.rt', rating.getRT())
    # the Routine "trial" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 1 repeats of 'word_loop'

# get names of stimulus parameters
if word_loop.trialList in ([], [None], None):  params = []
else:  params = word_loop.trialList[0].keys()
# save data for this loop
word_loop.saveAsExcel(filename + '.xlsx', sheetName='word_loop',
    stimOut=params,
    dataOut=['n','all_mean','all_std', 'all_raw'])

#------Prepare to start Routine "Goodbye"-------
t = 0
GoodbyeClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
exit_key = event.BuilderKeyResponse()  # create an object of type KeyResponse
exit_key.status = NOT_STARTED
# keep track of which components have finished
GoodbyeComponents = []
GoodbyeComponents.append(goodbye1)
GoodbyeComponents.append(exit_key)
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
    
    # *goodbye1* updates
    if t >= 0.0 and goodbye1.status == NOT_STARTED:
        # keep track of start time/frame for later
        goodbye1.tStart = t  # underestimates by a little under one frame
        goodbye1.frameNStart = frameN  # exact frame index
        goodbye1.setAutoDraw(True)
    
    # *exit_key* updates
    if t >= 0.0 and exit_key.status == NOT_STARTED:
        # keep track of start time/frame for later
        exit_key.tStart = t  # underestimates by a little under one frame
        exit_key.frameNStart = frameN  # exact frame index
        exit_key.status = STARTED
        # keyboard checking is just starting
        exit_key.clock.reset()  # now t=0
        event.clearEvents(eventType='keyboard')
    if exit_key.status == STARTED:
        theseKeys = event.getKeys(keyList=['return'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            exit_key.keys = theseKeys[-1]  # just the last key pressed
            exit_key.rt = exit_key.clock.getTime()
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
if exit_key.keys in ['', [], None]:  # No response was made
   exit_key.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('exit_key.keys',exit_key.keys)
if exit_key.keys != None:  # we had a response
    thisExp.addData('exit_key.rt', exit_key.rt)
thisExp.nextEntry()
# the Routine "Goodbye" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()
win.close()
core.quit()
