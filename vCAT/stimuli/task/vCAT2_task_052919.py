#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2020.2.4),
    on November 20, 2020, at 14:04
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2020.2.4'
expName = 'vCAT_task'  # from the Builder filename that created this script
expInfo = {'participant': ''}
dlg = gui.DlgFromDict(dictionary=expInfo, sort_keys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='G:\\My Drive\\MaDLab\\experiments\\vCAT2\\code\\vCAT2_task_052919.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=[1280, 800], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-0.5,-0.5,-0.5], colorSpace='rgb',
    blendMode='avg', useFBO=True)
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "instructions"
instructionsClock = core.Clock()
bkg_instructions = visual.ImageStim(
    win=win,
    name='bkg_instructions', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
text_instructions = visual.TextStim(win=win, name='text_instructions',
    text='The experiment will begin soon. \nPlease remain still and wait for the scanner to start.\n',
    font='Ariel',
    pos=(0, 0), height=0.08, wrapWidth=1.5, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);
start_key = keyboard.Keyboard()
expClock = core.Clock() #begins clock

import os
import numpy as np
from numpy.random import random, randint, normal, shuffle
import random
from math import ceil
from glob import glob

def fixed(v): #function to isolate fixed trials
    return v >= 0 and v < 4

def run_trialorder(): #function to create pseudo-random trial schedule for each run
    dummy_types = [0,1,2,3]*25 + [4]*50 + [-1]*30  #-1 = baseline, 0-3 = fixed, 4 = conditional
    rand_dummies = np.random.permutation(dummy_types)
    i = 0
    while i < len(rand_dummies[:-1]):
        v = rand_dummies[i]
        v2 = rand_dummies[i+1]
        if (v2 >= 4 and v >= 4) or (i == 0 and (v >=4 or v == -1)) or (v == -1 and (v2 == -1 or v2 >= 4)) or (v >= 4 and v2 == -1):
            for j, safe_val in enumerate(rand_dummies):
                v0 = None
                safe_v2 = None
                if j + 1 != len(rand_dummies):
                    safe_v2 = rand_dummies[j + 1]
                    if fixed(safe_v2) and fixed(safe_val):
                        temp = rand_dummies[j]
                        rand_dummies[j] = rand_dummies[i]
                        rand_dummies[i] = temp
                        i = -1
                        break
                if j > 0: #if not the first value
                    v0 = rand_dummies[j - 1]
                    if fixed(v0) and fixed(safe_val):
                        temp = rand_dummies[j]
                        rand_dummies[j] = rand_dummies[i]
                        rand_dummies[i] = temp
                        i = -1
                        break
        i += 1
    return rand_dummies

#function to place object images on sides
def sides_assign(corr_img, incorr_img):
    rand_side = np.random.permutation([1,2]) #randomize correct response button at the first index
    curr_corr = rand_side[0] #correct response based on side 
    if curr_corr == 1: #if correct side is left
        left_box_img = corr_img
        right_box_img = incorr_img
    else: #if correct side is right (corr_resp == 2)
        left_box_img = incorr_img
        right_box_img = corr_img
    return left_box_img, right_box_img, curr_corr

##Set Dependencies
#stims (x5)
all_kaleidos = glob(os.getcwd() + '/kaleido_img/*jpg') #import all images to be used
all_kaleidos = np.random.permutation(all_kaleidos) #randomize order to grab different images each time
#practice_kaleidos = all_kaleidos[:15] #5 for practice run
task_kaleidos = all_kaleidos[15:] #all remaining for actual task runs
#object (x4)
all_objects = glob(os.getcwd() + '/object_img/*jpg')
all_objects = np.random.permutation(all_objects)
#practice_objects = all_objects[:12] # 4 object images for practice run
task_objects = all_objects[12:]
#scene (x1)
all_scenes = glob(os.getcwd() + '/scene_img/*jpg')
all_scenes = np.random.permutation(all_scenes) 
#practice_scene = all_scenes[0] # 1 scene for practice choice
task_scenes = all_scenes[3:]
#face (x1)
all_faces = glob(os.getcwd() + '/face_img/*jpg')
all_faces = np.random.permutation(all_faces)
#practice_face = all_faces[0] # 1 face for practice choice
task_faces = all_faces[3:]

trialbox_color = 'gray'
color_1 = color_2 = trialbox_color
stim_opacity = 0

# Initialize components for Routine "stim_prep"
stim_prepClock = core.Clock()
##Set Dependencies
#stims (x5)
all_kaleidos = glob(os.getcwd() + '/kaleido_img/*jpg') #import all images to be used
all_kaleidos = np.random.permutation(all_kaleidos) #randomize order to grab different images each time
task_kaleidos = all_kaleidos[15:] #all remaining for actual task runs
#object (x4)
all_objects = glob(os.getcwd() + '/object_img/*jpg')
all_objects = np.random.permutation(all_objects)
task_objects = all_objects[12:]
#scene (x1)
all_scenes = glob(os.getcwd() + '/scene_img/*jpg')
all_scenes = np.random.permutation(all_scenes) 
task_scenes = all_scenes[3:]
#face (x1)
all_faces = glob(os.getcwd() + '/face_img/*/*bmp')
all_faces = np.random.permutation(all_faces)
task_faces = all_faces[3:]

all_kaleido_sets = []
all_object_sets = []
all_scene_sets = []
all_face_sets = []
for set in range(2):
    kaleido_stimset = task_kaleidos[set:30:5]
    object_stimset = task_objects[set:24:4]
    scene_stimset = task_scenes[set:60:3]
    face_stimset = task_faces[set:60:3]
    all_kaleido_sets.append(kaleido_stimset)
    all_object_sets.append(object_stimset)
    all_scene_sets.append(scene_stimset)
    all_face_sets.append(face_stimset)

# Initialize components for Routine "delay_init"
delay_initClock = core.Clock()
image_delay = visual.ImageStim(
    win=win,
    name='image_delay', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fixation = visual.TextStim(win=win, name='fixation',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "trial"
trialClock = core.Clock()
trialbox_color = 'gray'
color_1 = color_2 = trialbox_color
stim_opacity = 0
object_opacity = 0
isi_bkg_1 = visual.ImageStim(
    win=win,
    name='isi_bkg_1', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-1.0)
bkg_1 = visual.ImageStim(
    win=win,
    name='bkg_1', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color='black', colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
fixation_1 = visual.TextStim(win=win, name='fixation_1',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-3.0);
stim_img_1 = visual.ImageStim(
    win=win,
    name='stim_img_1', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(0.50, 0.75),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-4.0)
left_resp_1 = visual.Rect(
    win=win, name='left_resp_1',
    width=(0.55, 0.65)[0], height=(0.55, 0.65)[1],
    ori=0, pos=(-0.65, 0),
    lineWidth=1, lineColor='black', lineColorSpace='rgb',
    fillColor=1.0, fillColorSpace='rgb',
    opacity=1.0, depth=-5.0, interpolate=True)
right_resp_1 = visual.Rect(
    win=win, name='right_resp_1',
    width=(0.55, 0.65)[0], height=(0.55, 0.65)[1],
    ori=0, pos=(0.65, 0),
    lineWidth=1, lineColor='black', lineColorSpace='rgb',
    fillColor=1.0, fillColorSpace='rgb',
    opacity=1.0, depth=-6.0, interpolate=True)
left_obj = visual.ImageStim(
    win=win,
    name='left_obj', 
    image='sin', mask=None,
    ori=0, pos=(-0.65, 0), size=(0.5, 0.6),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-7.0)
right_obj = visual.ImageStim(
    win=win,
    name='right_obj', 
    image='sin', mask=None,
    ori=0, pos=(0.65, 0), size=(0.5, 0.6),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-8.0)
resp_trial_1 = keyboard.Keyboard()
cue_1 = visual.TextStim(win=win, name='cue_1',
    text='Go!',
    font='Arial',
    pos=(0, 0), height=0.4, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-10.0);
fixation_2 = visual.TextStim(win=win, name='fixation_2',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-11.0);

# Initialize components for Routine "feedback"
feedbackClock = core.Clock()
fb_bkg = visual.ImageStim(
    win=win,
    name='fb_bkg', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fb_text = visual.TextStim(win=win, name='fb_text',
    text='default text',
    font='Ariel',
    pos=(0, 0), height=0.5, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=0.8, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "delay_end"
delay_endClock = core.Clock()
image = visual.ImageStim(
    win=win,
    name='image', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fix_end = visual.TextStim(win=win, name='fix_end',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "run_break"
run_breakClock = core.Clock()
text_run_break = visual.TextStim(win=win, name='text_run_break',
    text='You have completed half of the first image set. \n\nPlease wait for the scanner to begin to complete the second half.',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);
key_run_break = keyboard.Keyboard()

# Initialize components for Routine "delay_init"
delay_initClock = core.Clock()
image_delay = visual.ImageStim(
    win=win,
    name='image_delay', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fixation = visual.TextStim(win=win, name='fixation',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "trial"
trialClock = core.Clock()
trialbox_color = 'gray'
color_1 = color_2 = trialbox_color
stim_opacity = 0
object_opacity = 0
isi_bkg_1 = visual.ImageStim(
    win=win,
    name='isi_bkg_1', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-1.0)
bkg_1 = visual.ImageStim(
    win=win,
    name='bkg_1', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color='black', colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
fixation_1 = visual.TextStim(win=win, name='fixation_1',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-3.0);
stim_img_1 = visual.ImageStim(
    win=win,
    name='stim_img_1', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(0.50, 0.75),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-4.0)
left_resp_1 = visual.Rect(
    win=win, name='left_resp_1',
    width=(0.55, 0.65)[0], height=(0.55, 0.65)[1],
    ori=0, pos=(-0.65, 0),
    lineWidth=1, lineColor='black', lineColorSpace='rgb',
    fillColor=1.0, fillColorSpace='rgb',
    opacity=1.0, depth=-5.0, interpolate=True)
right_resp_1 = visual.Rect(
    win=win, name='right_resp_1',
    width=(0.55, 0.65)[0], height=(0.55, 0.65)[1],
    ori=0, pos=(0.65, 0),
    lineWidth=1, lineColor='black', lineColorSpace='rgb',
    fillColor=1.0, fillColorSpace='rgb',
    opacity=1.0, depth=-6.0, interpolate=True)
left_obj = visual.ImageStim(
    win=win,
    name='left_obj', 
    image='sin', mask=None,
    ori=0, pos=(-0.65, 0), size=(0.5, 0.6),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-7.0)
right_obj = visual.ImageStim(
    win=win,
    name='right_obj', 
    image='sin', mask=None,
    ori=0, pos=(0.65, 0), size=(0.5, 0.6),
    color=[1,1,1], colorSpace='rgb', opacity=1.0,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-8.0)
resp_trial_1 = keyboard.Keyboard()
cue_1 = visual.TextStim(win=win, name='cue_1',
    text='Go!',
    font='Arial',
    pos=(0, 0), height=0.4, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-10.0);
fixation_2 = visual.TextStim(win=win, name='fixation_2',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-11.0);

# Initialize components for Routine "feedback"
feedbackClock = core.Clock()
fb_bkg = visual.ImageStim(
    win=win,
    name='fb_bkg', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fb_text = visual.TextStim(win=win, name='fb_text',
    text='default text',
    font='Ariel',
    pos=(0, 0), height=0.5, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=0.8, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "delay_end"
delay_endClock = core.Clock()
image = visual.ImageStim(
    win=win,
    name='image', 
    image=None, mask=None,
    ori=0, pos=(0, 0), size=(2,2),
    color=[-0.5,-0.5,-0.5], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)
fix_end = visual.TextStim(win=win, name='fix_end',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.3, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "set_break"
set_breakClock = core.Clock()
loop_text2 = None
text_set_break = visual.TextStim(win=win, name='text_set_break',
    text='default text',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);
set_break_key = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "instructions"-------
continueRoutine = True
# update component parameters for each repeat
start_key.keys = []
start_key.rt = []
_start_key_allKeys = []
# keep track of which components have finished
instructionsComponents = [bkg_instructions, text_instructions, start_key]
for thisComponent in instructionsComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
instructionsClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "instructions"-------
while continueRoutine:
    # get current time
    t = instructionsClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=instructionsClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *bkg_instructions* updates
    if bkg_instructions.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        bkg_instructions.frameNStart = frameN  # exact frame index
        bkg_instructions.tStart = t  # local t and not account for scr refresh
        bkg_instructions.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(bkg_instructions, 'tStartRefresh')  # time at next scr refresh
        bkg_instructions.setAutoDraw(True)
    
    # *text_instructions* updates
    if text_instructions.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_instructions.frameNStart = frameN  # exact frame index
        text_instructions.tStart = t  # local t and not account for scr refresh
        text_instructions.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_instructions, 'tStartRefresh')  # time at next scr refresh
        text_instructions.setAutoDraw(True)
    
    # *start_key* updates
    waitOnFlip = False
    if start_key.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        start_key.frameNStart = frameN  # exact frame index
        start_key.tStart = t  # local t and not account for scr refresh
        start_key.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(start_key, 'tStartRefresh')  # time at next scr refresh
        start_key.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(start_key.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(start_key.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if start_key.status == STARTED and not waitOnFlip:
        theseKeys = start_key.getKeys(keyList=['space', '5'], waitRelease=False)
        _start_key_allKeys.extend(theseKeys)
        if len(_start_key_allKeys):
            start_key.keys = _start_key_allKeys[-1].name  # just the last key pressed
            start_key.rt = _start_key_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "instructions"-------
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('bkg_instructions.started', bkg_instructions.tStartRefresh)
thisExp.addData('bkg_instructions.stopped', bkg_instructions.tStopRefresh)
thisExp.addData('text_instructions.started', text_instructions.tStartRefresh)
thisExp.addData('text_instructions.stopped', text_instructions.tStopRefresh)
# check responses
if start_key.keys in ['', [], None]:  # No response was made
    start_key.keys = None
thisExp.addData('start_key.keys',start_key.keys)
if start_key.keys != None:  # we had a response
    thisExp.addData('start_key.rt', start_key.rt)
thisExp.addData('start_key.started', start_key.tStartRefresh)
thisExp.addData('start_key.stopped', start_key.tStopRefresh)
thisExp.nextEntry()
expClock.reset() #rests clock at the end of instructions
# the Routine "instructions" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
set_loop = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('loop_schedule2.xlsx'),
    seed=None, name='set_loop')
thisExp.addLoop(set_loop)  # add the loop to the experiment
thisSet_loop = set_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisSet_loop.rgb)
if thisSet_loop != None:
    for paramName in thisSet_loop:
        exec('{} = thisSet_loop[paramName]'.format(paramName))

for thisSet_loop in set_loop:
    currentLoop = set_loop
    # abbreviate parameter names if possible (e.g. rgb = thisSet_loop.rgb)
    if thisSet_loop != None:
        for paramName in thisSet_loop:
            exec('{} = thisSet_loop[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "stim_prep"-------
    continueRoutine = True
    # update component parameters for each repeat
    kaleido_stims = all_kaleido_sets[loop_num] #kaleidoscopic images to be used in current set
    object_stims = all_object_sets[loop_num] #object images to be used in current set
    scene_stim = all_scene_sets[loop_num] #scene images to be used in current set
    face_stim = all_face_sets[loop_num] #face images to be used in current set
    
    current_run = 'run_1'
    # checks the trial order to make sure its balanced
    success = False
    while success == False:
        faces = 0
        run_schedule = run_trialorder()
        for i, item in enumerate(run_schedule):
            if i != 179:
                if item in [0,2] and run_schedule[i+1] == -1:
                    faces += 1
        if faces >= 12 and faces <= 18:
            success = True
    
    
    # keep track of which components have finished
    stim_prepComponents = []
    for thisComponent in stim_prepComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    stim_prepClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "stim_prep"-------
    while continueRoutine:
        # get current time
        t = stim_prepClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=stim_prepClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stim_prepComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "stim_prep"-------
    for thisComponent in stim_prepComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    expClock.reset()
    # the Routine "stim_prep" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "delay_init"-------
    continueRoutine = True
    routineTimer.add(5.500000)
    # update component parameters for each repeat
    # keep track of which components have finished
    delay_initComponents = [image_delay, fixation]
    for thisComponent in delay_initComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    delay_initClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "delay_init"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = delay_initClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=delay_initClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *image_delay* updates
        if image_delay.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image_delay.frameNStart = frameN  # exact frame index
            image_delay.tStart = t  # local t and not account for scr refresh
            image_delay.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image_delay, 'tStartRefresh')  # time at next scr refresh
            image_delay.setAutoDraw(True)
        if image_delay.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image_delay.tStartRefresh + 5.5-frameTolerance:
                # keep track of stop time/frame for later
                image_delay.tStop = t  # not accounting for scr refresh
                image_delay.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image_delay, 'tStopRefresh')  # time at next scr refresh
                image_delay.setAutoDraw(False)
        
        # *fixation* updates
        if fixation.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fixation.frameNStart = frameN  # exact frame index
            fixation.tStart = t  # local t and not account for scr refresh
            fixation.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fixation, 'tStartRefresh')  # time at next scr refresh
            fixation.setAutoDraw(True)
        if fixation.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fixation.tStartRefresh + 5.5-frameTolerance:
                # keep track of stop time/frame for later
                fixation.tStop = t  # not accounting for scr refresh
                fixation.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fixation, 'tStopRefresh')  # time at next scr refresh
                fixation.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in delay_initComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "delay_init"-------
    for thisComponent in delay_initComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    set_loop.addData('image_delay.started', image_delay.tStartRefresh)
    set_loop.addData('image_delay.stopped', image_delay.tStopRefresh)
    set_loop.addData('fixation.started', fixation.tStartRefresh)
    set_loop.addData('fixation.stopped', fixation.tStopRefresh)
    
    # set up handler to look after randomisation of conditions etc
    run_1 = data.TrialHandler(nReps=180, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='run_1')
    thisExp.addLoop(run_1)  # add the loop to the experiment
    thisRun_1 = run_1.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRun_1.rgb)
    if thisRun_1 != None:
        for paramName in thisRun_1:
            exec('{} = thisRun_1[paramName]'.format(paramName))
    
    for thisRun_1 in run_1:
        currentLoop = run_1
        # abbreviate parameter names if possible (e.g. rgb = thisRun_1.rgb)
        if thisRun_1 != None:
            for paramName in thisRun_1:
                exec('{} = thisRun_1[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        # update component parameters for each repeat
        runs = eval(current_run)
        runs.addData('StimOnset', expClock.getTime()+0.8)#add length of delay
        exec_run = False
        
        i = runs.thisN #counts the number of routine interation
        if run_schedule[i] == -1: #baseline trial
            trial_background = (np.random.rand(1200,800) > 0.85) * 1 #random static background image
            stim = kaleido_stims[0] #preset opacity of first stim to transparent (0)
            stim_opacity = 0
            object_opacity = 0
            trial_dur = 5.0
            curr_corr = np.random.randint(1,3) #pick random resp box as target
            color_1 = color_2 = 'white' #resp box color
            opacity_1 = opacity_2 = 0.55 #preset both to foil opacity    trial_dur = 8.0
            shape_onset = np.random.uniform(0.8, 3.3)
            shape_dur = 1.7
            resp_onset = shape_onset + 1.0
            resp_dur = 0.7
            fix_dur = shape_onset - 0.1
            fix2_onset = resp_onset + 0.7
            fix2_dur = 5.5 - fix2_onset
            exec('opacity_{0} = 0.65'.format(curr_corr)) #set increased opacity for target
            runs.addData('kaleido_img', 'NA') #if BL trial, no stim image
        else: #task trial
            trial_background = np.ones((1200,800)) * 0.5 #dark gray background
            stim_opacity = 1
            object_opacity = 1
            fix_dur = 1.3
            trial_dur = 2.0
            shape_onset = 0.8
            shape_dur = 1.7 #offset time = 2.5
            resp_onset = 1.8 #0.5 seconds 
            resp_dur = 0.7
            fix2_onset = 0
            fix2_dur = 0
            color_1 = color_2 = 'gray' #inside box is black
            opacity_1 = opacity_2 = 1 #response boxes to solid
            corr_img = None
            incorr_img = None
            #fixed trials
            if run_schedule[i] == 0: #S1F1 - face
                stim = kaleido_stims[0] #set fractal
                corr_img = object_stims[0]
                incorr_img = object_stims[1]
                stim_label = 'face1'
            elif run_schedule[i] == 1: #S1F2 - scene
                stim = kaleido_stims[1]
                corr_img = object_stims[1]
                incorr_img = object_stims[0]
                stim_label = 'scene1'
            elif run_schedule[i] == 2: #S2F1 - face
                stim = kaleido_stims[2]
                corr_img = object_stims[2]
                incorr_img = object_stims[3]
                stim_label = 'face2'
            elif run_schedule[i] == 3: #S2F2 - scene
                stim = kaleido_stims[3]
                corr_img = object_stims[3]
                incorr_img = object_stims[2]
                stim_label = 'scene2'
            #conditional trials
            else:
                stim = kaleido_stims[4]
                if run_schedule[i-1] in (0,2): #previous trial was S1F1 or S2F1 (face)
                    corr_img = random.choice(face_stim)
                    incorr_img = random.choice(scene_stim)
                elif run_schedule[i-1] in (1,3): #previous trial was S1F2 or S2F2 (scene)
                    corr_img = random.choice(scene_stim)
                    incorr_img = random.choice(face_stim)
                stim_label = 'COND'
            left_box_img, right_box_img, curr_corr = sides_assign(corr_img, incorr_img)
            runs.addData('kaleido_img', kaleido_stims[run_schedule[i]][-8:])
            runs.addData('stim_ID', stim_label)
        runs.addData('correct', curr_corr)
        bkg_1.setImage(trial_background)
        stim_img_1.setOpacity(stim_opacity)
        stim_img_1.setImage(stim)
        left_resp_1.setOpacity(opacity_1)
        right_resp_1.setOpacity(opacity_2)
        left_obj.setOpacity(object_opacity)
        left_obj.setImage(left_box_img)
        right_obj.setOpacity(object_opacity)
        right_obj.setImage(right_box_img)
        resp_trial_1.keys = []
        resp_trial_1.rt = []
        _resp_trial_1_allKeys = []
        # keep track of which components have finished
        trialComponents = [isi_bkg_1, bkg_1, fixation_1, stim_img_1, left_resp_1, right_resp_1, left_obj, right_obj, resp_trial_1, cue_1, fixation_2]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            if len(resp_trial_1.keys) != 0 and exec_run == False:
                which_box = resp_trial_1.keys[0]
                exec_run = True
                exec('color_{0} = "yellow"'.format(which_box))
            
            # *isi_bkg_1* updates
            if isi_bkg_1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                isi_bkg_1.frameNStart = frameN  # exact frame index
                isi_bkg_1.tStart = t  # local t and not account for scr refresh
                isi_bkg_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(isi_bkg_1, 'tStartRefresh')  # time at next scr refresh
                isi_bkg_1.setAutoDraw(True)
            if isi_bkg_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > isi_bkg_1.tStartRefresh + 0.6-frameTolerance:
                    # keep track of stop time/frame for later
                    isi_bkg_1.tStop = t  # not accounting for scr refresh
                    isi_bkg_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(isi_bkg_1, 'tStopRefresh')  # time at next scr refresh
                    isi_bkg_1.setAutoDraw(False)
            
            # *bkg_1* updates
            if bkg_1.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
                # keep track of start time/frame for later
                bkg_1.frameNStart = frameN  # exact frame index
                bkg_1.tStart = t  # local t and not account for scr refresh
                bkg_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(bkg_1, 'tStartRefresh')  # time at next scr refresh
                bkg_1.setAutoDraw(True)
            if bkg_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > bkg_1.tStartRefresh + trial_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    bkg_1.tStop = t  # not accounting for scr refresh
                    bkg_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(bkg_1, 'tStopRefresh')  # time at next scr refresh
                    bkg_1.setAutoDraw(False)
            
            # *fixation_1* updates
            if fixation_1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixation_1.frameNStart = frameN  # exact frame index
                fixation_1.tStart = t  # local t and not account for scr refresh
                fixation_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixation_1, 'tStartRefresh')  # time at next scr refresh
                fixation_1.setAutoDraw(True)
            if fixation_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixation_1.tStartRefresh + fix_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    fixation_1.tStop = t  # not accounting for scr refresh
                    fixation_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fixation_1, 'tStopRefresh')  # time at next scr refresh
                    fixation_1.setAutoDraw(False)
            
            # *stim_img_1* updates
            if stim_img_1.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                stim_img_1.frameNStart = frameN  # exact frame index
                stim_img_1.tStart = t  # local t and not account for scr refresh
                stim_img_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(stim_img_1, 'tStartRefresh')  # time at next scr refresh
                stim_img_1.setAutoDraw(True)
            if stim_img_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > stim_img_1.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    stim_img_1.tStop = t  # not accounting for scr refresh
                    stim_img_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(stim_img_1, 'tStopRefresh')  # time at next scr refresh
                    stim_img_1.setAutoDraw(False)
            
            # *left_resp_1* updates
            if left_resp_1.status == NOT_STARTED and tThisFlip >= shape_onset-frameTolerance:
                # keep track of start time/frame for later
                left_resp_1.frameNStart = frameN  # exact frame index
                left_resp_1.tStart = t  # local t and not account for scr refresh
                left_resp_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(left_resp_1, 'tStartRefresh')  # time at next scr refresh
                left_resp_1.setAutoDraw(True)
            if left_resp_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > left_resp_1.tStartRefresh + shape_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    left_resp_1.tStop = t  # not accounting for scr refresh
                    left_resp_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(left_resp_1, 'tStopRefresh')  # time at next scr refresh
                    left_resp_1.setAutoDraw(False)
            if left_resp_1.status == STARTED:  # only update if drawing
                left_resp_1.setFillColor(color_1, log=False)
            
            # *right_resp_1* updates
            if right_resp_1.status == NOT_STARTED and tThisFlip >= shape_onset-frameTolerance:
                # keep track of start time/frame for later
                right_resp_1.frameNStart = frameN  # exact frame index
                right_resp_1.tStart = t  # local t and not account for scr refresh
                right_resp_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(right_resp_1, 'tStartRefresh')  # time at next scr refresh
                right_resp_1.setAutoDraw(True)
            if right_resp_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > right_resp_1.tStartRefresh + shape_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    right_resp_1.tStop = t  # not accounting for scr refresh
                    right_resp_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(right_resp_1, 'tStopRefresh')  # time at next scr refresh
                    right_resp_1.setAutoDraw(False)
            if right_resp_1.status == STARTED:  # only update if drawing
                right_resp_1.setFillColor(color_2, log=False)
            
            # *left_obj* updates
            if left_obj.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                left_obj.frameNStart = frameN  # exact frame index
                left_obj.tStart = t  # local t and not account for scr refresh
                left_obj.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(left_obj, 'tStartRefresh')  # time at next scr refresh
                left_obj.setAutoDraw(True)
            if left_obj.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > left_obj.tStartRefresh + 1.7-frameTolerance:
                    # keep track of stop time/frame for later
                    left_obj.tStop = t  # not accounting for scr refresh
                    left_obj.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(left_obj, 'tStopRefresh')  # time at next scr refresh
                    left_obj.setAutoDraw(False)
            
            # *right_obj* updates
            if right_obj.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                right_obj.frameNStart = frameN  # exact frame index
                right_obj.tStart = t  # local t and not account for scr refresh
                right_obj.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(right_obj, 'tStartRefresh')  # time at next scr refresh
                right_obj.setAutoDraw(True)
            if right_obj.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > right_obj.tStartRefresh + 1.7-frameTolerance:
                    # keep track of stop time/frame for later
                    right_obj.tStop = t  # not accounting for scr refresh
                    right_obj.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(right_obj, 'tStopRefresh')  # time at next scr refresh
                    right_obj.setAutoDraw(False)
            
            # *resp_trial_1* updates
            waitOnFlip = False
            if resp_trial_1.status == NOT_STARTED and tThisFlip >= resp_onset-frameTolerance:
                # keep track of start time/frame for later
                resp_trial_1.frameNStart = frameN  # exact frame index
                resp_trial_1.tStart = t  # local t and not account for scr refresh
                resp_trial_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(resp_trial_1, 'tStartRefresh')  # time at next scr refresh
                resp_trial_1.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(resp_trial_1.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(resp_trial_1.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if resp_trial_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > resp_trial_1.tStartRefresh + resp_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    resp_trial_1.tStop = t  # not accounting for scr refresh
                    resp_trial_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(resp_trial_1, 'tStopRefresh')  # time at next scr refresh
                    resp_trial_1.status = FINISHED
            if resp_trial_1.status == STARTED and not waitOnFlip:
                theseKeys = resp_trial_1.getKeys(keyList=['1', '2'], waitRelease=False)
                _resp_trial_1_allKeys.extend(theseKeys)
                if len(_resp_trial_1_allKeys):
                    resp_trial_1.keys = _resp_trial_1_allKeys[0].name  # just the first key pressed
                    resp_trial_1.rt = _resp_trial_1_allKeys[0].rt
            
            # *cue_1* updates
            if cue_1.status == NOT_STARTED and tThisFlip >= resp_onset-frameTolerance:
                # keep track of start time/frame for later
                cue_1.frameNStart = frameN  # exact frame index
                cue_1.tStart = t  # local t and not account for scr refresh
                cue_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(cue_1, 'tStartRefresh')  # time at next scr refresh
                cue_1.setAutoDraw(True)
            if cue_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > cue_1.tStartRefresh + resp_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    cue_1.tStop = t  # not accounting for scr refresh
                    cue_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(cue_1, 'tStopRefresh')  # time at next scr refresh
                    cue_1.setAutoDraw(False)
            
            # *fixation_2* updates
            if fixation_2.status == NOT_STARTED and tThisFlip >= fix2_onset-frameTolerance:
                # keep track of start time/frame for later
                fixation_2.frameNStart = frameN  # exact frame index
                fixation_2.tStart = t  # local t and not account for scr refresh
                fixation_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixation_2, 'tStartRefresh')  # time at next scr refresh
                fixation_2.setAutoDraw(True)
            if fixation_2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixation_2.tStartRefresh + fix2_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    fixation_2.tStop = t  # not accounting for scr refresh
                    fixation_2.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fixation_2, 'tStopRefresh')  # time at next scr refresh
                    fixation_2.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        if resp_trial_1.keys == str(curr_corr):
            fb_color = 'teal'
            fb_txt = 'Yes'
            trial_acc = 1
        elif resp_trial_1.keys == None:
            fb_color = 'white'
            fb_txt = '?'
            trial_acc = 0
        else:
            fb_color = 'salmon'
            fb_txt = 'No'
            trial_acc = 0
        runs.addData('trial_acc', trial_acc)
        
        #color_1 = color_2 = trialbox_color
        
        run_1.addData('isi_bkg_1.started', isi_bkg_1.tStartRefresh)
        run_1.addData('isi_bkg_1.stopped', isi_bkg_1.tStopRefresh)
        run_1.addData('bkg_1.started', bkg_1.tStartRefresh)
        run_1.addData('bkg_1.stopped', bkg_1.tStopRefresh)
        run_1.addData('fixation_1.started', fixation_1.tStartRefresh)
        run_1.addData('fixation_1.stopped', fixation_1.tStopRefresh)
        run_1.addData('stim_img_1.started', stim_img_1.tStartRefresh)
        run_1.addData('stim_img_1.stopped', stim_img_1.tStopRefresh)
        run_1.addData('left_resp_1.started', left_resp_1.tStartRefresh)
        run_1.addData('left_resp_1.stopped', left_resp_1.tStopRefresh)
        run_1.addData('right_resp_1.started', right_resp_1.tStartRefresh)
        run_1.addData('right_resp_1.stopped', right_resp_1.tStopRefresh)
        run_1.addData('left_obj.started', left_obj.tStartRefresh)
        run_1.addData('left_obj.stopped', left_obj.tStopRefresh)
        run_1.addData('right_obj.started', right_obj.tStartRefresh)
        run_1.addData('right_obj.stopped', right_obj.tStopRefresh)
        # check responses
        if resp_trial_1.keys in ['', [], None]:  # No response was made
            resp_trial_1.keys = None
        run_1.addData('resp_trial_1.keys',resp_trial_1.keys)
        if resp_trial_1.keys != None:  # we had a response
            run_1.addData('resp_trial_1.rt', resp_trial_1.rt)
        run_1.addData('resp_trial_1.started', resp_trial_1.tStartRefresh)
        run_1.addData('resp_trial_1.stopped', resp_trial_1.tStopRefresh)
        run_1.addData('cue_1.started', cue_1.tStartRefresh)
        run_1.addData('cue_1.stopped', cue_1.tStopRefresh)
        run_1.addData('fixation_2.started', fixation_2.tStartRefresh)
        run_1.addData('fixation_2.stopped', fixation_2.tStopRefresh)
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # ------Prepare to start Routine "feedback"-------
        continueRoutine = True
        routineTimer.add(0.500000)
        # update component parameters for each repeat
        fb_text.setColor(fb_color, colorSpace='rgb')
        fb_text.setText(fb_txt)
        if run_schedule[i] == -1:
            continueRoutine = False
        else:
            continueRoutine = True
        # keep track of which components have finished
        feedbackComponents = [fb_bkg, fb_text]
        for thisComponent in feedbackComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        feedbackClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "feedback"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = feedbackClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=feedbackClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fb_bkg* updates
            if fb_bkg.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fb_bkg.frameNStart = frameN  # exact frame index
                fb_bkg.tStart = t  # local t and not account for scr refresh
                fb_bkg.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fb_bkg, 'tStartRefresh')  # time at next scr refresh
                fb_bkg.setAutoDraw(True)
            if fb_bkg.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fb_bkg.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fb_bkg.tStop = t  # not accounting for scr refresh
                    fb_bkg.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fb_bkg, 'tStopRefresh')  # time at next scr refresh
                    fb_bkg.setAutoDraw(False)
            
            # *fb_text* updates
            if fb_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fb_text.frameNStart = frameN  # exact frame index
                fb_text.tStart = t  # local t and not account for scr refresh
                fb_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fb_text, 'tStartRefresh')  # time at next scr refresh
                fb_text.setAutoDraw(True)
            if fb_text.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fb_text.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fb_text.tStop = t  # not accounting for scr refresh
                    fb_text.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fb_text, 'tStopRefresh')  # time at next scr refresh
                    fb_text.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in feedbackComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "feedback"-------
        for thisComponent in feedbackComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        run_1.addData('fb_bkg.started', fb_bkg.tStartRefresh)
        run_1.addData('fb_bkg.stopped', fb_bkg.tStopRefresh)
        run_1.addData('fb_text.started', fb_text.tStartRefresh)
        run_1.addData('fb_text.stopped', fb_text.tStopRefresh)
        thisExp.nextEntry()
        
    # completed 180 repeats of 'run_1'
    
    
    # ------Prepare to start Routine "delay_end"-------
    continueRoutine = True
    routineTimer.add(6.000000)
    # update component parameters for each repeat
    # keep track of which components have finished
    delay_endComponents = [image, fix_end]
    for thisComponent in delay_endComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    delay_endClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "delay_end"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = delay_endClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=delay_endClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *image* updates
        if image.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image.frameNStart = frameN  # exact frame index
            image.tStart = t  # local t and not account for scr refresh
            image.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
            image.setAutoDraw(True)
        if image.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image.tStartRefresh + 6-frameTolerance:
                # keep track of stop time/frame for later
                image.tStop = t  # not accounting for scr refresh
                image.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image, 'tStopRefresh')  # time at next scr refresh
                image.setAutoDraw(False)
        
        # *fix_end* updates
        if fix_end.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fix_end.frameNStart = frameN  # exact frame index
            fix_end.tStart = t  # local t and not account for scr refresh
            fix_end.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fix_end, 'tStartRefresh')  # time at next scr refresh
            fix_end.setAutoDraw(True)
        if fix_end.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fix_end.tStartRefresh + 6-frameTolerance:
                # keep track of stop time/frame for later
                fix_end.tStop = t  # not accounting for scr refresh
                fix_end.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fix_end, 'tStopRefresh')  # time at next scr refresh
                fix_end.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in delay_endComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "delay_end"-------
    for thisComponent in delay_endComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    set_loop.addData('image.started', image.tStartRefresh)
    set_loop.addData('image.stopped', image.tStopRefresh)
    set_loop.addData('fix_end.started', fix_end.tStartRefresh)
    set_loop.addData('fix_end.stopped', fix_end.tStopRefresh)
    
    # ------Prepare to start Routine "run_break"-------
    continueRoutine = True
    # update component parameters for each repeat
    current_run = 'run_2'
    # checks the trial order to make sure its balanced
    success = False
    while success == False:
        faces = 0
        run_schedule = run_trialorder()
        for i, item in enumerate(run_schedule):
            if i != 179:
                if item in [0,2] and run_schedule[i+1] == -1:
                    faces += 1
        if faces >= 12 and faces <= 18:
            success = True
    key_run_break.keys = []
    key_run_break.rt = []
    _key_run_break_allKeys = []
    # keep track of which components have finished
    run_breakComponents = [text_run_break, key_run_break]
    for thisComponent in run_breakComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    run_breakClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "run_break"-------
    while continueRoutine:
        # get current time
        t = run_breakClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=run_breakClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_run_break* updates
        if text_run_break.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_run_break.frameNStart = frameN  # exact frame index
            text_run_break.tStart = t  # local t and not account for scr refresh
            text_run_break.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_run_break, 'tStartRefresh')  # time at next scr refresh
            text_run_break.setAutoDraw(True)
        
        # *key_run_break* updates
        waitOnFlip = False
        if key_run_break.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            key_run_break.frameNStart = frameN  # exact frame index
            key_run_break.tStart = t  # local t and not account for scr refresh
            key_run_break.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(key_run_break, 'tStartRefresh')  # time at next scr refresh
            key_run_break.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(key_run_break.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(key_run_break.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if key_run_break.status == STARTED and not waitOnFlip:
            theseKeys = key_run_break.getKeys(keyList=['space', '5'], waitRelease=False)
            _key_run_break_allKeys.extend(theseKeys)
            if len(_key_run_break_allKeys):
                key_run_break.keys = _key_run_break_allKeys[-1].name  # just the last key pressed
                key_run_break.rt = _key_run_break_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in run_breakComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "run_break"-------
    for thisComponent in run_breakComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    expClock.reset()
    set_loop.addData('text_run_break.started', text_run_break.tStartRefresh)
    set_loop.addData('text_run_break.stopped', text_run_break.tStopRefresh)
    # check responses
    if key_run_break.keys in ['', [], None]:  # No response was made
        key_run_break.keys = None
    set_loop.addData('key_run_break.keys',key_run_break.keys)
    if key_run_break.keys != None:  # we had a response
        set_loop.addData('key_run_break.rt', key_run_break.rt)
    set_loop.addData('key_run_break.started', key_run_break.tStartRefresh)
    set_loop.addData('key_run_break.stopped', key_run_break.tStopRefresh)
    # the Routine "run_break" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "delay_init"-------
    continueRoutine = True
    routineTimer.add(5.500000)
    # update component parameters for each repeat
    # keep track of which components have finished
    delay_initComponents = [image_delay, fixation]
    for thisComponent in delay_initComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    delay_initClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "delay_init"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = delay_initClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=delay_initClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *image_delay* updates
        if image_delay.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image_delay.frameNStart = frameN  # exact frame index
            image_delay.tStart = t  # local t and not account for scr refresh
            image_delay.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image_delay, 'tStartRefresh')  # time at next scr refresh
            image_delay.setAutoDraw(True)
        if image_delay.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image_delay.tStartRefresh + 5.5-frameTolerance:
                # keep track of stop time/frame for later
                image_delay.tStop = t  # not accounting for scr refresh
                image_delay.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image_delay, 'tStopRefresh')  # time at next scr refresh
                image_delay.setAutoDraw(False)
        
        # *fixation* updates
        if fixation.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fixation.frameNStart = frameN  # exact frame index
            fixation.tStart = t  # local t and not account for scr refresh
            fixation.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fixation, 'tStartRefresh')  # time at next scr refresh
            fixation.setAutoDraw(True)
        if fixation.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fixation.tStartRefresh + 5.5-frameTolerance:
                # keep track of stop time/frame for later
                fixation.tStop = t  # not accounting for scr refresh
                fixation.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fixation, 'tStopRefresh')  # time at next scr refresh
                fixation.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in delay_initComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "delay_init"-------
    for thisComponent in delay_initComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    set_loop.addData('image_delay.started', image_delay.tStartRefresh)
    set_loop.addData('image_delay.stopped', image_delay.tStopRefresh)
    set_loop.addData('fixation.started', fixation.tStartRefresh)
    set_loop.addData('fixation.stopped', fixation.tStopRefresh)
    
    # set up handler to look after randomisation of conditions etc
    run_2 = data.TrialHandler(nReps=180, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='run_2')
    thisExp.addLoop(run_2)  # add the loop to the experiment
    thisRun_2 = run_2.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRun_2.rgb)
    if thisRun_2 != None:
        for paramName in thisRun_2:
            exec('{} = thisRun_2[paramName]'.format(paramName))
    
    for thisRun_2 in run_2:
        currentLoop = run_2
        # abbreviate parameter names if possible (e.g. rgb = thisRun_2.rgb)
        if thisRun_2 != None:
            for paramName in thisRun_2:
                exec('{} = thisRun_2[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        # update component parameters for each repeat
        runs = eval(current_run)
        runs.addData('StimOnset', expClock.getTime()+0.8)#add length of delay
        exec_run = False
        
        i = runs.thisN #counts the number of routine interation
        if run_schedule[i] == -1: #baseline trial
            trial_background = (np.random.rand(1200,800) > 0.85) * 1 #random static background image
            stim = kaleido_stims[0] #preset opacity of first stim to transparent (0)
            stim_opacity = 0
            object_opacity = 0
            trial_dur = 5.0
            curr_corr = np.random.randint(1,3) #pick random resp box as target
            color_1 = color_2 = 'white' #resp box color
            opacity_1 = opacity_2 = 0.55 #preset both to foil opacity    trial_dur = 8.0
            shape_onset = np.random.uniform(0.8, 3.3)
            shape_dur = 1.7
            resp_onset = shape_onset + 1.0
            resp_dur = 0.7
            fix_dur = shape_onset - 0.1
            fix2_onset = resp_onset + 0.7
            fix2_dur = 5.5 - fix2_onset
            exec('opacity_{0} = 0.65'.format(curr_corr)) #set increased opacity for target
            runs.addData('kaleido_img', 'NA') #if BL trial, no stim image
        else: #task trial
            trial_background = np.ones((1200,800)) * 0.5 #dark gray background
            stim_opacity = 1
            object_opacity = 1
            fix_dur = 1.3
            trial_dur = 2.0
            shape_onset = 0.8
            shape_dur = 1.7 #offset time = 2.5
            resp_onset = 1.8 #0.5 seconds 
            resp_dur = 0.7
            fix2_onset = 0
            fix2_dur = 0
            color_1 = color_2 = 'gray' #inside box is black
            opacity_1 = opacity_2 = 1 #response boxes to solid
            corr_img = None
            incorr_img = None
            #fixed trials
            if run_schedule[i] == 0: #S1F1 - face
                stim = kaleido_stims[0] #set fractal
                corr_img = object_stims[0]
                incorr_img = object_stims[1]
                stim_label = 'face1'
            elif run_schedule[i] == 1: #S1F2 - scene
                stim = kaleido_stims[1]
                corr_img = object_stims[1]
                incorr_img = object_stims[0]
                stim_label = 'scene1'
            elif run_schedule[i] == 2: #S2F1 - face
                stim = kaleido_stims[2]
                corr_img = object_stims[2]
                incorr_img = object_stims[3]
                stim_label = 'face2'
            elif run_schedule[i] == 3: #S2F2 - scene
                stim = kaleido_stims[3]
                corr_img = object_stims[3]
                incorr_img = object_stims[2]
                stim_label = 'scene2'
            #conditional trials
            else:
                stim = kaleido_stims[4]
                if run_schedule[i-1] in (0,2): #previous trial was S1F1 or S2F1 (face)
                    corr_img = random.choice(face_stim)
                    incorr_img = random.choice(scene_stim)
                elif run_schedule[i-1] in (1,3): #previous trial was S1F2 or S2F2 (scene)
                    corr_img = random.choice(scene_stim)
                    incorr_img = random.choice(face_stim)
                stim_label = 'COND'
            left_box_img, right_box_img, curr_corr = sides_assign(corr_img, incorr_img)
            runs.addData('kaleido_img', kaleido_stims[run_schedule[i]][-8:])
            runs.addData('stim_ID', stim_label)
        runs.addData('correct', curr_corr)
        bkg_1.setImage(trial_background)
        stim_img_1.setOpacity(stim_opacity)
        stim_img_1.setImage(stim)
        left_resp_1.setOpacity(opacity_1)
        right_resp_1.setOpacity(opacity_2)
        left_obj.setOpacity(object_opacity)
        left_obj.setImage(left_box_img)
        right_obj.setOpacity(object_opacity)
        right_obj.setImage(right_box_img)
        resp_trial_1.keys = []
        resp_trial_1.rt = []
        _resp_trial_1_allKeys = []
        # keep track of which components have finished
        trialComponents = [isi_bkg_1, bkg_1, fixation_1, stim_img_1, left_resp_1, right_resp_1, left_obj, right_obj, resp_trial_1, cue_1, fixation_2]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            if len(resp_trial_1.keys) != 0 and exec_run == False:
                which_box = resp_trial_1.keys[0]
                exec_run = True
                exec('color_{0} = "yellow"'.format(which_box))
            
            # *isi_bkg_1* updates
            if isi_bkg_1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                isi_bkg_1.frameNStart = frameN  # exact frame index
                isi_bkg_1.tStart = t  # local t and not account for scr refresh
                isi_bkg_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(isi_bkg_1, 'tStartRefresh')  # time at next scr refresh
                isi_bkg_1.setAutoDraw(True)
            if isi_bkg_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > isi_bkg_1.tStartRefresh + 0.6-frameTolerance:
                    # keep track of stop time/frame for later
                    isi_bkg_1.tStop = t  # not accounting for scr refresh
                    isi_bkg_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(isi_bkg_1, 'tStopRefresh')  # time at next scr refresh
                    isi_bkg_1.setAutoDraw(False)
            
            # *bkg_1* updates
            if bkg_1.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
                # keep track of start time/frame for later
                bkg_1.frameNStart = frameN  # exact frame index
                bkg_1.tStart = t  # local t and not account for scr refresh
                bkg_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(bkg_1, 'tStartRefresh')  # time at next scr refresh
                bkg_1.setAutoDraw(True)
            if bkg_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > bkg_1.tStartRefresh + trial_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    bkg_1.tStop = t  # not accounting for scr refresh
                    bkg_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(bkg_1, 'tStopRefresh')  # time at next scr refresh
                    bkg_1.setAutoDraw(False)
            
            # *fixation_1* updates
            if fixation_1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixation_1.frameNStart = frameN  # exact frame index
                fixation_1.tStart = t  # local t and not account for scr refresh
                fixation_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixation_1, 'tStartRefresh')  # time at next scr refresh
                fixation_1.setAutoDraw(True)
            if fixation_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixation_1.tStartRefresh + fix_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    fixation_1.tStop = t  # not accounting for scr refresh
                    fixation_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fixation_1, 'tStopRefresh')  # time at next scr refresh
                    fixation_1.setAutoDraw(False)
            
            # *stim_img_1* updates
            if stim_img_1.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                stim_img_1.frameNStart = frameN  # exact frame index
                stim_img_1.tStart = t  # local t and not account for scr refresh
                stim_img_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(stim_img_1, 'tStartRefresh')  # time at next scr refresh
                stim_img_1.setAutoDraw(True)
            if stim_img_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > stim_img_1.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    stim_img_1.tStop = t  # not accounting for scr refresh
                    stim_img_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(stim_img_1, 'tStopRefresh')  # time at next scr refresh
                    stim_img_1.setAutoDraw(False)
            
            # *left_resp_1* updates
            if left_resp_1.status == NOT_STARTED and tThisFlip >= shape_onset-frameTolerance:
                # keep track of start time/frame for later
                left_resp_1.frameNStart = frameN  # exact frame index
                left_resp_1.tStart = t  # local t and not account for scr refresh
                left_resp_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(left_resp_1, 'tStartRefresh')  # time at next scr refresh
                left_resp_1.setAutoDraw(True)
            if left_resp_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > left_resp_1.tStartRefresh + shape_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    left_resp_1.tStop = t  # not accounting for scr refresh
                    left_resp_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(left_resp_1, 'tStopRefresh')  # time at next scr refresh
                    left_resp_1.setAutoDraw(False)
            if left_resp_1.status == STARTED:  # only update if drawing
                left_resp_1.setFillColor(color_1, log=False)
            
            # *right_resp_1* updates
            if right_resp_1.status == NOT_STARTED and tThisFlip >= shape_onset-frameTolerance:
                # keep track of start time/frame for later
                right_resp_1.frameNStart = frameN  # exact frame index
                right_resp_1.tStart = t  # local t and not account for scr refresh
                right_resp_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(right_resp_1, 'tStartRefresh')  # time at next scr refresh
                right_resp_1.setAutoDraw(True)
            if right_resp_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > right_resp_1.tStartRefresh + shape_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    right_resp_1.tStop = t  # not accounting for scr refresh
                    right_resp_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(right_resp_1, 'tStopRefresh')  # time at next scr refresh
                    right_resp_1.setAutoDraw(False)
            if right_resp_1.status == STARTED:  # only update if drawing
                right_resp_1.setFillColor(color_2, log=False)
            
            # *left_obj* updates
            if left_obj.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                left_obj.frameNStart = frameN  # exact frame index
                left_obj.tStart = t  # local t and not account for scr refresh
                left_obj.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(left_obj, 'tStartRefresh')  # time at next scr refresh
                left_obj.setAutoDraw(True)
            if left_obj.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > left_obj.tStartRefresh + 1.7-frameTolerance:
                    # keep track of stop time/frame for later
                    left_obj.tStop = t  # not accounting for scr refresh
                    left_obj.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(left_obj, 'tStopRefresh')  # time at next scr refresh
                    left_obj.setAutoDraw(False)
            
            # *right_obj* updates
            if right_obj.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
                # keep track of start time/frame for later
                right_obj.frameNStart = frameN  # exact frame index
                right_obj.tStart = t  # local t and not account for scr refresh
                right_obj.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(right_obj, 'tStartRefresh')  # time at next scr refresh
                right_obj.setAutoDraw(True)
            if right_obj.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > right_obj.tStartRefresh + 1.7-frameTolerance:
                    # keep track of stop time/frame for later
                    right_obj.tStop = t  # not accounting for scr refresh
                    right_obj.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(right_obj, 'tStopRefresh')  # time at next scr refresh
                    right_obj.setAutoDraw(False)
            
            # *resp_trial_1* updates
            waitOnFlip = False
            if resp_trial_1.status == NOT_STARTED and tThisFlip >= resp_onset-frameTolerance:
                # keep track of start time/frame for later
                resp_trial_1.frameNStart = frameN  # exact frame index
                resp_trial_1.tStart = t  # local t and not account for scr refresh
                resp_trial_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(resp_trial_1, 'tStartRefresh')  # time at next scr refresh
                resp_trial_1.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(resp_trial_1.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(resp_trial_1.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if resp_trial_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > resp_trial_1.tStartRefresh + resp_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    resp_trial_1.tStop = t  # not accounting for scr refresh
                    resp_trial_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(resp_trial_1, 'tStopRefresh')  # time at next scr refresh
                    resp_trial_1.status = FINISHED
            if resp_trial_1.status == STARTED and not waitOnFlip:
                theseKeys = resp_trial_1.getKeys(keyList=['1', '2'], waitRelease=False)
                _resp_trial_1_allKeys.extend(theseKeys)
                if len(_resp_trial_1_allKeys):
                    resp_trial_1.keys = _resp_trial_1_allKeys[0].name  # just the first key pressed
                    resp_trial_1.rt = _resp_trial_1_allKeys[0].rt
            
            # *cue_1* updates
            if cue_1.status == NOT_STARTED and tThisFlip >= resp_onset-frameTolerance:
                # keep track of start time/frame for later
                cue_1.frameNStart = frameN  # exact frame index
                cue_1.tStart = t  # local t and not account for scr refresh
                cue_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(cue_1, 'tStartRefresh')  # time at next scr refresh
                cue_1.setAutoDraw(True)
            if cue_1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > cue_1.tStartRefresh + resp_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    cue_1.tStop = t  # not accounting for scr refresh
                    cue_1.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(cue_1, 'tStopRefresh')  # time at next scr refresh
                    cue_1.setAutoDraw(False)
            
            # *fixation_2* updates
            if fixation_2.status == NOT_STARTED and tThisFlip >= fix2_onset-frameTolerance:
                # keep track of start time/frame for later
                fixation_2.frameNStart = frameN  # exact frame index
                fixation_2.tStart = t  # local t and not account for scr refresh
                fixation_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixation_2, 'tStartRefresh')  # time at next scr refresh
                fixation_2.setAutoDraw(True)
            if fixation_2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixation_2.tStartRefresh + fix2_dur-frameTolerance:
                    # keep track of stop time/frame for later
                    fixation_2.tStop = t  # not accounting for scr refresh
                    fixation_2.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fixation_2, 'tStopRefresh')  # time at next scr refresh
                    fixation_2.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        if resp_trial_1.keys == str(curr_corr):
            fb_color = 'teal'
            fb_txt = 'Yes'
            trial_acc = 1
        elif resp_trial_1.keys == None:
            fb_color = 'white'
            fb_txt = '?'
            trial_acc = 0
        else:
            fb_color = 'salmon'
            fb_txt = 'No'
            trial_acc = 0
        runs.addData('trial_acc', trial_acc)
        
        #color_1 = color_2 = trialbox_color
        
        run_2.addData('isi_bkg_1.started', isi_bkg_1.tStartRefresh)
        run_2.addData('isi_bkg_1.stopped', isi_bkg_1.tStopRefresh)
        run_2.addData('bkg_1.started', bkg_1.tStartRefresh)
        run_2.addData('bkg_1.stopped', bkg_1.tStopRefresh)
        run_2.addData('fixation_1.started', fixation_1.tStartRefresh)
        run_2.addData('fixation_1.stopped', fixation_1.tStopRefresh)
        run_2.addData('stim_img_1.started', stim_img_1.tStartRefresh)
        run_2.addData('stim_img_1.stopped', stim_img_1.tStopRefresh)
        run_2.addData('left_resp_1.started', left_resp_1.tStartRefresh)
        run_2.addData('left_resp_1.stopped', left_resp_1.tStopRefresh)
        run_2.addData('right_resp_1.started', right_resp_1.tStartRefresh)
        run_2.addData('right_resp_1.stopped', right_resp_1.tStopRefresh)
        run_2.addData('left_obj.started', left_obj.tStartRefresh)
        run_2.addData('left_obj.stopped', left_obj.tStopRefresh)
        run_2.addData('right_obj.started', right_obj.tStartRefresh)
        run_2.addData('right_obj.stopped', right_obj.tStopRefresh)
        # check responses
        if resp_trial_1.keys in ['', [], None]:  # No response was made
            resp_trial_1.keys = None
        run_2.addData('resp_trial_1.keys',resp_trial_1.keys)
        if resp_trial_1.keys != None:  # we had a response
            run_2.addData('resp_trial_1.rt', resp_trial_1.rt)
        run_2.addData('resp_trial_1.started', resp_trial_1.tStartRefresh)
        run_2.addData('resp_trial_1.stopped', resp_trial_1.tStopRefresh)
        run_2.addData('cue_1.started', cue_1.tStartRefresh)
        run_2.addData('cue_1.stopped', cue_1.tStopRefresh)
        run_2.addData('fixation_2.started', fixation_2.tStartRefresh)
        run_2.addData('fixation_2.stopped', fixation_2.tStopRefresh)
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # ------Prepare to start Routine "feedback"-------
        continueRoutine = True
        routineTimer.add(0.500000)
        # update component parameters for each repeat
        fb_text.setColor(fb_color, colorSpace='rgb')
        fb_text.setText(fb_txt)
        if run_schedule[i] == -1:
            continueRoutine = False
        else:
            continueRoutine = True
        # keep track of which components have finished
        feedbackComponents = [fb_bkg, fb_text]
        for thisComponent in feedbackComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        feedbackClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "feedback"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = feedbackClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=feedbackClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fb_bkg* updates
            if fb_bkg.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fb_bkg.frameNStart = frameN  # exact frame index
                fb_bkg.tStart = t  # local t and not account for scr refresh
                fb_bkg.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fb_bkg, 'tStartRefresh')  # time at next scr refresh
                fb_bkg.setAutoDraw(True)
            if fb_bkg.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fb_bkg.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fb_bkg.tStop = t  # not accounting for scr refresh
                    fb_bkg.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fb_bkg, 'tStopRefresh')  # time at next scr refresh
                    fb_bkg.setAutoDraw(False)
            
            # *fb_text* updates
            if fb_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fb_text.frameNStart = frameN  # exact frame index
                fb_text.tStart = t  # local t and not account for scr refresh
                fb_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fb_text, 'tStartRefresh')  # time at next scr refresh
                fb_text.setAutoDraw(True)
            if fb_text.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fb_text.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fb_text.tStop = t  # not accounting for scr refresh
                    fb_text.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fb_text, 'tStopRefresh')  # time at next scr refresh
                    fb_text.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in feedbackComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "feedback"-------
        for thisComponent in feedbackComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        run_2.addData('fb_bkg.started', fb_bkg.tStartRefresh)
        run_2.addData('fb_bkg.stopped', fb_bkg.tStopRefresh)
        run_2.addData('fb_text.started', fb_text.tStartRefresh)
        run_2.addData('fb_text.stopped', fb_text.tStopRefresh)
        thisExp.nextEntry()
        
    # completed 180 repeats of 'run_2'
    
    
    # ------Prepare to start Routine "delay_end"-------
    continueRoutine = True
    routineTimer.add(6.000000)
    # update component parameters for each repeat
    # keep track of which components have finished
    delay_endComponents = [image, fix_end]
    for thisComponent in delay_endComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    delay_endClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "delay_end"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = delay_endClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=delay_endClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *image* updates
        if image.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image.frameNStart = frameN  # exact frame index
            image.tStart = t  # local t and not account for scr refresh
            image.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
            image.setAutoDraw(True)
        if image.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image.tStartRefresh + 6-frameTolerance:
                # keep track of stop time/frame for later
                image.tStop = t  # not accounting for scr refresh
                image.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image, 'tStopRefresh')  # time at next scr refresh
                image.setAutoDraw(False)
        
        # *fix_end* updates
        if fix_end.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fix_end.frameNStart = frameN  # exact frame index
            fix_end.tStart = t  # local t and not account for scr refresh
            fix_end.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fix_end, 'tStartRefresh')  # time at next scr refresh
            fix_end.setAutoDraw(True)
        if fix_end.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fix_end.tStartRefresh + 6-frameTolerance:
                # keep track of stop time/frame for later
                fix_end.tStop = t  # not accounting for scr refresh
                fix_end.frameNStop = frameN  # exact frame index
                win.timeOnFlip(fix_end, 'tStopRefresh')  # time at next scr refresh
                fix_end.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in delay_endComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "delay_end"-------
    for thisComponent in delay_endComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    set_loop.addData('image.started', image.tStartRefresh)
    set_loop.addData('image.stopped', image.tStopRefresh)
    set_loop.addData('fix_end.started', fix_end.tStartRefresh)
    set_loop.addData('fix_end.stopped', fix_end.tStopRefresh)
    
    # ------Prepare to start Routine "set_break"-------
    continueRoutine = True
    # update component parameters for each repeat
    current_run = 'run_2'
    
    if set_loop.thisN <= 0:
        loop_text2 = 'Please take a short break. \n\nPlease wait for the scanner to begin.'
    else:
        loop_text2 = 'You have completed the final run. \n\nPlease remain still for the next scanning session.'
    text_set_break.setText(loop_text2)
    set_break_key.keys = []
    set_break_key.rt = []
    _set_break_key_allKeys = []
    # keep track of which components have finished
    set_breakComponents = [text_set_break, set_break_key]
    for thisComponent in set_breakComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    set_breakClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "set_break"-------
    while continueRoutine:
        # get current time
        t = set_breakClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=set_breakClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_set_break* updates
        if text_set_break.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_set_break.frameNStart = frameN  # exact frame index
            text_set_break.tStart = t  # local t and not account for scr refresh
            text_set_break.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_set_break, 'tStartRefresh')  # time at next scr refresh
            text_set_break.setAutoDraw(True)
        
        # *set_break_key* updates
        waitOnFlip = False
        if set_break_key.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            set_break_key.frameNStart = frameN  # exact frame index
            set_break_key.tStart = t  # local t and not account for scr refresh
            set_break_key.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(set_break_key, 'tStartRefresh')  # time at next scr refresh
            set_break_key.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(set_break_key.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(set_break_key.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if set_break_key.status == STARTED and not waitOnFlip:
            theseKeys = set_break_key.getKeys(keyList=['space', '5'], waitRelease=False)
            _set_break_key_allKeys.extend(theseKeys)
            if len(_set_break_key_allKeys):
                set_break_key.keys = _set_break_key_allKeys[-1].name  # just the last key pressed
                set_break_key.rt = _set_break_key_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in set_breakComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "set_break"-------
    for thisComponent in set_breakComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    expClock.reset()
    set_loop.addData('text_set_break.started', text_set_break.tStartRefresh)
    set_loop.addData('text_set_break.stopped', text_set_break.tStopRefresh)
    # check responses
    if set_break_key.keys in ['', [], None]:  # No response was made
        set_break_key.keys = None
    set_loop.addData('set_break_key.keys',set_break_key.keys)
    if set_break_key.keys != None:  # we had a response
        set_loop.addData('set_break_key.rt', set_break_key.rt)
    set_loop.addData('set_break_key.started', set_break_key.tStartRefresh)
    set_loop.addData('set_break_key.stopped', set_break_key.tStopRefresh)
    # the Routine "set_break" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 1 repeats of 'set_loop'


# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
