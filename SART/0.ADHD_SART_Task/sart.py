#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.74.00), Mon Jun 23 12:50:59 2014
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import division #so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, logging, gui
from psychopy.constants import * #things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle, seed
from random import choice
import os #handy system and path functions
from itertools import islice
from math import floor

#store info about the experiment session
expName='SART_RT_VAR'#from the Builder filename that created this script
expInfo={'participant':'', 'schedule':'1'}
dlg=gui.DlgFromDict(dictionary=expInfo,title=expName)
if dlg.OK==False: core.quit() #user pressed cancel
expInfo['date']=data.getDateStr()#add a simple timestamp
expInfo['expName']=expName
#setup files for saving
if not os.path.isdir('data'):
    os.makedirs('data') #if this fails (e.g. permissions) we will get error
filename='data' + os.path.sep + '%s_schedule%s' %(expInfo['participant'], expInfo['schedule'])
logFile=logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)#this outputs to the screen, not a file

#an ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
    
sart_datFile=open(filename+'.txt','a')
sart_datFile.write('nRep\tTrial\tStim\tResponse\tRT\tHit\tCorrRej\tCommission\tOmission\n')
probe_datFile=open(filename+'_probe.txt','a')
probe_datFile.write('nRep\tTrial\tStim1\tStim2\tResponse\tRT\n')

#setup the Window
win = visual.Window(size=(1920, 1080), fullscr=True, screen=0, allowGUI=False, allowStencil=False,
    monitor='testMonitor', color='white', colorSpace='rgb')
    
def last_four(seq, curr_trl, n):
    "Returns the last n items of a list"
    if curr_trl < n:
        curr_last_four = seq
        return curr_last_four
    else:
        curr_last_four = seq[curr_trl - n:]
        return curr_last_four

#Initialise components for Routine "waiting"
waitingClock=core.Clock()
scanner=visual.TextStim(win=win, ori=0, name='scanner',
    text='Press Return key to start.',
    font='Arial',
    pos=[0, 0], height=0.1,wrapWidth=None,
    color='black', colorSpace='rgb', opacity=1,
    depth=0.0)

#Initialise components for Routine "fix1"
fix1Clock=core.Clock()
initial_fix=visual.TextStim(win=win, ori=0, name='initial_fix',
    text=u'+',
    font=u'Arial',
    pos=[0, 0], height=0.3,wrapWidth=None,
    color=u'black', colorSpace=u'rgb', opacity=1,
    depth=-4.0)

#Initialise components for Routine "trial"
trialClock=core.Clock()

text=visual.TextStim(win=win, ori=0, name='text',
    text='nonsense',
    font=u'Arial',
    pos=[0, 0], height=0.3,wrapWidth=None,
    color=u'black', colorSpace=u'rgb', opacity=1,
    depth=-2.0)
probe=visual.TextStim(win=win, ori=0, name='probe',
    text='nonsense',
    font=u'Arial',
    pos=[0, 0], height=0.2,wrapWidth=1.5,
    color=u'black', colorSpace=u'rgb', opacity=1,
    depth=-4.0)
fix=visual.TextStim(win=win, ori=0, name='fix',
    text=u'+',
    font=u'Arial',
    pos=[0, 0], height=0.3,wrapWidth=None,
    color=u'black', colorSpace=u'rgb', opacity=1,
    depth=-4.0)

# Create some handy timers
globalClock=core.Clock() #to track the time since experiment started
routineTimer=core.CountdownTimer() #to track time remaining of each (non-slip) routine 

#------Prepare to start Routine"waiting"-------
t=0; waitingClock.reset() #clock 
frameN=-1
#update component parameters for each repeat
sync_pulse = event.BuilderKeyResponse() #create an object of type KeyResponse
sync_pulse.status=NOT_STARTED
#keep track of which components have finished
waitingComponents=[]
waitingComponents.append(scanner)
waitingComponents.append(sync_pulse)
for thisComponent in waitingComponents:
    if hasattr(thisComponent,'status'): thisComponent.status = NOT_STARTED
#-------Start Routine "waiting"-------
continueRoutine=True
while continueRoutine:
    #get current time
    t=waitingClock.getTime()
    frameN=frameN+1#number of completed frames (so 0 in first frame)
    #update/draw components on each frame
    
    #*scanner* updates
    if t>=0.0 and scanner.status==NOT_STARTED:
        #keep track of start time/frame for later
        scanner.tStart=t#underestimates by a little under one frame
        scanner.frameNStart=frameN#exact frame index
        scanner.setAutoDraw(True)
    
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
        theseKeys = event.getKeys(keyList=['+', 'num_add', 'return'])
        if len(theseKeys)>0:#at least one key was pressed
            sync_pulse.keys=theseKeys[-1]#just the last key pressed
            sync_pulse.rt = sync_pulse.clock.getTime()
            #abort routine on response
            continueRoutine=False
    
    #check if all components have finished
    if not continueRoutine: #a component has requested that we end
        routineTimer.reset() #this is the new t0 for non-slip Routines
        break
    continueRoutine=False#will revert to True if at least one component still running
    for thisComponent in waitingComponents:
        if hasattr(thisComponent,"status") and thisComponent.status!=FINISHED:
            continueRoutine=True; break#at least one component has not yet finished
    
    #check for quit (the [Esc] key)
    if event.getKeys(["escape"]):
        core.quit()
    
    #refresh the screen
    if continueRoutine:#don't flip if this routine is over or we'll get a blank screen
        win.flip()

#End of Routine "waiting"
for thisComponent in waitingComponents:
    if hasattr(thisComponent,"setAutoDraw"): thisComponent.setAutoDraw(False)


#------Prepare to start Routine"fix1"-------
t=0; fix1Clock.reset() #clock 
frameN=-1
#update component parameters for each repeat
sync_pulse = event.BuilderKeyResponse() #create an object of type KeyResponse
sync_pulse.status=NOT_STARTED
#keep track of which components have finished
fix1Components=[]
fix1Components.append(initial_fix)
for thisComponent in fix1Components:
    if hasattr(thisComponent,'status'): thisComponent.status = NOT_STARTED
#-------Start Routine "waiting"-------
continueRoutine=True
while continueRoutine:
    #get current time
    t=fix1Clock.getTime()
    frameN=frameN+1#number of completed frames (so 0 in first frame)
    #update/draw components on each frame
    
    #*initial_fix* updates
    if t>=0.0 and initial_fix.status==NOT_STARTED:
        #keep track of start time/frame for later
        initial_fix.tStart=t#underestimates by a little under one frame
        initial_fix.frameNStart=frameN#exact frame index
        initial_fix.setAutoDraw(True)
    elif initial_fix.status==STARTED and t>=(0.0+2):
        initial_fix.setAutoDraw(False)
    
    #check if all components have finished
    if not continueRoutine: #a component has requested that we end
        routineTimer.reset() #this is the new t0 for non-slip Routines
        break
    continueRoutine=False#will revert to True if at least one component still running
    for thisComponent in fix1Components:
        if hasattr(thisComponent,"status") and thisComponent.status!=FINISHED:
            continueRoutine=True; break#at least one component has not yet finished
    
    #check for quit (the [Esc] key)
    if event.getKeys(["escape"]):
        core.quit()
    
    #refresh the screen
    if continueRoutine:#don't flip if this routine is over or we'll get a blank screen
        win.flip()

#End of Routine "fix1"
for thisComponent in fix1Components:
    if hasattr(thisComponent,"setAutoDraw"): thisComponent.setAutoDraw(False)

#set up handler to look after randomisation of conditions etc
schedule = np.genfromtxt('schedule_%s.txt'%expInfo['schedule'],delimiter=', ',dtype=object)
probe_word_pairs = np.genfromtxt('word_pairs.txt', dtype=object)
probe_word_pairs = probe_word_pairs.tolist()
non_targs = [1,2,4,5,6,7,8,9]
stims = []
for i in schedule:
    if i == '0':
        stims.append(choice(non_targs))
    elif i == '1':
        stims.append(3)
    elif i == '2':
        stims.append(choice(non_targs))
        
#stims = range(1,4)+['probe']+['probe2']+range(1,6)
myarray = []
for i in range(len(stims)):
    myarray.append({'stim': stims[i]}) #puts data into an array of dictionaries that the TrialHandler function will accept
trials=data.TrialHandler(nReps=2, method='sequential', 
    extraInfo=expInfo, originPath=None,
    trialList=myarray,
    seed=None, name='trials')
thisExp.addLoop(trials)#add the loop to the experiment
thisTrial=trials.trialList[0]#so we can initialise stimuli with some values
#abbreviate parameter names if possible (e.g. rgb=thisTrial.rgb)
if thisTrial!=None:
    for paramName in thisTrial.keys():
        exec(paramName+'=thisTrial.'+paramName)

all_rts = []
all_mean_rts = []
all_stdev_rts = []
for thisTrial in trials:
    currentLoop = trials
    #abbrieviate parameter names if possible (e.g. rgb=thisTrial.rgb)
    if thisTrial!=None:
        for paramName in thisTrial.keys():
            exec(paramName+'=thisTrial.'+paramName)
    
    #------Prepare to start Routine"trial"-------
    t=0; trialClock.reset() #clock 
    frameN=-1
    #update component parameters for each repeat
    if trials.thisTrialN+1 > 12:
        if all_stdev_rts[-1] > all_stdev_rts[-2]*2.0:
            dur = 4
            curr_prob_words = probe_word_pairs.pop(int(floor(np.random.rand() * len(probe_word_pairs))))
            probe_word1 = curr_prob_words.pop(int(round(np.random.rand())))
            probe_word2 = curr_prob_words[0]
            probe.text = '%s        %s'%(probe_word1, probe_word2)
        else:
            dur = 2
    else:
        dur = 2
    text.text = stim
    key_resp = event.BuilderKeyResponse() #create an object of type KeyResponse
    key_resp.status=NOT_STARTED
    #keep track of which components have finished
    trialComponents=[]
    trialComponents.append(text)
    trialComponents.append(probe)
    trialComponents.append(key_resp)
    trialComponents.append(fix)
    for thisComponent in trialComponents:
        if hasattr(thisComponent,'status'): thisComponent.status = NOT_STARTED
    #-------Start Routine "trial"-------
    continueRoutine=True
    while continueRoutine:
        #get current time
        t=trialClock.getTime()
        frameN=frameN+1#number of completed frames (so 0 in first frame)
        #update/draw components on each frame
        
        #*text* updates
        if trials.thisTrialN+1 > 12:
            if all_stdev_rts[-1] > all_stdev_rts[-2]*2.0:
                stim = 'probe'
                if t>=0.0 and probe.status==NOT_STARTED:
                    #keep track of start time/frame for later
                    probe.tStart=t#underestimates by a little under one frame
                    probe.frameNStart=frameN#exact frame index
                    probe.setAutoDraw(True)
                elif probe.status==STARTED and t>=(0.0+dur):
                    probe.setAutoDraw(False)
                text.status=FINISHED
            else:
                if t>=0.0 and text.status==NOT_STARTED:
                    #keep track of start time/frame for later
                    text.tStart=t#underestimates by a little under one frame
                    text.frameNStart=frameN#exact frame index
                    text.setAutoDraw(True)
                elif text.status==STARTED and t>=(0.0+dur):
                    text.setAutoDraw(False)
                probe.status=FINISHED
        else:
            if t>=0.0 and text.status==NOT_STARTED:
                #keep track of start time/frame for later
                text.tStart=t#underestimates by a little under one frame
                text.frameNStart=frameN#exact frame index
                text.setAutoDraw(True)
            elif text.status==STARTED and t>=(0.0+dur):
                text.setAutoDraw(False)
            probe.status=FINISHED
        
        #*key_resp* updates
        if t>=0.0 and key_resp.status==NOT_STARTED:
            #keep track of start time/frame for later
            key_resp.tStart=t#underestimates by a little under one frame
            key_resp.frameNStart=frameN#exact frame index
            key_resp.status=STARTED
            #keyboard checking is just starting
            key_resp.clock.reset() # now t=0
            event.clearEvents()
        elif key_resp.status==STARTED and t>=(0.0+2):#dur2):
            key_resp.status=STOPPED
        if key_resp.status==STARTED:#only update if being drawn
            theseKeys = event.getKeys(keyList=['1', '2'])
            if len(theseKeys)>0:#at least one key was pressed
                if key_resp.keys==[]:#then this was the first keypress
                    key_resp.keys=theseKeys[0]#just the first key pressed
                    key_resp.rt = key_resp.clock.getTime()
        
        if not isinstance(key_resp.rt, float):
            key_resp.rt = 0.0
            
        #*fix* updates
        if t>=dur and fix.status==NOT_STARTED:
            #keep track of start time/frame for later
            fix.tStart=t#underestimates by a little under one frame
            fix.frameNStart=frameN#exact frame index
            fix.setAutoDraw(True)
        elif fix.status==STARTED and t>=(dur+1.7):
            fix.setAutoDraw(False)
        
        #check if all components have finished
        if not continueRoutine: #a component has requested that we end
            routineTimer.reset() #this is the new t0 for non-slip Routines
            break
        continueRoutine=False#will revert to True if at least one component still running
        for thisComponent in trialComponents:
            if hasattr(thisComponent,"status") and thisComponent.status!=FINISHED:
                continueRoutine=True; break#at least one component has not yet finished
        
        #check for quit (the [Esc] key)
        if event.getKeys(["escape"]):
            core.quit()
        
        #refresh the screen
        if continueRoutine:#don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    #End of Routine "trial"
    for thisComponent in trialComponents:
        if hasattr(thisComponent,"setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    #check responses
    if len(key_resp.keys)==0: #No response was made
       key_resp.keys=None
    #store data for trials (TrialHandler)
    trials.addData('key_resp.keys',key_resp.keys)
    if key_resp.keys != None:#we had a response
        trials.addData('key_resp.rt',key_resp.rt)
    thisExp.nextEntry()
    
    if trials.thisTrialN+1 > 12:
        if all_stdev_rts[-1] > all_stdev_rts[-2]*2.0:
            probe_datFile.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(trials.thisRepN, trials.thisTrialN+1,probe_word1,probe_word2,key_resp.keys,key_resp.rt))
    if stim==3:
        if trials.thisRepN==1:
            key_resp.rt=average(last_four(all_rts, len(all_rts), 4))
            all_rts.append(key_resp.rt)
    else: 
        all_rts.append(key_resp.rt)
    if len(all_rts) >= 4:
        last4_rts = last_four(all_rts, len(all_rts), 4)
        all_mean_rts.append(average(last4_rts))
        all_stdev_rts.append(std(last4_rts))
    print all_stdev_rts
    if stim=='probe' or key_resp.keys=='2':
        hit,cr,com,omis=[0,0,0,0]
    elif stim in non_targs and key_resp.keys=='1':
        hit,cr,com,omis=[1,0,0,0]
    elif stim in non_targs and key_resp.keys==None:
        hit,cr,com,omis=[0,0,0,1]
    elif stim==3 and key_resp.keys=='1':
        hit,cr,com,omis=[0,0,1,0]
    elif stim==3 and key_resp.keys==None:
        hit,cr,com,omis=[0,1,0,0]
    else:
        print 'there was an error with the output for stimulus %s, stim type%s'%(stim,type(stim))
    sart_datFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(trials.thisRepN, trials.thisTrialN+1,stim,key_resp.keys,key_resp.rt,hit,cr,com,omis))

#get names of stimulus parameters
if trials.trialList in ([], [None], None):  params=[]
else:  params = trials.trialList[0].keys()
#save data for this loop
#trials.saveAsPickle(filename+'trials', fileCollisionMethod='rename')
#trials.saveAsExcel(filename+'.xlsx', sheetName='trials',
#    stimOut=params,
#    dataOut=['n','all_mean','all_std', 'all_raw'])



#Shutting down:
win.close()
core.quit()
