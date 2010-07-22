#!/usr/bin/env python
import os, sys
import Tkinter
import Tkinter as ttk
from Tkinter import N,E,S,W
import ScrolledText
import tkSimpleDialog
import tkFileDialog, tkMessageBox
import re
import subprocess
import glob
import pyfits
import datetime
import time

import logging


import gpiifs


#from IPython.Debugger import Tracer; stop = Tracer()


# Helpful flag:
DEBUG=False
#if DEBUG:
    #filenames=['pipeline_'+version+'.zip'] 

####################################


####################################

class GPIMechanism(object):
    """ a wrapper class for a mechanism. Called to create each mech part of the GUI """
    def __init__(self, root, parent, name, positions, axis=0, keyword="KEYWORD",**formatting):
        self.name = name
        self.axis=axis
        self.keyword=keyword

        self.frame = ttk.Frame(root, padx=2, pady=2) #padding=(10,10,10,10)) #, padx=10, pady=10)
        ttk.Label(self.frame,text=" "+self.name+":    ", **formatting).grid(row=0,column=0)

        self.value=' -Unknown- '
        if positions is None:
            self.type='Continuous'
            self.entry = ttk.Entry(self.frame, width=8 , **formatting)
            self.entry.insert(0,'0')
            self.entry.grid(row=0,column=1)

        else:  # wheel or stage
            self.type='Discrete'
            self.positions = positions # this should be a hash of values
            self.pos_names = self.positions.keys()

            self.list=Tkinter.Listbox(self.frame, height=len(self.pos_names), **formatting)
            for s in self.pos_names:
                self.list.insert(Tkinter.END, s)
            self.list.grid(row=0,column=1)

        ttk.Button(self.frame, text="Move", command=self.move, **formatting).grid(row=0,column=2)
        self.pos=ttk.Label(self.frame,text=self.value, **formatting)
        self.pos.grid(row=0,column=3)

        self.parent=parent


        
    def move(self):
        if self.type == 'Discrete':
            sel = self.list.curselection()
            if len(sel) == 0:
                self.parent.log( "Nothing selected!")
            else:
                val = self.pos_names[int(self.list.curselection()[0])]
    
                pos = self.positions[val]
                self.parent.log( "MOVE! motor: %s \t position %d \titem: %s " % (self.name,pos, val))
                self.value=val
                self.pos.config(text="%10s" % (val))
        else:
            pos = int(self.entry.get())
            self.value = '%d' % pos
            self.pos.config(text="%10d" % (pos))
            self.parent.log("MOVE! motor: %s \tposition %d" % (self.name, pos))

        self.parent.runcmd('$TLC_ROOT/scripts/gpMcdMove.csh 1 $MACHINE_NAME %d %d' % (self.axis, pos) )


    def printkey(self):
        self.parent.log("%08s = '%s'" % (self.keyword, self.value))

    def updatekey(self):
        return (self.keyword, self.value)

####################################

class GPI_TempGUI(object):
    def __init__(self, logger=None):

        self.DEBUG=DEBUG


        self._make_widgets()

        self.oldfitslist=[]

        self.IFSController = gpiifs.IFSController(parent=self)

        self.datadir = self.IFSController.datadir # '/net/hydrogen/data/projects/gpi/Detector/TestData'+os.sep+self.entry_dir.get()
        #if not os.path.isdir(self.datadir):
        #    self.log("Creating directory "+self.datadir)
        #    os.mkdir(self.datadir)
        #    os.chmod(self.datadir, 0755)
        #os.chdir(self.datadir)

        self.file_logger=None
        self.logger = ViewLog(self.root, geometry='550x600+650+50',logpath=self.datadir)

        self.log('Starting the GPI IFS Temporary GUI')
        self.log('Informative log messages will appear in this window.')
        self.log("    Log written to file "+self.logger.filename)
        self.log(' ')

        self._set_initial_filename()
        self.watchid=None
        self.watchDir(init=True)


        self.root.protocol("WM_DELETE_WINDOW", self.quit)
        self.root.update()
        
    def _make_widgets(self):

        self.root = Tkinter.Tk()
        self.root.title("GPITempGUI")
        self.root.geometry('+100+50')
        formatting = {"background": "#CCCCCC"}
        formatting={}
        frame = ttk.Frame(self.root, padx=10, pady=10) #padding=(10,10,10,10)) #, padx=10, pady=10)

        r=0
        ttk.Label(frame, text="GPI IFS Temporary GUI", font="Helvetica 20 bold italic",  **formatting).grid(row=r, columnspan=4)
        r=r+1

        #----
        mframe = ttk.LabelFrame(frame, text="Motors:", padx=2, pady=2, bd=2, relief=Tkinter.GROOVE) #padding=(10,10,10,10)) #, padx=10, pady=10)

        ttk.Label(mframe, justify=Tkinter.LEFT, text=" Click on desired position, then press MOVE. ", 
              **formatting).grid(row=0, columnspan=3, stick=W+E)

        self.motors=[]
        
        positions = {"Y": 800, "J":400, "H": 00, "K1":1200, "K2":1600}
        self.motors.append( GPIMechanism(mframe, self, "   Filter", positions, axis=3, keyword='FILTER', **formatting) )

        positions = {"Spectral": 9200, "Wollaston":00, "None": 4600}
        self.motors.append( GPIMechanism(mframe, self, "    Prism", positions, axis=4, keyword='PRISM', **formatting) )
        
        positions = {"Inserted": -300, "Removed":3400}
        self.motors.append(  GPIMechanism(mframe, self, "PupilCam", positions, axis=2, keyword='PUPILMIR', **formatting) )

        positions = {"L1": 00, "Dark":200}
        self.motors.append(  GPIMechanism(mframe, self, "%15s" % ("Lyot"), positions, axis=1, keyword='LYOT', **formatting) )

        self.motors.append(  GPIMechanism(mframe, self, "%15s" % ("Focus"), None, axis=5, keyword='FOCUS', **formatting) )
        
        for i in range(len(self.motors)):
            self.motors[i].frame.grid(row=i+1, columnspan=4, stick=W)

        mframe.grid(row=r,stick=E+W)
        r=r+1
        #----
        mframe = ttk.LabelFrame(frame, text="Header Keywords:", padx=2, pady=2, bd=2, relief=Tkinter.GROOVE) #padding=(10,10,10,10)) #, padx=10, pady=10)

                                                                            
        mr=0
        ttk.Label(mframe,text="TARGET:",  **formatting).grid(row=mr,column=1, stick=N+E+S+W)
        self.entry_target = ttk.Entry(mframe, )
        self.entry_target.grid(row=mr,column=2,stick=N+E+S+W)
        mr=mr+1
   
        ttk.Label(mframe,text="COMMENTS:",  **formatting).grid(row=mr,column=1, stick=N+E+S+W)
        self.entry_comment = ttk.Entry(mframe, )
        self.entry_comment.grid(row=mr,column=2,stick=N+E+S+W)
        mr=mr+1
   
        ttk.Label(mframe,text="OBSERVER:",  **formatting).grid(row=mr,column=1, stick=N+E+S+W)
        self.entry_obs= ttk.Entry(mframe, )
        self.entry_obs.grid(row=mr,column=2,stick=N+E+S+W)
        mr=mr+1
 

        mframe.grid(row=r,stick=E+W) 
        r=r+1
        #----

        mframe = ttk.LabelFrame(frame, text="Actions:", padx=2, pady=2, bd=2, relief=Tkinter.GROOVE) #padding=(10,10,10,10)) #, padx=10, pady=10)

        mr=0

        ttk.Label(mframe, justify=Tkinter.LEFT, text=" Modes: 1=single, 2=CDS, 3=MCDS, 4=UTR ", 
              **formatting).grid(row=mr, columnspan=3, stick=W+E)
        mr=mr+1

        ttk.Label(mframe,text="Mode nreads ncoadds:",  **formatting).grid(row=mr,column=1, stick=W)
        self.entry_mode = ttk.Entry(mframe,  width=10, )
        self.entry_mode.grid(row=mr,column=2,stick=N+E+S+W)
        self.entry_mode.insert(0,"2 2 1")
        mr=mr+1


        ttk.Label(mframe,text="ITIME:",  **formatting).grid(row=mr,column=1, stick=W)
        self.entry_itime = ttk.Entry(mframe,  width=10, )
        self.entry_itime.grid(row=mr,column=2,stick=N+E+S+W)
        self.entry_itime.insert(0,"1.5")
        mr=mr+1

        ttk.Label(mframe,text=" ",  **formatting).grid(row=mr,column=1, stick=W)
        mr=mr+1

        ttk.Label(mframe,text="dir on Y:\\\\",  **formatting).grid(row=mr,column=1, stick=W)
        self.entry_dir = ttk.Entry(mframe,  width=10, )
        self.entry_dir.grid(row=mr,column=2,stick=N+E+S+W)
        today = datetime.date.today()
        self.entry_dir.insert(0,"%02d%02d%02d" % (today.year%100, today.month,today.day))
        mr=mr+1

        ttk.Label(mframe,text="filename:",  **formatting).grid(row=mr,column=1, stick=W)
        self.entry_fn = ttk.Entry(mframe,  width=10, )
        self.entry_fn.grid(row=mr,column=2,stick=N+E+S+W)
        self.entry_fn.insert(0,'test0001')
        ttk.Label(mframe,text=".fits",  **formatting).grid(row=mr,column=3, stick=N+E+S+W)
        mr=mr+1
 
        buttonbar = ttk.Frame(mframe, **formatting)
        buttonbar.grid(row=mr,column=0, columnspan=3)
        ttk.Button(buttonbar, text="Start & Init", command=self.startSW, **formatting).grid(row=0,column=0)
        ttk.Button(buttonbar, text="Configure", command=self.configDet, **formatting).grid(row=0,column=1)
        ttk.Button(buttonbar, text="Take Image", command=self.takeImage, **formatting).grid(row=0,column=2)
        ttk.Button(buttonbar, text="Print Keywords", command=self.printKeywords, **formatting).grid(row=0,column=3)
        ttk.Button(buttonbar, text="GMB", command=self.checkGMB, **formatting).grid(row=0,column=4)
        ttk.Button(buttonbar, text="Quit", command=self.quit, **formatting).grid(row=0,column=5)

        mframe.grid(row=r,stick=E+W) 
        r=r+1

        #----

        frame.grid(row=0) 
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def mainloop(self):
        self.root.mainloop()
    
    def quit(self):
        if  tkMessageBox.askokcancel(title="Confirm Quit",message="Are you sure you want to quit?"):
            self.log(' User quit GPITempGUI.py ')
            self.log(' ' )

            # NO - don't do this - leave the detector server running so we can just restart the GUI
            #self.runcmd('gpIfDetector_tester localhost shutdown 1')

            if self.watchid is not None:
                self.root.after_cancel(self.watchid)
            #self.logger.root.destroy()
            del self.logger
            self.root.destroy()


    def startSW(self):
        if  tkMessageBox.askokcancel(title="Confirm Start & Init",message="Are you sure you want to start the servers and initialize the hardware?"):
            self.IFSController.initialize()


    def configDet(self):
        parts = self.entry_mode.get().split()

        self.IFSController.configureDetector( float(self.entry_itime.get())*1000000., parts[0], parts[1], parts[2] )
        #self.runcmd('gpIfDetector_tester localhost configureExposure 0 %d %s' % (  float(self.entry_itime.get())*1000000., ) )

    def takeImage(self):

        self.IFSController.takeExposure(self.entry_fn.get() )
        #outpath = "Y:\\\\%s\\%s.fits" % (self.entry_dir.get(), self.entry_fn.get())
        #self.log("Taking an exposure to %s" % outpath)

        #self.runcmd('gpIfDetector_tester localhost startExposure "%s" 0' % outpath)

        #increment the filename:
        m =  re.search('(.*?)([0-9]+)',self.entry_fn.get())
        if m is not None:
            formatstr = "%0"+str(len(m.group(2)))+"d"
            newfn = m.group(1) + formatstr % ( int(m.group(2))+1)
            self.entry_fn.delete(0,Tkinter.END)
            self.entry_fn.insert(0,newfn)

#
#    def is_camera_exposing(self):
#        """ check GMB to see if exposure is in progress. 
#    returns 1 if camera is exposing, 0 if not.  """
#    
#        vars = ['ifs.observation.cameraReady', 'ifs.observation.exposureInProgress']
#
#        results = []
#    for v in vars: 
#            cmd = '$TLC_ROOT/bin/linux64/gpUtGmbReadVal -sim -connectAs IFS -name '+v
#            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            stdout, stderr = p.communicate()
#            r = re.search(' = <([0-1])>',stdout)
#            if r is not None:
#                results.append(int(r.group(1)))
#                #self.log(" %s = %d" % (v, int(r.group(1))) )
#        return results[1]
#    
 
#    def checkGMB(self):
#        cmd = '$TLC_ROOT/bin/linux64/gpUtGmbReadVal -sim -connectAs IFS -name ifs.observation.cameraReady'
#
#        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        stdout, stderr = p.communicate()
#        self.log('>>    '+cmd)
#        self.log('      '+stdout)
#        r = re.search(' = <([0-1])>',stdout)
#        if r is not None:
#            self.log("    VALUE = %s" % r.group(1)) 
#        
#

    def printKeywords(self):
        keywords = ['TARGET','COMMENTS','OBSERVER']
        entries = [self.entry_target, self.entry_comment, self.entry_obs]
        self.log('-- keywords --')
        for k, e in zip(keywords, entries):
            self.log("%8s = '%s'" % (k, e.get()))
        
        for m in self.motors:
            m.printkey()
        self.log(' ')

    def runcmd(self,cmdstring, **kwargs):
        self.log(">    " +cmdstring)
        sys.stdout.flush()   #make sure to flush output before starting subproc.
        result = subprocess.call(cmdstring, shell=True, **kwargs)
        self.log("     return code= %d " %result)


    def watchDir(self,init=False):
        self.watchid = self.root.after(1000,self.watchDir)


        if self.IFSController.is_camera_exposing(): 
            self.log(' still exposing - waiting.')
            return

        directory = os.curdir
        curlist = set(glob.glob(directory+os.sep+"*.fits"))
        if init:
            self.log("Now watching directory %s for new FITS files." % (self.datadir))
            self.oldfitslist = curlist
            #self.log("Ignoring already present FITS files" +' '.join(curlist))
            self.log("Ignoring %d already present FITS files." % len(curlist))
            return


        newfiles = curlist - self.oldfitslist
        if len(newfiles) > 0:
            self.log("New FITS files!: "+ ' '.join(newfiles))
            for fn in newfiles: self.updateKeywords(fn)

        self.oldfitslist=curlist


    def updateKeywords(self, filename):
        #try:
        if 1:
            # need to watch out for network latency/long compute times meaning that
            # the files are not fully written yet. In that case, wait for the write to finish.
            # FIXME - this code could be more robust, to avoid infinite loop potential etc
            size = os.path.getsize(filename)
            while ((size != 8392320 ) & (size != 16784640)):
                self.log('file %s not fully written to disk yet - waiting.' % filename)
                for i in range(5):
                    self.root.update_idletasks() # process events!
                    time.sleep(0.1)
                size = os.path.getsize(filename)

            
            f = pyfits.open(filename,mode='update')
            if 'TARGET' not in f[0].header.keys():  # check if it's already been updated
                keywords = ['TARGET','COMMENTS','OBSERVER']
                entries = [self.entry_target, self.entry_comment, self.entry_obs]
                f[0].header.add_history(' Header keywords updated by GPITempGUI.py')
                for k, e in zip(keywords, entries):
                    f[0].header.update(k,e.get())
                
                for m in self.motors:
                    f[0].header.update(*m.updatekey() )
                self.log("Updated FITS keywords in %s" % filename) 
            f.close()
     

        #except:
            #self.log("ERROR: could not update FITS keywords in %s" % filename)

    def _set_initial_filename(self):
        i=1
        fn = "test%04d" % i
        while (os.path.exists(self.datadir+os.sep+fn+".fits")):
            i+=1
            fn = "test%04d" % i
        self.entry_fn.delete(0,Tkinter.END)
        self.entry_fn.insert(0,fn)
        self.log("Next filename set to %s.fits"% fn)


    def log(self, message):
        if self.logger is not None:
            self.logger.log(message)

        

#######################################################################     

class ViewLog(tkSimpleDialog.Dialog):
    """ Display log messages of a program """

    def __init__(self, parent, geometry='+950+50',timestep=100, logpath=None):
        """ Create and display window. Log is CumulativeLogger. """

        # open a file logger
        if logpath is None:
            logpath = os.curdir
        today = datetime.date.today()
        log_fn = logpath+os.sep+ "ifslog_%02d%02d%02d.txt" % (today.year%100, today.month,today.day)

        self.file_logger=logging.getLogger('GPITempGUI Log')
        hdlr = logging.FileHandler(log_fn)
        formatter=logging.Formatter('%(asctime)s %(levelname)s\t%(message)s')
        hdlr.setFormatter(formatter)
        self.file_logger.addHandler(hdlr)
        self.file_logger.setLevel(logging.INFO)
        self.file_logger.log(logging.INFO,'----------------------------------------------------------')


        self.hdlr = hdlr # save for deleting later? 

        # create window  for log
        self.root = Tkinter.Tk()
        self.root.geometry(geometry)
        self.root.title("GPITempGUI Log")

        frame = ttk.Frame(self.root)#, padding=(10,10,10,10))
        frame.pack_configure(fill=Tkinter.BOTH, expand=1)
        self.t = ScrolledText.ScrolledText(frame, width=60, height=45) 
        self.t.insert(Tkinter.END, "")
        self.t.configure(state=Tkinter.DISABLED)
        self.t.see(Tkinter.END)
        self.t.pack(fill=Tkinter.BOTH)

        self.currtext=""

        self.root.update()


        self.filename = log_fn  

        #self.timestep=timestep
        #self.aw_alarm = self.root.after(self.timestep,self.checkLog)

        #self.dontclose=True
        #self.root.protocol("WM_DELETE_WINDOW", self.close)

    def close(self):
        if self.dontclose:
            pass
        else:
            self.root.destroy()

    def __del__(self):
        try:
            logging.shutdown()
            self.hdlr.flush()
            self.hdlr.close()
            self.file_logger.removeHandler(self.hdlr)
        except:
            print "problem closing log?"
        #self.dontclose=False
        self.root.destroy()
        print "logger closed"


    def log(self, message, log_level=logging.INFO):
        # write to file
        if self.file_logger is not None:
            self.file_logger.log(log_level, message)
        # and to screen
        self.addText(message)


    def addText(self,newtext):
        self.t.configure(state=Tkinter.NORMAL)
        self.t.insert(Tkinter.END, newtext+"\n")
        self.t.configure(state=Tkinter.DISABLED)
        self.t.see(Tkinter.END)

#######################################################################     

if __name__ == "__main__":

    g = GPI_TempGUI()


    g.mainloop()

    