'''
 *  \brief      Plot of GP results
 *              
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HyPro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. '''

import re
from matplotlib import pyplot, cm
import matplotlib
import os
import numpy as np
from PIL.ImageOps import fit
import sys, getopt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def readGPdata(fileName):
    with open(fileName, 'r') as GPfile:
        data = []
        for row in GPfile:
            g = re.match('(^[0-9]*) GP: (.*)$', row)
            thisGP = {'id':g.group(1), 'genome':g.group(2)}
            
            row = GPfile.next()
            g = re.match('^[ \t]*([0-9.eE+\-inf]*) (.*)$', row)
            thisGP['fitness'] = float(g.group(1))
            thisGP['message'] = g.group(2)
            
            data.append(thisGP)
                
    return data

def plotGP(testDir, h=None, col=None, zoomVal=0.1):
    GPfiles = [f for f in os.listdir(testDir) if re.match('^Generation-[0-9]*\.dat$', f)]
    GPfiles.sort(key=lambda fil: int(re.match('^Generation-([0-9]*)\.dat$', fil).group(1)))
    
    I = [int(re.match('^Generation-([0-9]*)\.dat$', fil).group(1)) for fil in GPfiles]
    minFit = []
    aveFit = []
    sigFit = []
    fitHist = [[], []]
    infFrac = []
    i = 0
    for f in GPfiles:
        gen = readGPdata(os.path.join(testDir, f))
        fit = np.array([g['fitness'] for g in gen if g['fitness']!=float('Inf')])
        minFit.append(min(fit))
        aveFit.append(np.mean(fit))
        sigFit.append(np.std(fit))
        infFrac.append(1 - float(len(fit))/float(len(gen)))
        
        fitHist[0] = np.concatenate((fitHist[0], fit))
        fitHist[1] = np.concatenate((fitHist[1], [i]*len(fit)))
        
        i += 1
      
    if h is None:
        pyplot.figure()  
        pyplot.hist2d(np.log10(fitHist[0]), fitHist[1], bins=(100, len(I)))
        
        pyplot.figure()  
        pyplot.hist2d(fitHist[0][np.where(fitHist[0]<zoomVal)],
                      fitHist[1][np.where(fitHist[0]<zoomVal)], bins=(100, len(I)))
    
        pyplot.figure()
        pyplot.plot(infFrac)
        
        pyplot.figure()
        h = pyplot.axes()
    if col:
        h.plot(minFit, label=testDir, color=col)
    else:
        h.plot(minFit, label=testDir)
        
    return minFit

if __name__ == '__main__':
    testDir = None
    zoomVal = 0.1
    paper = False
    helpStr = 'plotGP.py [-h|--help] [-t|--test <test directory>] [-z|--zoom <zoom value>] [--paper]'

    try:
        opts, args = getopt.getopt(sys.argv[1:],"ht:z:",["test=","zoom=","paper"])
    except getopt.GetoptError:
        print helpStr
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
           print helpStr
           sys.exit()
        elif opt in ("-t", "--test"):
           testDir = arg
        elif opt in ("-z", "--zoom"):
           zoomVal = float(arg)
        elif opt == '--paper':
           paper = True
    
    if testDir:
        plotGP(testDir, zoomVal=zoomVal)
    else:
        h = pyplot.axes()
        
        GPtests = [f for f in os.listdir('.') if re.match('^Test-[0-9]*$', f)]
        GPtests.sort(key=lambda fil: int(re.match('^Test-([0-9]*)$', fil).group(1)))
        
        fitMin = []
        for f, i in zip(GPtests, range(len(GPtests))):
            if paper:
                if i==0:
                    color = 'k'
                else:
                    color = 'gray'
            else:
                color = cm.get_cmap('jet')(float(i+1)/11.0)
            fit = plotGP(f, h, color)
            fitMin.append(fit[-1])
        
        table = zip(GPtests, fitMin)
        table.sort(key=lambda f: f[1])
        for tb in table:
            print tb[0] + '\t' + repr(tb[1])
    
    pyplot.xlabel('Generation #')
    pyplot.ylabel('Fitness [1/s]')
    pyplot.grid(True)
    pyplot.legend()
    pyplot.show()
        