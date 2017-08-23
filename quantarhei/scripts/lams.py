# -*- coding: utf-8 -*-
#
###############################################################################
#
#
#  lams (Load and manipulate spectra) command line tool
#
#  author: Tomas Mancal
#          Charles University in Prague
#  email : mancal@karlov.mff.cuni.cz
#
#  The script performs some simple manipulations of 
#
#
#  The tool requires Python 3 interpreter and Numpy and Quantarhei packages
#  Use 
#
#       $ lams -h 
#
#  or
#
#       $ lams --help 
#
#  to read about the command line options
#
#

###############################################################################
#
#  Imports from standard library
#
###############################################################################
import os
import sys
import argparse
import glob

try:
    from blessings import Terminal
except:
    print("Blessings package not installed")
    quit()
    
###############################################################################
#
#  Imports from Quantarhei and numerics 
#
###############################################################################
import numpy
import scipy
from quantarhei import AbsSpectContainer, AbsSpectrumBase
from quantarhei import TimeAxis, FrequencyAxis, energy_units

def main():

    ###############################################################################
    #
    #  Parsing command line arguments
    #
    ###############################################################################
    parser = argparse.ArgumentParser(description='Loads and manupilates spectra.')
    #
    #
    #
    parser.add_argument("files", metavar='file', type=str, nargs='+',
                        help='files to be processed')
    parser.add_argument("-o", "--outfile", metavar="outfile", 
                        help="output file", default="tot.dat")
    #parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                    const=sum, default=max,
    #                    help='sum the integers (default: find the max)')
    parser.add_argument("-op", "--operation", default="sum", 
                        choices=["sum", "norm", "max", "convert"], 
                        help="operation to be performed on the spectra")
    parser.add_argument("-t", "--type", default="abs", 
                        help="type of spectra", choices=["abs", "2d"])
    parser.add_argument("-n", "--nosplash", action="store_true",
                        help="no script info (a 'spash screen')")
    
    #
    #
    #
    parser.add_argument("-e", "--exe", help="prints out the location of python"+
                        " executable", action="store_true")
    
    #
    #
    #
    args = parser.parse_args()
    
    ###############################################################################
    #
    #
    #  Results of command line parsing
    #
    #
    ###############################################################################
    files = args.files
    oper  = args.operation
    stype = args.type
    show_exe = args.exe
    splash = not args.nosplash
        
    ###############################################################################
    #
    #
    #  "Splash screen"
    #
    #
    ###############################################################################
    term = Terminal()
    t_height = term.height
    t_width = term.width
    
    class StarLine:
        
        def __init__(self, term):
            self.term = term
            
        def full(self):
            print('*'*self.term.width)
            
        def empty(self):
            print('*'+(' '*(self.term.width-2))+'*')
                
        def center_text(self, text, bold=False):
            tlen = len(text)
            spaces = self.term.width - tlen - 2
            lspaces = spaces//2
            last = self.term.width - (tlen + 2 + 2*lspaces)
            rspaces = lspaces + last
            if bold:
                tx = '*'+(' '*lspaces)+self.term.bold(text)+(' '*rspaces)+'*'
            else:
                tx = '*'+(' '*lspaces)+text+(' '*rspaces)+'*'
            print(tx)
                
        def left_text(self, text):
            tlen = len(text)
            spaces = self.term.width - tlen - 2
            print('*'+text+(' '*spaces)+'*')
            
    indent = 4
            
    if splash:
        star_line = StarLine(term)
        print("")
        star_line.full()
        star_line.empty()
        star_line.center_text("lams (Loading and manipulating spectra)"+
                              " command line tool", bold=True)
        star_line.center_text("Based on Quantarhei package "+
                              "(github.com/tmancal74/quantarhei)")
        star_line.empty()
        star_line.left_text("   Author: Tomas Mancal")
        star_line.left_text("           Charles University")
        star_line.left_text("   Email : mancal@karlov.mff.cuni.cz ")
        star_line.empty()
        star_line.left_text("   try 'lams -h' for basic usage information")
        star_line.empty()
        star_line.full()
        print("")
    
    if show_exe:
        with term.location(indent, term.height-1):
            print("Location of Python executable: ", sys.executable, "\n")
        quit()
        
    with term.location(indent, term.height-1):
        print("Lams task summary:")    
    with term.location(indent, term.height-1):
        print(oper, stype, ": outfile =", args.outfile)
        print("")
        
    fname, ext = os.path.splitext(args.outfile)
    ext = ext[1:]
    
    ###############################################################################
    #
    #  Identification of the files to work on
    #
    ###############################################################################
    files_abs = []
    for aaa in files:
        faaa = glob.glob(aaa)
        for f in faaa:
            files_abs.append(f)
    
            
    abs_sum = (stype == "abs") and (oper == "sum")
    abs_norm = (stype == "abs") and (oper == "norm")
    abs_max = (stype == "abs") and (oper == "max")
    abs_conv = (stype == "abs") and (oper == "convert")
    
    
    if abs_sum:
    
        with term.location(indent, term.height-1):
            print("Summing absorption spectra")
            
        k = 0
        spect = None
        for fl in files_abs:
    
            a = AbsSpectrumBase()
            a.load(fl)
        
            #
            # sum them up
            #
    
            if spect is None:
                spect = a
            else:
                spect.add_to_data(a)
    
            k += 1
    
            
        spect.save(fname, ext=ext)
    
        with term.location(indent, term.height-1):
            print("... finished with ", k, " files\n")    
    
    
    if abs_norm:
    
        with term.location(indent, term.height-1):    
            print("Normalizing absorption spectra")
    
        k = 0
    
        fl = files_abs[0]
        spect = AbsSpectrumBase()
        spect.load(fl)
        
        #
        # Normalize spectra
        #
    
        mx = numpy.max(spect.data)
        spect.data = spect.data/mx
        
        k += 1
    
        spect.save(fname, ext=ext)
    
        with term.location(indent, term.height-1):    
            print("... finished with ", k, " files\n")    
        
    if abs_max:
    
        with term.location(indent, term.height-1):    
            print("Looking for absorption maxima")
    
        k = 0
    
        fl = files_abs[0]
        spect = AbsSpectrumBase()
        spect.load(fl)
        
        #
        # Find location of maximum
        #
    
        mxi = numpy.argmax(spect.data)
        #with energy_units("1/cm"):
        mx = spect.axis.data[mxi]
        
        k += 1
    
        with term.location(indent, term.height-1):    
            print("{t.bold}Maximum at ".format(t=term), mx, "1/cm")
        
        with term.location(indent, term.height-1):    
            print("... finished with ", k, " files\n")    
        
    if abs_conv:
 
        with term.location(indent, term.height-1):    
            print("Converting from nm to 1/cm")
            
        fl = files_abs[0]
        data = numpy.genfromtxt(fl, converters={0: lambda s: 1.0e7/float(s)})
        
        x = data[:,0]
        ab = data[:,1]
        
        _spline_r = scipy.interpolate.UnivariateSpline(x,ab,s=0)
        
        wa = FrequencyAxis(10000.0, 5000, 1.0)
        data = numpy.zeros(wa.data.shape)
        i = 0
        for w in wa.data:
            data[i] = _spline_r(w)
            i += 1
            
        spect = AbsSpectrumBase(axis = wa, data = data)
            
        spect.save(fname, ext=ext)            
        