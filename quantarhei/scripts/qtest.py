# -*- coding: utf-8 -*-

import os
import subprocess
import shlex

def run_unit_test(toexec):

    cmd = "nosetests --nocapture -v "
    return subprocess.call(shlex.split(cmd+toexec))

def run_doc_test(toexec):

    cmd = "nosetests -v --with-doctest "
    return subprocess.call(shlex.split(cmd+toexec))
    

def main():
    
    docpath = "quantarhei"
    unitpath = "tests/unit"
    
    
    #
    # FIXME: test definition will go to a file
    #
    testlist_abs = dict(doctest=["spectroscopy/abscalculator.py"], 
                        unit=["spectroscopy/abs_test.py", 
                              "spectroscopy/abs2_test.py"])
    testlist_prop = dict(doctest=["qm/propagators/rdmpropagator.py"],
                         unit=["qm/propagators/rdmpropagator_test.py",
                               "qm/propagators/test_rdm_propagator_usage.py"])
    
    testlist_lab = dict(doctest=[],
                        unit=["spectroscopy/test_labsetup.py"])
    
    testdef = dict(abscalc=testlist_abs, propag=testlist_prop,
                   labsetup=testlist_lab)
    
    #
    #
    #
    
    print("\n-------------------------------------------------")
    print("qtest script: Running limited functionality test.")
    print("-------------------------------------------------\n")
     
    testname = "propag"
    
    which = ["unit"] #, "doc"]
    
    tests = testdef[testname]
    
    print("Test set:", testname)
    
    failures = 0
    
    if "doc" in which:
        for file in tests["doctest"]:
            
            toexec = os.path.join(docpath,file)
            print("\n*** Running doc tests:", toexec)
            if run_doc_test(toexec) == 0:
                print("\n... doc tests finished 0K") 
            else:
                failures += 1
    
    if failures == 0:
        
        if "unit" in which:
            for file in tests["unit"]:
                
                toexec = os.path.join(unitpath, file)
                print("\n*** Running unit tests:", toexec)
                if run_unit_test(toexec) == 0:
                    print("\n... unit tests finished 0K")
                else:
                    failures += 1
    
    if failures == 0:
        print("\n *** ALL TESTS FINISHED OK ***\n")
    else:
        print("\n *** Finished with", failures, "failure(s) ***\n")