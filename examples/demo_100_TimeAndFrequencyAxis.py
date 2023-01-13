# -*- coding: utf-8 -*-

import quantarhei as qr

print("""
      
Example of TimeAxis and FrequencyAxis usage:
    
    Quantarhei defines two time of Time and Frequency Axes. The types are
    `complete` and `upper-half`. See documentation for details of their
    behaviours:
        
    https://quantarhei.readthedocs.io/en/latest/classes/time.html
    https://quantarhei.readthedocs.io/en/latest/classes/frequency.html

      
""")
print("Default behaviour")
print("-----------------")

Ns = [100, 101]

print("\nStarting with TimeAxis:")
for N in Ns:
    print("\nN =", N)
    tmaxs1 = qr.TimeAxis(0.0, N, 10.0)
    print("TimeAxis type       :", tmaxs1.atype)
    print("start, length, step :", tmaxs1.start, tmaxs1.length, tmaxs1.step)
    fraxs1 = tmaxs1.get_FrequencyAxis()
    print("FrequencyAxis type  :",fraxs1.atype)
    print("start, length, step :",fraxs1.start, fraxs1.length, fraxs1.step)

print("Starting with FrequencyAxis:")
for N in Ns:
    print("\nN =", N)
    fraxs2 = qr.FrequencyAxis(0.0, N, 0.01)
    print("FrequencyAxis type  :", fraxs2.atype)
    print("start, length, step :", fraxs2.start, fraxs2.length, fraxs2.step)
    tmaxs2 = fraxs2.get_TimeAxis()
    print("TimeAxis type       :", tmaxs2.atype)
    print("start, length, step :", tmaxs2.start, tmaxs2.length, tmaxs2.step)    

    
print("\n")
print("Non-default behaviour")
print("---------------------")

print("\nStarting with TimeAxis of `complete` type:")
for N in Ns:
    print("\nN =", N)
    tmaxs1 = qr.TimeAxis(0.0, N, 10.0, atype="complete")
    print("TimeAxis type       :", tmaxs1.atype)
    print("start, length, step :", tmaxs1.start, tmaxs1.length, tmaxs1.step)
    fraxs1 = tmaxs1.get_FrequencyAxis()
    print("FrequencyAxis type  :",fraxs1.atype)
    print("start, length, step :",fraxs1.start, fraxs1.length, fraxs1.step)

print("\nStarting with FrequencyAxis of `upper-half` type:")
for N in Ns:
    print("\nN =", N)
    fraxs2 = qr.FrequencyAxis(0.0, N, 0.01, atype="upper-half")

    print("FrequencyAxis type  :", fraxs2.atype)
    print("start, length, step :", fraxs2.start, fraxs2.length, fraxs2.step)
    try:
        tmaxs2 = fraxs2.get_TimeAxis()
        print("TimeAxis type       :", tmaxs2.atype)
    except:
        print("*** Exception occurs when you try to get TimeAxis from "+
              "FrequencyAxis of the upper-half type and odd number of points")
