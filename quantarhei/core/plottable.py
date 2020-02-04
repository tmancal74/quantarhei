# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


class Plottable():
    
    
    def __init__(self):
        
        pass
    
    
    def get_figure(self):
        
        fig = plt.figure(constrained_layout=True)
        return fig