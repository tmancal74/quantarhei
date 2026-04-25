
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


class Plottable:


    def __init__(self) -> None:

        pass


    def get_figure(self) -> Figure:

        fig = plt.figure(constrained_layout=True)
        return fig
