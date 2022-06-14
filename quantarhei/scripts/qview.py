"""

    Quantarhei file viewer

    Loads *.qrp, *.png, and other files and shows their content


"""
from tkinter import Tk
from tkinter import Label, Button, Frame, Menu
from tkinter import YES, BOTH, SUNKEN, BOTTOM, RIGHT, X
from tkinter.messagebox import showerror, askyesno
from tkinter.filedialog import askopenfilename
import matplotlib
#matplotlib.use("TkAgg")
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
#from matplotlib.figure import Figure

import quantarhei as qr
from quantarhei import printlog as print

class Viewer(Frame):
    def __init__(self, parent=None, numrow=5, numcol=5):
        Frame.__init__(self, parent)
        self.pack(expand=YES,fill=BOTH)
        self.makeWidgets()
        self.master.title("Quantarhei File Viewer")


    def makeWidgets(self):
        self.makeMenuBar()
        self.makeToolBar()
        L = Label(self, text="Quantarhei version "+str(qr.Manager().version))
        L.config(relief=SUNKEN, width=40, height=5, bg="white")
        L.pack(expand=YES, fill=BOTH)

    def makeMenuBar(self):
        self.menubar= Menu(self.master)
        self.master.config(menu=self.menubar)
        self.fileMenu()

    def makeToolBar(self):
        toolbar = Frame(self, cursor="hand2", relief=SUNKEN, bd=2)
        toolbar.pack(side=BOTTOM, fill=X)
        Button(toolbar, text="Quit", command=self.quit).pack(side=RIGHT)

    def fileMenu(self):
        pulldown = Menu(self.menubar)
        pulldown.add_command(label="Open...", command=self.onOpen)
        pulldown.add_command(label="Quit", command=self.quit)
        self.menubar.add_cascade(label="File", underline=0, menu=pulldown)

    def notdone(self):
        showerror("Not implemented", "Functionality not yet available")

    def quit(self):
        if askyesno("Verify Quit", "Do you really want to quit?"):
            Frame.quit(self)

    def onOpen(self):
        import matplotlib.pyplot as plt
        file = askopenfilename()
        #print(file, loglevel=qr.LOG_REPORT)
        obj = qr.load_parcel(file)
        self.figure(obj)

    def figure(self, obj):
        import matplotlib.pyplot as plt
        f = plt.figure(figsize=(5,5), dpi=100)
        #a = f.add_subplot(111)
        #a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])
        with qr.energy_units("1/cm"):
            obj.plot(fig=f, spart=qr.part_ABS)

        canvas = FigureCanvasTkAgg(f, self)
        canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)



def main():
    import sys
    #root = Tk()
    #root.title("Open file")
    #Viewer().mainloop()
    pass

