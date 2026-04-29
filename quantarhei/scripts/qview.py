"""Quantarhei file viewer

Loads *.qrp, *.png, and other files and shows their content


"""

from __future__ import annotations

from tkinter import (
    BOTH,
    BOTTOM,
    RIGHT,
    SUNKEN,
    YES,
    Button,
    Frame,
    Label,
    Menu,
    Misc,
    X,
)
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import askyesno, showerror
from typing import Any

# matplotlib.use("TkAgg")
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
# from matplotlib.figure import Figure
import quantarhei as qr


class Viewer(Frame):
    def __init__(
        self, parent: Misc | None = None, numrow: int = 5, numcol: int = 5
    ) -> None:
        Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)
        self.makeWidgets()
        master: Any = self.master
        master.title("Quantarhei File Viewer")

    def makeWidgets(self) -> None:
        self.makeMenuBar()
        self.makeToolBar()
        L = Label(self, text="Quantarhei version " + str(qr.Manager().version))
        L.config(relief=SUNKEN, width=40, height=5, bg="white")
        L.pack(expand=YES, fill=BOTH)

    def makeMenuBar(self) -> None:
        self.menubar = Menu(self.master)
        master: Any = self.master
        master.config(menu=self.menubar)
        self.fileMenu()

    def makeToolBar(self) -> None:
        toolbar = Frame(self, cursor="hand2", relief=SUNKEN, bd=2)
        toolbar.pack(side=BOTTOM, fill=X)
        Button(toolbar, text="Quit", command=self.quit).pack(side=RIGHT)

    def fileMenu(self) -> None:
        pulldown = Menu(self.menubar)
        pulldown.add_command(label="Open...", command=self.onOpen)
        pulldown.add_command(label="Quit", command=self.quit)
        self.menubar.add_cascade(label="File", underline=0, menu=pulldown)

    def notdone(self) -> None:
        showerror("Not implemented", "Functionality not yet available")

    def quit(self) -> None:
        if askyesno("Verify Quit", "Do you really want to quit?"):
            Frame.quit(self)

    def onOpen(self) -> None:
        file = askopenfilename()
        # print(file, loglevel=qr.LOG_REPORT)
        obj = qr.load_parcel(file)
        self.figure(obj)

    def figure(self, obj: Any) -> None:
        from tkinter import TOP

        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg,
            NavigationToolbar2Tk,
        )

        f = plt.figure(figsize=(5, 5), dpi=100)
        # a = f.add_subplot(111)
        # a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])
        with qr.energy_units("1/cm"):
            obj.plot(fig=f, spart=qr.part_ABS)

        canvas = FigureCanvasTkAgg(f, self)
        canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)


def main() -> None:
    # root = Tk()
    # root.title("Open file")
    # Viewer().mainloop()
    pass
