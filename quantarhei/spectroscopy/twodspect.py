"""The class representing 2D spectrum. It is obtained by the chain:


Option 1) --- Impulsive spectrum only ---

TwoDSpectrumCalculator.calculate()
-> TwoDSpectrumContainer

TwoDSpectrumContainer.get_TwoDSpectrum()
-> TwoDSpectrum

Option 2) --- Allows convoluting with the fields

TwoDResponseCalculator.calculate()
-> TwoDResponseContainer

TwoDResponseContainer.get_TwoDSpectrumContainer  # here the fields get convoluted
-> TwoDSpectrumContainer

TwoDSpectrumContainer.get_TwoDSpectrum
-> TwoDSpectum



TwoDResponseCalculator.

"""

from __future__ import annotations

import warnings
from typing import Any

import matplotlib.pyplot as plt
import numpy

from .. import TWOD_SIGNALS, part_ABS, part_IMAGINARY, part_REAL, signal_TOTL
from ..core.datasaveable import DataSaveable
from ..core.dfunction import DFunction
from ..core.frequency import FrequencyAxis
from ..core.saveable import Saveable
from ..core.valueaxis import ValueAxis
from ..exceptions import QuantarheiError


class TwoDSpectrum(DataSaveable, Saveable):
    """A single two-dimensional Fourier-transform spectrum.

    Stores the 2D spectral data indexed by (omega_1, omega_3) at a fixed
    waiting time ``t2``. The data may be rephasing, non-rephasing, or the
    total sum of both.
    """

    dtypes = TWOD_SIGNALS

    def __init__(self) -> None:

        self.xaxis: ValueAxis | None = None
        self.yaxis: ValueAxis | None = None
        self.data: numpy.ndarray | None = None  # type: ignore[explicit-any]
        self.dtype: str | None = None

        self.t2 = -1.0

        self.params: Any = None  # type: ignore[explicit-any]

    def set_axis_1(self, axis: Any) -> None:  # type: ignore[explicit-any]
        """Sets the x-axis of te spectrum (omega_1 axis)"""
        self.xaxis = axis

    def set_axis_3(self, axis: Any) -> None:  # type: ignore[explicit-any]
        """Sets the y-axis of te spectrum (omega_3 axis)"""
        self.yaxis = axis

    def set_data_type(self, dtype: str = signal_TOTL) -> None:  # "Tot"):
        """Set the data type for this 2D spectrum.

        Parameters
        ----------
        dtype : str, optional
            Type of data stored in this object. Must be one of the values
            in ``TwoDSpectrum.dtypes``. Default is ``qr.signal_TOTL``.

        Raises
        ------
        Exception
            If ``dtype`` is not a recognised data type.
        """
        if dtype in self.dtypes.values():
            self.dtype = dtype
        else:
            raise QuantarheiError("Unknown data type for TwoDSpectrum object")

    def get_spectrum_type(self) -> Any:  # type: ignore[explicit-any]

        return self.dtype

    def set_data(  # type: ignore[explicit-any]
        self, data: numpy.ndarray, dtype: str = signal_TOTL
    ) -> None:  # "Tot"):
        """Set the data of the 2D spectrum.

        Parameters
        ----------
        data : numpy.ndarray
            2D array of spectral values (float or complex).
        dtype : str, optional
            Type of the data stored. Recognised values are ``qr.signal_REPH``
            for rephasing spectra, ``qr.signal_NONR`` for non-rephasing
            spectra, and ``qr.signal_TOTL`` for the total spectrum (sum of
            both). Default is ``qr.signal_TOTL``.

        Raises
        ------
        Exception
            If ``dtype`` is unknown or the axes have not been set yet.
        """
        self.set_data_type(dtype)

        if (self.xaxis is None) or (self.yaxis is None):
            raise QuantarheiError("Axes of the 2D spectrum are not set")

        if (self.xaxis.length == data.shape[0]) and (
            self.yaxis.length == data.shape[1]
        ):
            self.data = data

        else:
            raise QuantarheiError(
                "Axes not compatible with data\n"
                "xaxis.length = "
                + str(self.xaxis.length)
                + "\n"
                + "yaxis.length = "
                + str(self.yaxis.length)
                + "\n"
                + "data shape = "
                + str(data.shape)
            )

    def add_data(self, data: numpy.ndarray) -> None:  # type: ignore[explicit-any]
        """Sets the data of the 2D spectrum"""
        if self.data is None:
            raise QuantarheiError("Data is not initialized: use set_data method.")
        self.data += data

    def overlay_pulses(self, lab: Any) -> None:  # type: ignore[explicit-any]
        """Use labsetup class to overlay pulse spectra over this 2D spectrum"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        # first two pulses are on omega_1 axis
        ome1 = self.xaxis.data
        spect1 = lab.get_pulse_spectrum(0, ome1)
        spect2 = lab.get_pulse_spectrum(1, ome1)

        # third pulse
        ome3 = self.yaxis.data
        spect3 = lab.get_pulse_spectrum(2, ome3)

        self.data = self.data * (
            spect1 * spect2
        )  # multiplying the second (x) axis by E^2
        self.data = (
            self.data * spect3[:, numpy.newaxis]
        )  # multiplying the first (y) axis by E

        # do we need to add also the detection pulse ?

    def set_t2(self, t2: float) -> None:
        """Sets the t2 (waiting time) of the spectrum"""
        self.t2 = t2

    def get_t2(self) -> float:
        """Returns the t2 (waiting time) of the spectrum"""
        return self.t2

    def log_params(self, params: Any = None) -> None:  # type: ignore[explicit-any]
        """Store custom information about the spectrum


        This can be used e.g. for storage of parameters of the system for
        which the spectrum was calculated

        """
        self.params = params

    def get_log_params(self) -> Any:  # type: ignore[explicit-any]
        return self.params

    def get_value_at(self, x: float, y: float) -> Any:  # type: ignore[explicit-any]
        """Returns value of the spectrum at a given coordinate"""
        if self.dtype is None:
            raise QuantarheiError("Data type not set")

        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        (ix, dist) = self.xaxis.locate(x)
        (iy, dist) = self.yaxis.locate(y)

        return self.data[iy, ix]

    def get_cut_along_x(self, y0: float) -> DFunction:
        """Returns a DFunction with the cut of the spectrum along the x axis"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        (iy, dist) = self.yaxis.locate(y0)

        ax = self.xaxis
        vals = numpy.zeros(ax.length, dtype=self.data.dtype)
        for ii in range(ax.length):
            vals[ii] = self.data[iy, ii]

        return DFunction(ax, vals)

    def get_cut_along_y(self, x0: float) -> DFunction:
        """Returns a DFunction with the cut of the spectrum along the y axis"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        (ix, dist) = self.xaxis.locate(x0)

        ay = self.yaxis
        vals = numpy.zeros(ay.length, dtype=self.data.dtype)
        for ii in range(ay.length):
            vals[ii] = self.data[ii, ix]

        return DFunction(ay, vals)

    def get_cut_along_line(
        self,
        point1: list,
        point2: list,
        which_step: str | None = None,
        step: float | None = None,
    ) -> DFunction:
        """Returns a cut along a line specified by two points"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        vx1 = point1[0]
        vy1 = point1[1]
        vx2 = point2[0]
        vy2 = point2[1]

        (x1, dist) = self.xaxis.locate(vx1)
        (y1, dist) = self.yaxis.locate(vy1)
        (x2, dist) = self.xaxis.locate(vx2)
        (y2, dist) = self.yaxis.locate(vy2)

        length = numpy.sqrt((vx1 - vx2) ** 2 + (vy1 - vy2) ** 2)

        if which_step is None:
            wstep = "x"
        else:
            wstep = which_step

        if step is None:
            if wstep == "x":
                dx = self.xaxis.step
            elif wstep == "y":
                dx = self.yaxis.step
            else:
                raise QuantarheiError("Step along the cut line is not defined")

        else:
            dx = step

        Nstep = int(length / dx)

        axis = ValueAxis(0.0, Nstep + 1, dx)

        vals = numpy.zeros(Nstep + 1, dtype=self.data.dtype)
        ii = 0
        for val in axis.data:
            vx1 = self.xaxis.data[x1]
            vx2 = self.xaxis.data[x2]
            vy1 = self.yaxis.data[y1]
            vy2 = self.yaxis.data[y2]
            x = vx1 + val * (vx2 - vx1) / (Nstep * dx)
            y = vy1 + val * (vy2 - vy1) / (Nstep * dx)
            vals[ii] = self.get_value_at(x, y)
            ii += 1

        return DFunction(axis, vals)

    def get_diagonal_cut(self) -> DFunction:
        """Returns cut of the spectrum along the diagonal"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        point1 = [self.xaxis.min, self.yaxis.min]
        point2 = [self.xaxis.max, self.yaxis.max]

        fce = self.get_cut_along_line(point1, point2, which_step="x")

        fce.axis.data += point1[0]

        return fce

    def get_anti_diagonal_cut(self, point: Any) -> None:  # type: ignore[explicit-any]

        pass

    def get_max_value(self, dpart: str = part_REAL) -> float:
        """Maximum value of the real part of the spectrum"""
        assert self.data is not None
        if dpart == part_REAL:
            return numpy.amax(numpy.real(self.data))
        if dpart == part_IMAGINARY:
            return numpy.amax(numpy.imag(self.data))
        if dpart == part_ABS:
            return numpy.amax(numpy.abs(self.data))
        raise QuantarheiError("Unknown data part")

    def get_min_value(self, dpart: str = part_REAL) -> float:
        """Minimum value of the real part of the spectrum"""
        assert self.data is not None
        if dpart == part_REAL:
            return numpy.amin(numpy.real(self.data))
        if dpart == part_IMAGINARY:
            return numpy.min(numpy.imag(self.data))
        if dpart == part_ABS:
            return numpy.amin(numpy.abs(self.data))
        raise QuantarheiError("Unknown data part")

    def get_area_integral(self, area: Any, dpart: str = part_REAL) -> Any:  # type: ignore[explicit-any]
        """Returns an integral of a given area in the 2D spectrum"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None

        def integral_square(  # type: ignore[explicit-any]
            x1: float,
            x2: float,
            y1: float,
            y2: float,
            data: numpy.ndarray,
            dx: float,
            dy: float,
        ) -> float:
            (n1, n2) = data.shape
            data.reshape(n1 * n2)
            return numpy.sum(data) * dy * dy

        area_shape = area[0]
        x1 = area[1][0]
        x2 = area[1][1]
        y1 = area[1][2]
        y2 = area[1][3]

        dx = self.xaxis.step
        dy = self.yaxis.step

        (nx1, derr) = self.xaxis.locate(x1)
        (nx2, derr) = self.xaxis.locate(x2)
        (ny1, derr) = self.yaxis.locate(y1)
        (ny2, derr) = self.yaxis.locate(y2)

        if area_shape == "square":
            int_fce = integral_square
        else:
            raise QuantarheiError("Unknown area type: " + area_shape)

        data = self.data[nx1:nx2, ny1:ny2]

        if dpart == part_REAL:
            return int_fce(x1, x2, y1, y2, numpy.real(data), dx, dy)
        if dpart == part_IMAGINARY:
            return int_fce(x1, x2, y1, y2, numpy.imag(data), dx, dy)
        if dpart == part_ABS:
            return int_fce(x1, x2, y1, y2, numpy.abs(data), dx, dy)
        raise QuantarheiError("Unknown data part")

    def get_area_max(  # type: ignore[explicit-any]
        self, area: Any, dpart: str = part_REAL, loc: list | None = None
    ) -> Any:
        """Returns a max value in a given area in the 2D spectrum"""
        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None

        def find_in_square(  # type: ignore[explicit-any]
            x1: float,
            x2: float,
            y1: float,
            y2: float,
            data: numpy.ndarray,
            dx: float,
            dy: float,
        ) -> float:
            return numpy.amax(data)

        def loc_in_square(  # type: ignore[explicit-any]
            x1: float,
            x2: float,
            y1: float,
            y2: float,
            data: numpy.ndarray,
            dx: float,
            dy: float,
        ) -> tuple:
            loc = numpy.unravel_index(numpy.argmax(data), data.shape)
            x = x1 + loc[1] * dx
            y = y1 + loc[0] * dy
            return (x, y)

        area_shape = area[0]
        x1 = area[1][0]
        x2 = area[1][1]
        y1 = area[1][2]
        y2 = area[1][3]

        dx = self.xaxis.step
        dy = self.yaxis.step

        (nx1, derr) = self.xaxis.locate(x1)
        (nx2, derr) = self.xaxis.locate(x2)
        (ny1, derr) = self.yaxis.locate(y1)
        (ny2, derr) = self.yaxis.locate(y2)

        x1 = self.xaxis.data[nx1]
        x2 = self.xaxis.data[nx2]
        y1 = self.yaxis.data[ny1]
        y2 = self.yaxis.data[ny2]

        if area_shape == "square":
            int_fce = find_in_square
        else:
            raise QuantarheiError("Unknown area type: " + area_shape)

        data = self.data[ny1:ny2, nx1:nx2]

        if dpart == part_REAL:
            if loc is not None:
                loc.append(loc_in_square(x1, x2, y1, y2, data, dx, dy))
            return int_fce(x1, x2, y1, y2, numpy.real(data), dx, dy)
        if dpart == part_IMAGINARY:
            return int_fce(x1, x2, y1, y2, numpy.imag(data), dx, dy)
        if dpart == part_ABS:
            return int_fce(x1, x2, y1, y2, numpy.abs(data), dx, dy)
        raise QuantarheiError("Unknown data part")

    def normalize2(
        self,
        norm: float = 1.0,
        dpart: str = part_REAL,
        nmax: list | None = None,
        use_max: bool = False,
    ) -> None:
        """Normalizes the spectrum to the given maximum.

        Normalizes the spectrum to the maximum of a given part of the spectrum.
        Parts are real, imaginary or absolute values. Normalization is to
        a supplied norm. If `nmax` is not specified, `use_max` is not
        taken into account. If `nmax` is a list-like object, the maximum
        of the current spectrum is appended to it. If `use_max` is set to False
        we use the first element of the list-like object `nmax` as the
        maximum of the present data.


        Parameters
        ----------
        norm : float
            Norm to which we normalize. Default value is 1.

        dpart : string
            Part of the spectrum we normalize against.

        namx : list-like or None
            If `nmax` is a list-like object, the maximum of the current
            spectrum is appended to it.

        use_max : bool
            If `use_max` is set to False we use the first element of the
            list-like object `nmax` as the  maximum of the present data.


        """
        mx = self.get_max_value(dpart=dpart)
        if nmax is not None:
            nmax.append(mx)
            if not use_max:
                mx = nmax[0]
        if mx != 0.0:
            self.data = (self.data / mx) * norm

    def devide_by(self, val: float | int) -> None:
        """Devides the total spectrum by a value


        Parameters
        ----------
        val : float, int
            Value by which we devide the spectrum

        """
        assert self.data is not None
        self.data = self.data / val

    def _interpolate(self) -> None:
        """Interpolate the spectrum by splines"""
        pass

    def zeropad(self, fac: int = 1) -> None:
        if fac == 1:
            return

        assert self.xaxis is not None
        assert self.yaxis is not None
        assert self.data is not None
        Nxa = self.xaxis.length
        Nya = self.yaxis.length

        print("Start: ", self.xaxis.start)
        print("Step:  ", self.xaxis.step)
        print("end:   ", self.xaxis.step * self.xaxis.length + self.xaxis.start)

        nNxa = fac * Nxa
        nNya = fac * Nya

        print("New length x: ", nNxa)
        print("New length y: ", nNya)

    def shift_energy(self, dE: float, interpolation: str = "linear") -> None:
        """Shift the spectrum in both frequency axis by certain amount

        The shift is specified as the energy that has to be added to the
        2D spectrum to shifted to its expected position. As a result
        we have a new spectrum

        S(E_3, E_1) = S(E_3-dE, E_1-dE)

        Because of the negative sign above, we have to change the sign of dE
        as a first thing in the code. This is for user convenience.

        """
        dE = -dE

        ndata = numpy.zeros(self.data.shape, dtype=self.data.dtype)
        ndata1 = numpy.zeros(self.data.shape, dtype=self.data.dtype)

        if interpolation == "linear":
            print("Shifting spectrum: linear interpolation")

            data = self.data

            N1 = self.data.shape[0]
            N2 = self.data.shape[1]

            Dx = self.xaxis.step

            A = numpy.abs(dE)
            n = int(numpy.floor(A / Dx))
            dx = A - n * Dx  # dx is positive

            s = int(numpy.sign(dE))

            # dimesion 1
            # make sure we stay within the defined array
            if s == 1:
                n1 = 0
                n2 = N1 - (n + 1)
            else:
                n1 = n + 1
                n2 = N1

            for i1 in range(n1, n2):
                ndata1[i1, :] = data[i1 + s * n, :] + (
                    data[i1 + s * (n + 1), :] - data[i1 + s * n, :]
                ) * (dx / Dx)

            # dimension 2
            # make sure we stay within the defined array
            if s == 1:
                n1 = 0
                n2 = N2 - (n + 1)
            else:
                n1 = n + 1
                n2 = N2

            for i2 in range(n1, n2):
                ndata[:, i2] = ndata1[:, i2 + s * n] + (
                    ndata1[:, i2 + s * (n + 1)] - ndata1[:, i2 + s * n]
                ) * (dx / Dx)

        elif interpolation == "spline":
            pass

        elif interpolation == "fft":
            print("Shifting spectrum: fft interpolation")

            # inverse FFT
            ndata = numpy.fft.fftshift(self.data)
            ndata1 = numpy.fft.fft2(ndata)
            ndata = numpy.fft.fftshift(ndata1)

            timex = numpy.fft.fftfreq(self.xaxis.length, self.xaxis.step)
            timex = numpy.fft.fftshift(timex)
            timey = numpy.fft.fftfreq(self.yaxis.length, self.yaxis.step)
            timey = numpy.fft.fftshift(timey)

            # multiply by exponentials
            etx = numpy.exp(1j * dE * timex)
            ety = numpy.exp(1j * dE * timey)

            for k in range(ndata.shape[0]):
                ndata1[k, :] = etx[k] * ndata[k, :] * ety[:]

            # back FFT
            ndata = numpy.fft.ifft2(ndata1)
            ndata = numpy.fft.fftshift(ndata)

        else:
            raise QuantarheiError("Unknown interpolation type")

        self.data[:, :] = ndata[:, :]

    # FIXME: implement this
    def get_PumpProbeSpectrum(self) -> Any:  # type: ignore[explicit-any]
        """Returns a PumpProbeSpectrum corresponding to the 2D spectrum"""
        # from .pumpprobe import PumpProbeSpectrumCalculator
        from . import pumpprobe as pp

        # from ..core.time import TimeAxis
        # fake_t = TimeAxis(0,1,1.0)
        # ppc = PumpProbeSpectrumCalculator(fake_t, fake_t, fake_t)
        # return ppc.calculate_from_2D(self)
        return pp.calculate_from_2D(self)

    def plot(  # type: ignore[explicit-any]
        self,
        fig: Any = None,
        window: list | None = None,
        stype: str | None = None,
        spart: str = part_REAL,
        vmax: float | None = None,
        vmin_ratio: float = 0.5,
        colorbar: bool = True,
        colorbar_loc: str = "right",
        cmap: Any = None,
        plot_type: str = "both",
        Npos_contours: int = 10,
        positive_contour_color: Any = "k",
        negative_contour_color: Any = "k",
        zero_contour_color: Any = "b",
        zero_contour: bool = True,
        show_states: Any = None,
        text_loc: list | None = None,
        fontsize: str = "20",
        label: Any = None,
        show: bool = False,
        show_diagonal: Any = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        axis_label_font: Any = None,
    ) -> None:
        """Plots the 2D spectrum

        Parameters
        ----------
        fig : matplotlib.figure
            Figure into which plotting will be done. This is used e.g. when
            making a movie using moview writter (may be obsolete).
            If fig is None, we create a new figure object.

        window : list
            Specifies the plotted window in current energy units. When axes
            are x and y, the window is specified as window=[x_min,x_max,
            y_min,y_max]

        spart : {part_REAL, part_IMAGINARY, part_ABS}
            part of the spectrum to be plotted, constants such as part_REAL are
            defined at the highest import level of quantarhei.

        plot_type : {'both', 'image', 'contour'}
            Whether to plot the image with contour overlay ('both'), only the
            image ('image'), or only the contour lines ('contour').

        positive_contour_color : color
            Color used for positive contour levels.

        negative_contour_color : color
            Color used for negative contour levels.

        zero_contour_color : color
            Color used for the zero contour.

        zero_contour : bool
            If True, draw the zero contour line.

        vmax : float
            max of the plotting range in the z-direction. If vmax is None,
            maximum of the real part of the spectrum is used to determine
            the values of `vmax`




        """
        if text_loc is None:
            text_loc = [0.05, 0.9]

        if plot_type not in ("both", "image", "contour"):
            raise ValueError("plot_type must be 'both', 'image', or 'contour'")

        show_image = plot_type in ("both", "image")
        show_contours = plot_type in ("both", "contour")

        spect2D = self.data

        #
        # What part of the spectrum to plot
        #
        if spart == part_REAL:
            spect2D = numpy.real(spect2D)
        elif spart == part_IMAGINARY:
            spect2D = numpy.imag(spect2D)
        elif spart == part_ABS:
            spect2D = numpy.abs(spect2D)
        else:
            raise QuantarheiError("Undefined part of the spectrum: " + spart)

        if window is not None:
            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)

        else:
            i1_min = 0
            i1_max = self.xaxis.length
            i3_min = 0
            i3_max = self.yaxis.length

        #
        # Plotting with given units on axes
        #

        realout = spect2D[i3_min:i3_max, i1_min:i1_max]
        xdata = self.xaxis.data[i1_min:i1_max]
        ydata = self.yaxis.data[i3_min:i3_max]

        #
        #  How to treat the figures
        #
        if fig is None:
            fig, ax = plt.subplots(1, 1)
        else:
            fig.clear()
            ax = fig.add_subplot(1, 1, 1)

        #
        # Color map
        #
        if cmap is None:
            cmap = plt.cm.rainbow  # type: ignore[attr-defined]

        #
        # Actual plotting
        #
        if vmax is None:
            vmax = numpy.amax(realout)

        vmin = numpy.amin(realout)
        if vmin < -vmax * vmin_ratio:
            vmax = -vmin
        else:
            vmin = -vmax * vmin_ratio

        Npos = max(1, Npos_contours)
        poslevels = [i * vmax / Npos for i in range(1, Npos)]
        neglevels = [-i * vmax / Npos for i in range(Npos, 1, -1)]

        levo = float(xdata[0])
        prvo = float(xdata[-1])
        dole = float(ydata[0])
        hore = float(ydata[-1])

        cm = None
        if show_image:
            cm = ax.imshow(
                realout,
                extent=(levo, prvo, dole, hore),
                origin="lower",
                vmax=vmax,
                vmin=vmin,
                interpolation="bilinear",
                cmap=cmap,
            )

        #
        # Label
        #
        if label is not None:
            if isinstance(label, str) and len(text_loc) == 2:
                text_loc = [text_loc]
                label = [label]
                ln_t = 1
                ln_l = 1
            else:
                try:
                    ln_t = len(text_loc)
                    ln_l = len(label)
                except TypeError:
                    raise QuantarheiError(
                        "text_loc and label parameters must be"
                        " lists of the same lengths"
                    )

            if ln_t != ln_l:
                raise QuantarheiError(
                    "text_loc and label parameters have to have the"
                    " same number of members"
                )

            for pos, lbl in zip(text_loc, label):
                if lbl is not None:
                    ax.text(
                        (prvo - levo) * pos[0] + levo,
                        (hore - dole) * pos[1] + dole,
                        lbl,
                        fontsize=str(fontsize),
                    )

        #
        # Contours
        #
        if show_contours:
            warnings.filterwarnings("error")

            try:
                ax.contour(
                    xdata,
                    ydata,
                    realout,
                    levels=poslevels,
                    colors=positive_contour_color,
                )
            except (UserWarning, ValueError):
                pass

            if spart != part_ABS:
                if zero_contour:
                    try:
                        ax.contour(
                            xdata,
                            ydata,
                            realout,
                            levels=[0],
                            colors=zero_contour_color,
                        )
                    except (UserWarning, ValueError):
                        pass

                try:
                    ax.contour(
                        xdata,
                        ydata,
                        realout,
                        levels=neglevels,
                        colors=negative_contour_color,
                    )
                except (UserWarning, ValueError):
                    pass

            warnings.resetwarnings()

        #
        # Color bar presence
        #
        if colorbar and show_image and cm is not None:
            try:
                fig.colorbar(cm, ax=ax, location=colorbar_loc)
            except TypeError:
                fig.colorbar(cm, ax=ax)
        elif colorbar and not show_image:
            warnings.warn(
                "colorbar is ignored for contour-only plots",
                UserWarning,
            )

        #
        # Plot lines denoting positions of selected transitions
        #
        if show_states is not None:
            for en in show_states:
                try:
                    en1 = en[0]
                    co1 = en[1]
                except (TypeError, IndexError):
                    en1 = en
                    co1 = "--k"

                if en1 >= levo and en1 <= prvo:
                    ax.plot([en1, en1], [dole, hore], co1, linewidth=1.0)
                if en1 >= dole and en1 <= hore:
                    ax.plot([levo, prvo], [en1, en1], co1, linewidth=1.0)

        #
        # show diagonal line
        #
        if show_diagonal is not None:
            ax.plot([levo, prvo], [dole, hore], "-k", linewidth=1.0)

        #
        # axis labels
        #
        if axis_label_font is not None:
            font = axis_label_font
        else:
            font = {"size": 20}

        if xlabel is None:
            xl = ""
        if ylabel is None:
            yl = ""

        if xlabel is not None:
            xl = r"$\omega$ [fs$^{-1}$]"

        if isinstance(self.xaxis, FrequencyAxis):
            units = self.xaxis.unit_repr_latex()
            xl = r"$\omega_{1}$ [" + units + "]"
            yl = r"$\omega_{3}$ [" + units + "]"
        #        if isinstance(self.axis, TimeAxis):
        #            xl = r'$t$ [fs]'
        #            yl = r'$f(t)$'

        if xlabel is not None:
            xl = xlabel
        if ylabel is not None:
            yl = ylabel

        ax.set_xlabel(xl, **font)
        ax.set_ylabel(yl, **font)

        #
        # Should the spectra be showed now?
        #
        if show:
            self.show()

    def show(self) -> None:
        """Show the plot of 2D spectrum

        By default, plots are not shown. It is waited until explicit show()
        is called

        """
        plt.show()

    def savefig(self, filename: str) -> None:
        """Saves the fige of the plot into a file"""
        plt.savefig(filename, bbox_inches="tight")

    def trim_to(self, window: list | None = None) -> None:
        """Trims the 2D spectrum to a specified region

        Parameters
        ----------
        window : array
            Spectral window to which the present 2D spectrum
            should be trimmed. The window is specified as
            an array [w1_min, w1_max, w3_min, w3_max], where
            w1_min is the lower bound of the w1 axis (the x-axis),
            w1_max is the upper bound of the w1 axis and similarly
            for the w3 axis (y-axis) of the spectrum.


        """
        if window is not None:
            assert self.xaxis is not None
            assert self.yaxis is not None
            xaxis_fa: Any = self.xaxis  # type: ignore[explicit-any]
            yaxis_fa: Any = self.yaxis  # type: ignore[explicit-any]
            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)

            # create minimal off-set
            i1_min -= 1
            i1_max += 1
            i3_min -= 1
            i3_max += 1

            # reconstruct xaxis
            start_1 = self.xaxis.data[i1_min]
            length_1 = i1_max - i1_min
            step_1 = self.xaxis.step
            atype = xaxis_fa.atype

            xaxis = FrequencyAxis(
                start_1, length_1, step_1, atype=atype, time_start=xaxis_fa.time_start
            )
            self.xaxis = xaxis

            # reconstruct yaxis
            start_3 = self.yaxis.data[i3_min]
            length_3 = i3_max - i3_min
            step_3 = self.yaxis.step
            yaxis = FrequencyAxis(
                start_3,
                length_3,
                step_3,
                atype=yaxis_fa.atype,
                time_start=yaxis_fa.time_start,
            )
            self.yaxis = yaxis

            # reconstruct data
            if self.data is not None:
                data = self.data[i1_min:i1_max, i3_min:i3_max]
                self.data = data

        else:
            # some automatic trimming in the future
            pass
