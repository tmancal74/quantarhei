from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy

from ...core.time import TimeAxis


def _normalize_integrals(integrals: Any) -> list[str]:
    """Return a stable list of integral-time labels from config metadata."""
    if integrals is None:
        return []
    if isinstance(integrals, str):
        return [integrals]
    if isinstance(integrals, set):
        return sorted(integrals)
    return list(integrals)


def _label_from_terms(terms: list[tuple[int, str]]) -> str:
    """Build a compact label such as ``t1+t2-s2`` from signed terms."""
    label = ""
    for kk, (sign, name) in enumerate(terms):
        if sign not in (-1, 1):
            raise ValueError("Only unit signs are supported in time labels")
        if kk == 0:
            if sign < 0:
                label += "-"
        else:
            label += "+" if sign > 0 else "-"
        label += name
    return label


def _parse_time_label(label: str) -> list[tuple[int, str]]:
    """Parse labels made of variable names joined by ``+`` and ``-`` signs."""
    terms = []
    sign = 1
    token = ""
    for char in label:
        if char in "+-":
            if token:
                terms.append((sign, token))
                token = ""
            sign = 1 if char == "+" else -1
        elif not char.isspace():
            token += char
    if token:
        terms.append((sign, token))
    if not terms:
        raise ValueError("Empty time label")
    return terms


def _generate_time_argument_labels(response_times: dict[str, Any]) -> list[str]:
    """Generate all supported line-shape-function time arguments.

    The ordinary arguments are all non-empty ordered combinations of response
    times. Integral variables attached to reset times generate two additional
    variants for every argument containing their parent: one with the parent
    replaced by the integration time, and one with it replaced by the remaining
    interval ``parent-integral``.
    """
    labels = []
    names = list(response_times)

    reset_integrals = {}
    for name, info in response_times.items():
        if info.get("reset", False):
            integrals = _normalize_integrals(info.get("integral", None))
            if integrals:
                reset_integrals[name] = integrals

    for length in range(1, len(names) + 1):
        for combo in combinations(names, length):
            labels.append("+".join(combo))

            replacement_groups = []
            for name in combo:
                variants = [[(1, name)]]
                for integral in reset_integrals.get(name, []):
                    variants.append([(1, integral)])
                    variants.append([(1, name), (-1, integral)])
                replacement_groups.append(variants)

            generated: list[list[tuple[int, str]]] = [[]]
            for variants in replacement_groups:
                generated = [
                    base + variant for base in generated for variant in variants
                ]

            for terms in generated[1:]:
                labels.append(_label_from_terms(terms))

    unique_labels = []
    for label in labels:
        if label not in unique_labels:
            unique_labels.append(label)
    return unique_labels


def _default_storage_config(one_jump: bool = False) -> dict[str, Any]:
    """Return the standard third-order response storage configuration."""
    t2_integral = {"s2"} if one_jump else None
    return {
        "response_times": {
            "t1": {"reset": False, "integral": None, "axis": 0},
            "t2": {"reset": True, "integral": t2_integral},
            "t3": {"reset": False, "integral": None, "axis": 1},
        }
    }


class FunctionStorage:
    """Data storage for discrete representation of correlation functions
    and lineshape functions.


    Parameters
    ----------
    N: int
        Number of functions to be stored.


    """

    import numpy

    def __init__(
        self,
        N: int = 1,
        timeaxis: Any = None,
        dtype: Any = numpy.complex64,
        show_config: bool = False,
        config: dict | int | None = None,
    ) -> None:

        #
        # Storage configuration
        #
        # This dictionary describes what g(t) values will be stored and how
        #
        self.config: dict[str, Any]
        if config is None or config == 0:
            self.config = _default_storage_config(one_jump=False)
        elif config == 1:
            self.config = _default_storage_config(one_jump=True)
        elif isinstance(config, dict):
            self.config = config
        else:
            raise ValueError("FunctionStorage config has to be None, 0, 1, or a dict")

        #######################################################################
        #
        # Analyze the configuration
        #
        #######################################################################

        #
        # generate the list of argument combinations
        #
        values = self.config["response_times"]
        if self.config.get("arguments", "auto") == "auto":
            self.time_combination_strings = _generate_time_argument_labels(values)
        else:
            self.time_combination_strings = list(self.config["arguments"])
        self.time_combination_tuples = [
            tuple(name for _, name in _parse_time_label(label))
            for label in self.time_combination_strings
        ]

        if show_config:
            print("\nStorage configuration")
            print("Time argument combinations: ", self.time_combination_strings)

        #
        # creating time_index
        #
        time_index: dict[str, int | tuple[int, ...]] = {}
        variable_index: dict[str, int] = {}

        # find all reset times
        reset_times = []
        self.integral_parent = {}
        for rst_key in self.config["response_times"]:
            if self.config["response_times"][rst_key]["reset"]:
                reset_times.append(rst_key)
                variable_index[rst_key] = -1
                for integral in _normalize_integrals(
                    self.config["response_times"][rst_key].get("integral", None)
                ):
                    variable_index[integral] = -1
                    self.integral_parent[integral] = rst_key

        # find all 1D times
        oned_times = []
        for rst_key in self.config["response_times"]:
            if not self.config["response_times"][rst_key]["reset"]:
                axis = self.config["response_times"][rst_key]["axis"]
                oned_times.append(rst_key)
                variable_index[rst_key] = axis

        # is the number of 1D times matching the number of submitted time axes?
        if isinstance(timeaxis, TimeAxis):
            tlen = 1
        else:
            tlen = len(timeaxis)
        if len(oned_times) != tlen:
            print("Storage configuration:")
            print(self.config)
            raise Exception(
                "The number of time axes has to be consistent with storage conguration."
            )

        self.oned_times = oned_times

        self.variable_index = variable_index
        self.time_expressions = {}
        self.time_dependencies: dict[str, set[str]] = {}

        # find all time arguments and calculate their dimension
        for tstr in self.time_combination_strings:
            terms = _parse_time_label(tstr)
            axes = []
            dependencies: set[str] = set()
            for _, tm in terms:
                if tm not in variable_index:
                    raise Exception("Unknown time variable in storage config: " + tm)
                dependencies.add(tm)
                axis = variable_index[tm]
                if axis >= 0 and axis not in axes:
                    axes.append(axis)

            self.time_expressions[tstr] = terms
            self.time_dependencies[tstr] = dependencies
            if len(axes) == 0:
                time_index[tstr] = -1
            elif len(axes) == 1:
                time_index[tstr] = axes[0]
            else:
                time_index[tstr] = tuple(axes)

        # order time_index
        # Sorting criteria:
        # - First, sort by type: Integers should come before tuples
        # - Second, for integers: Sort numerically
        # - Third, for tuples: Sort by length
        sorted_items = sorted(
            time_index.items(),
            key=lambda x: (
                isinstance(x[1], tuple),
                x[1] if isinstance(x[1], int) else len(x[1]),
            ),
        )

        # Convert back to a dictionary (optional)
        time_index = dict(sorted_items)

        if show_config:
            print("Property 'time_index': ", time_index)
            # print(self.time_combination_tuples)

        ##########################################################################
        #
        ##########################################################################

        # Number of unique functions stored
        self.N = N

        # The max value of the index identifying the stored functions for the user
        self.Nmax = self.N - 1

        # One or more time axes on which the functions should be represented
        if timeaxis is None:
            raise Exception("At least one time axis has to be specified")
        self.timeaxis = timeaxis

        # Storage dimensions; they correspond to the submitted time axes
        if isinstance(timeaxis, TimeAxis):
            self.Ndim = 1
            dim_arr: numpy.ndarray = numpy.zeros(self.Ndim, dtype=numpy.int32)
            dim_arr[0] = timeaxis.length
        else:
            self.Ndim = len(timeaxis)
            dim_arr = numpy.zeros(self.Ndim, dtype=numpy.int32)
            for ii, ta in enumerate(timeaxis):
                if not isinstance(ta, TimeAxis):
                    raise Exception(
                        "'timeaxis' argument must be a list/tuple of Quantarhei TimeAxis objects."
                    )
                dim_arr[ii] = ta.length
        self.dim = dim_arr

        # Data type of the storage (default is complex64)
        self.dtype = dtype

        # Dictionary of time arguments and indices of the 'dim' argument defining the corresponding
        # length of the discrete representation. The mapping of the time arguments on the storage integer
        # indices is done in a natural order.
        #
        # The default setting means that the value stored under the index 0 is the one of g(t2) and
        # its representation has a size of 1 element (this is the meaning of -1). The value stored
        # under index 6 is the one of g(t1+t2+t3) with a length of dim[0]*dim[1] (this what the values
        # in the tuple define).
        # self.time_index = {"t2":-1,"t1":0,"t1+t2":0, "t3":1,"t2+t3":1,"t1+t3":(0,1),"t1+t2+t3":(0,1)}
        self.time_index = time_index

        reshapes: dict[str, list[Any]] = dict()
        for dms in self.time_index:
            val = self.time_index[dms]
            if isinstance(val, int):
                if val == -1:
                    reshapes[dms] = []
                else:
                    reshapes[dms] = [self.dim[val]]

            else:
                rlist = []
                dcount = 0
                kl = 0
                for kk, dims in enumerate(val):
                    if dims >= 0:
                        rlist.append(self.dim[dims])
                        dcount += 1
                        kl += 1

                reshapes[dms] = rlist

        self.reshapes_dic = reshapes
        self.reshapes = []
        for key in self.reshapes_dic:
            self.reshapes.append(self.reshapes_dic[key])

        # The number of time arguments
        self.Nt = len(self.time_index)

        # mapping of times to indices
        self.time_mapping = {}
        kk = 0
        for tm in self.time_index:
            self.time_mapping[tm] = kk
            kk += 1

        #
        # Extract the sizes from 'dim' argument
        #
        # self._N stores the sizes of individual index representations
        self._N = numpy.zeros(self.Nt, dtype=numpy.int32)

        # loop over time indices
        kk = 0
        for ts in self.time_index:
            ival = self.time_index[ts]
            if isinstance(ival, tuple):
                self._N[kk] = 1
                for sz in ival:
                    if sz >= 0:
                        self._N[kk] *= dim_arr[sz]
            elif isinstance(ival, int):
                if ival == -1:
                    self._N[kk] = 1
                else:
                    self._N[kk] = dim_arr[ival]
            kk += 1

        # one function takes one dimensional array of length 'data_stride'
        self.data_stride = 0
        for kk in range(self.Nt):
            self.data_stride += self._N[kk]

        # total number of values
        self.data_dim = self.N * self.data_stride

        # data size in MB
        # print(self.data_dim)
        # print(numpy.dtype(self.dtype))
        self.data_size = (
            self.data_dim * numpy.dtype(self.dtype).itemsize / (1024 * 1024)
        )
        self.size_units = "MB"

        #
        # Here the data are stored
        #
        self.data = numpy.zeros(self.data_dim, dtype=self.dtype)

        # 2D view of the data for fast access
        self._data2d = self.data.reshape((self.N, self.data_stride))
        self._last_reset_values: dict[str, Any] | None = None
        self._data_initialized = False

        #
        # Here the g(t) functions are stored
        #
        self.funcs: dict[Any, Any] = {}

        # The number of stored functions
        self.Nf = 0

        #
        # These arrays store the starts and ends of the g(t) arrays in a give stride
        #
        self.start = numpy.zeros(self.Nt, dtype=numpy.int32)
        self.end = numpy.zeros(self.Nt, dtype=numpy.int32)

        # filling the starts
        for kk in range(self.Nt):
            if kk == 0:
                self.start[kk] = 0
            else:
                self.start[kk] = self.start[kk - 1] + self._N[kk - 1]

        # filling the ends
        for ii in range(self.Nt):
            self.end[ii] = self.start[ii] + self._N[ii]

        #
        # As an alternative, one can use the string to access the data
        #
        self.start_dic = {}
        self.end_dic = {}
        kk = 0
        for tm in self.time_mapping:
            self.start_dic[tm] = self.start[kk]
            self.end_dic[tm] = self.end[kk]
            kk += 1

        #
        # Mapping between the actual stored functions and the index representing
        # the functions for the outside world.
        #
        self.mapping = numpy.zeros(self.Nmax + 1, dtype=numpy.int32)
        for kk in range(self.N):
            self.mapping[kk] = -1

        if show_config:
            print("\n")

    def __str__(self) -> str:
        """String representation of the storage"""
        return (
            "Correlation function storage\n"
            " dimension: "
            + str(self.dim)
            + "\n"
            + " data size: "
            + str(self.data_size)
            + " "
            + self.size_units
        )

    def __setitem__(self, index: Any, value: Any) -> None:
        """Set the value(s) of the storage"""
        if isinstance(index, tuple):
            if len(index) == 2:
                i, j = index
                start = i * self.data_stride

                if isinstance(j, int):
                    sta = start + self.start[j]
                    end = start + self.end[j]

                else:
                    sta = start + self.start_dic[j]
                    end = start + self.end_dic[j]

                self.data[sta:end] = numpy.array(value, dtype=self.dtype)

            else:
                raise Exception()

        else:
            raise Exception()

    def __getitem__(self, index: Any) -> Any:
        """Returns the stored function"""
        if isinstance(index, tuple):
            if len(index) == 2:
                i, j = index

                # if i is an integer, we return one array
                if isinstance(i, int):
                    start = self.mapping[i] * self.data_stride

                    if isinstance(j, int):
                        sta = start + self.start[j]
                        end = start + self.end[j]

                    else:
                        sta = start + self.start_dic[j]
                        end = start + self.end_dic[j]

                    if isinstance(j, int):
                        tpl = self.reshapes[j]
                    else:
                        tpl = self.reshapes_dic[j]

                    if sta + 1 == end and len(tpl) == 0:
                        return self.data[sta]

                    return self.data[sta:end].reshape(tpl)

                # if i is a slice :, we return a view an all arrays
                if isinstance(i, slice) and i == slice(None):
                    if isinstance(j, int):
                        sta = self.start[j]
                        end = self.end[j]

                    else:
                        sta = self.start_dic[j]
                        end = self.end_dic[j]

                    if isinstance(j, int):
                        tpl = [self._data2d.shape[0]] + self.reshapes[j]  # type: ignore[misc]
                    else:
                        tpl = [self._data2d.shape[0]] + self.reshapes_dic[j]  # type: ignore[misc]

                    if sta + 1 == end and len(tpl) == 1:
                        return self._data2d[i, sta]

                    return self._data2d[i, sta:end].reshape(tpl)

            else:
                raise Exception()

        else:
            raise Exception()

    def set_goft(self, N: Any, func: Any = None) -> None:
        """Sets the values for a given stored function."""
        #
        # The argument N must be an integer or a list of indices by which the function should
        # represented for the outside world
        #
        if not isinstance(N, (int, list, tuple, numpy.ndarray)):
            raise Exception(
                "Argument N has to be an integer or an array (list, tuple) of integers"
            )

        #
        # Check the consistency of the submitted time axes
        #
        self.ta: list[Any] = []  # This is a synonym for the self.timeaxis
        if isinstance(self.timeaxis, (tuple, list)):
            self.ta = list(self.timeaxis)
            # self.t1a = self.timeaxis[0]
            # self.t3a = self.timeaxis[1]
        else:
            self.ta = [self.timeaxis]
            # self.t1a = self.timeaxis
            # self.t3a = self.timeaxis

        # for an integer
        if isinstance(N, int):
            self._check_and_make_space(N)

            # FIXME: mapping has to be merged with the previous mappings
            if self.mapping[N] > 0:
                raise Exception("Function already set for " + str(N))
            # set the function to the first available position
            if self.Nf > self.N:
                raise Exception("The storage full. Unable to add additional functions")

            if func is not None:
                pos = self._fce_already_stored(func)
                if pos == -1:
                    self.funcs[self.Nf] = func
                    self.mapping[N] = self.Nf
                    self.Nf += 1
                else:
                    self.mapping[N] = pos

            else:
                raise Exception("Function not submitted")

        else:
            # In this case N is a list of mappings to the submitted function
            added = False
            kk = 0
            for nn in N:
                self._check_and_make_space(nn)

                # if the position is free
                if self.mapping[nn] < 0:
                    # check that the same function was not submitted before
                    pos = self._fce_already_stored(func)
                    if pos == -1:
                        if kk == 0:
                            self.funcs[self.Nf] = func
                            added = True
                        self.mapping[nn] = self.Nf
                        kk += 1
                    else:
                        self.mapping[nn] = pos

                else:
                    raise Exception("Function already set for " + str(nn))

            # we added just one function
            if added:
                self.Nf += 1

        self._data_initialized = False

    def _check_and_make_space(self, nn: int) -> None:
        """Check if the required position is outside the allocated mapping
        and if so, make more space.

        """
        if nn > self.Nmax:
            # create more space
            mapping = numpy.zeros(nn + 1, dtype=numpy.int32)
            mapping[:] = -1
            # copy earlier values
            # print(nn, self.Nmax, mapping.shape, self.mapping.shape)
            mapping[: self.Nmax + 1] = self.mapping
            # set new Nmax
            self.Nmax = nn
            # make larger variables the object properties
            self.mapping = mapping

    def _fce_already_stored(self, func: Any) -> int:
        """If the function was already stored in the object, the function returns its position,
        otherwise it returns -1.

        """
        ret = -1
        if len(self.funcs) == 0:
            return ret

        for key in self.funcs:
            fc = self.funcs[key]
            if fc is func:
                ret = key

        return ret

    def create_data(self, reset: dict | None = None, force: bool = False) -> None:
        """We create data for all submitted functions"""
        if reset is None:
            reset = {"t2": 0.0}
        reset = dict(reset)
        reset_defaults = getattr(self, "_reset_defaults", {})
        for integral, parent in self.integral_parent.items():
            if integral not in reset:
                if integral in reset_defaults:
                    reset[integral] = reset_defaults[integral]
                elif parent not in reset:
                    raise Exception("Parent reset time not specified: " + parent)
                else:
                    reset[integral] = reset[parent] / 2.0

        changed_variables = self._changed_reset_variables(reset)
        if force or not self._data_initialized or self._last_reset_values is None:
            labels_to_update = list(self.time_index)
        else:
            labels_to_update = [
                label
                for label in self.time_index
                if self.time_dependencies[label].intersection(changed_variables)
            ]

        if not labels_to_update:
            self._last_reset_values = dict(reset)
            return

        tt_matrix = []

        _colon_ = slice(None, None, None)

        for tm in labels_to_update:
            dms = self.time_index[tm]
            if isinstance(dms, int):
                axes = [] if dms == -1 else [dms]
            else:
                axes = [axis for axis in dms if axis >= 0]

            tt_dim = [self.dim[axis] for axis in axes]
            if len(tt_dim) == 0:
                tt = numpy.zeros(1, dtype=numpy.float64)
            else:
                tt = numpy.zeros(tt_dim, dtype=numpy.float64)

            for sign, name in self.time_expressions[tm]:
                axis = self.variable_index[name]
                if axis == -1:
                    if name not in reset:
                        raise Exception("Reset time not specified: " + name)
                    tt += sign * reset[name]
                elif len(tt_dim) == 1:
                    tt += sign * self.ta[axis].data
                else:
                    index_list: list[slice | None] = [None] * len(tt_dim)
                    index_list[axes.index(axis)] = _colon_
                    tt += sign * self.ta[axis].data.__getitem__(tuple(index_list))

            tt_matrix.append(tt)

        #
        # loop over all stored functions
        #
        for key in self.funcs:
            func = self.funcs[key]
            # start = key*self.data_stride

            for kk, tm in enumerate(labels_to_update):
                tm_index = self.time_mapping[tm]
                self.__setitem__(
                    (key, tm), func(tt_matrix[kk]).reshape(self._N[tm_index])
                )

        self._last_reset_values = dict(reset)
        self._data_initialized = True

    def _changed_reset_variables(self, reset: dict) -> set[str]:
        """Return reset variables whose values differ from the previous call."""
        if self._last_reset_values is None:
            return set(reset)

        changed: set[str] = set()
        for name, value in reset.items():
            if name not in self._last_reset_values:
                changed.add(name)
            elif not numpy.array_equal(value, self._last_reset_values[name]):
                changed.add(name)

        for name in self._last_reset_values:
            if name not in reset:
                changed.add(name)

        return changed

    def effective_size(self) -> int:
        """Effective size of the storage. Some stored functions have the same functional form"""
        return len(self.mapping)

    def get_number_of_functions(self) -> int:
        """Returns the number of different stored functions"""
        return int(numpy.max(self.mapping)) + 1

    def get_number_of_sites(self) -> int:
        """Returns the number of assigned sites"""
        return (self.mapping >= 0).sum()

    def get_mapping_matrix(self) -> numpy.ndarray:
        """Returns a matrix that maps site index on the function index


        M[n,i] equals 1 if the site n has the correlation function stored
        at the position i in the storage

        Returns
        -------
        Mapping matrix


        """
        Nstore = self.get_number_of_functions()
        Nsites = self.get_number_of_sites()

        if Nsites != len(self.mapping):
            raise Exception(
                "Mapping has to be complete."
                " All of its elements have to be set without gaps."
            )

        mtrx = numpy.zeros((Nsites, Nstore), dtype=numpy.int32)
        for ii in range(Nsites):
            jj_ind = self.mapping[ii]
            mtrx[ii, jj_ind] = 1.0

        return mtrx

    def get_reorganization_energies(self) -> numpy.ndarray:
        """Returns the estimate of the reoganization energies of the stored functions"""
        label = "t1"
        axis = self.ta[0]
        for candidate in self.oned_times:
            axis_index = self.variable_index[candidate]
            candidate_axis = self.ta[axis_index]
            if not numpy.isclose(candidate_axis.data[-1], 0.0):
                label = candidate
                axis = candidate_axis
                break

        index = (slice(None, None, None), label)
        igg = -numpy.imag(self.__getitem__(index))
        tal = axis.length - 1
        if numpy.isclose(axis.data[tal], 0.0):
            return numpy.zeros(igg.shape[0], dtype=igg.dtype)
        lam = igg[:, tal] / axis.data[tal]

        return lam


class FastFunctionStorage(FunctionStorage):
    """Function storage with minimal-overhead array retrieval.

    A subclass of :class:`FunctionStorage` optimised for fast random access
    via direct array slicing.  Supports integer or full-slice ``:``) as the
    first index.

    Parameters
    ----------
    N : int, optional
        Number of functions to store. Default is ``1``.
    timeaxis : TimeAxis or list of TimeAxis
        One or more time axes on which functions are discretised.
    dtype : numpy.dtype, optional
        Data type of the stored arrays. Default is ``numpy.complex64``.
    show_config : bool, optional
        If ``True``, print the storage configuration at construction time.
        Default is ``False``.
    config : dict or None, optional
        Custom storage configuration dictionary. If ``None``, the default
        three-time (t1, t2, t3) configuration is used.
    """

    def __getitem__(self, index: Any) -> Any:
        i, j = index

        # if the index is integer
        if isinstance(i, int):
            start = self.mapping[i] * self.data_stride
            sta = start + self.start[j]
            end = start + self.end[j]
            if sta + 1 == end:
                return self.data[sta]
            return self.data[sta:end]

        # i is a slice :
        if isinstance(i, slice) and i == slice(None):
            sta = self.start[j]
            end = self.end[j]
            if sta + 1 == end:
                return self._data2d[i, sta]
            return self._data2d[i, sta:end]

        raise Exception("Integer or slice : required as a first index")


class SingleGoft(FunctionStorage):
    """Storage of a single lineshape function"""

    def __init__(
        self,
        timeaxis: Any = None,
        dtype: Any = numpy.complex64,
        show_config: bool = False,
        config: dict | None = None,
    ) -> None:
        super().__init__(
            N=1, timeaxis=timeaxis, dtype=dtype, show_config=show_config, config=config
        )

    def __getitem__(self, index: Any) -> Any:
        """Returns the stored function"""
        j = index
        start = 0

        if isinstance(j, int):
            sta = start + self.start[j]
            end = start + self.end[j]

        else:
            sta = start + self.start_dic[j]
            end = start + self.end_dic[j]

        return self.data[sta:end]
