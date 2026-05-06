"""This is the class representing tightly organized molecular aggregates such
as photosynthetic antenna and light-harvesting complexes in
Quantarhei. Appart from represeting data, this class also provides a
simplified interface to much of Quantarhei's functionality, such as
calculation of spectra and dynamics. In order to make the core more
organized, the class `Aggregate` is the tip of series of mutually
inheriting classes. They start with AggregateBase, a class which implements
some of the core functionality and add functionality in classes like
`AggregateSpectroscopy`, `AggregateExcitonAnalysis` etc.

Inheritance in Aggregate class
------------------------------

The dependency of the classes is the following

AggregateBase :
basic functionality of the Aggregate

AggregateSpectroscopy :
adds Liouville pathway generation

AggregateExcitonAnalysis :
adds analysis of excitons

AggregatePureDephasing :
adds calculation of effective pure dephasing rates

Aggregate :
wraps everything up


Class Details
-------------

"""

from .aggregate_pdeph import AggregatePureDephasing


class Aggregate(AggregatePureDephasing):
    """Tightly organized molecular aggregate such as a photosynthetic antenna.

    Combines a collection of :class:`~quantarhei.builders.molecules.Molecule`
    objects into a coupled system and provides a unified interface for
    building Hamiltonians, computing spectra, and simulating dynamics.
    This class is the final layer in a hierarchy of inheriting classes that
    starts with ``AggregateBase``.

    Notes
    -----
    Call :meth:`build` after adding all molecules and setting all parameters
    before using any calculation methods.
    """

    pass
