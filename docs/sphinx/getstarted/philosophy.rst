Quantarhei Philosophy
=====================
Physics, Chemistry and Biology ultimately form a single broad subject.
With physical laws forming the basis of chemical laws and chemical and
physical laws forming the basis of biology, it seems that it the physicist
who is best positioned to answer all questions of fundamental importance.

This logic has one fundamental flaw, however. The problem is that there is
nothing like "a physicist" - physics itself became so broad that every
physicist is by her training a specialist in a limited sub-field of physics.
With every sub-field, there comes a different background knowledge, a different
set of standard mathematical techniques, and very often a different set
of approximation taken for granted. This should not be a surprise! The physics
of atoms and molecules is done on a completely different spatial scale than
astrophysics, but even when we stay in the realm of small, we see fundamental
differences between problems from, say, solid state and molecular physics.

These facts do not invalidate the uninity of physics as a science. They only
demonstrate practical human limits in attaining knowledge. It is very important
to keep in mind that every sub-field of physical science has accummulated 
its knowledge in only a lose contact with others subfields. Every sub-field
has its background knowledge which is often not very well documented in 
literature, as it is only passed down between generations by direct teaching.

A general definition of open quantum systems, which would satisfy experts
across a broad range of physical science, seems to be possible. Any system
which can be somehow defined as different from its environment, and which is
in contact with this environment, can be called open. It is also clear that
most, if not all, realistic systems are open. More on 
:ref:`open quantum systems <open-quantum-systems>`
in a separate section. However, when particular open quantum systems are
studied, it becomes clear, how different the descriptions can be. Techniques
and approximations suitable for descrition of transport phenomena in 
solid state systems are very different than in molecular physics. 

**Quantarhei aims at implementation and documentation of the shared know-how
of the molecular open quantum systems community, which meets over the problems
of the charge and energy transfer and spectroscopy of light-harvesting 
molecular aggregates**. This know-how can be very useful outside this
community, especially for studying other open quantum systems, but we are
aware of the limits of our endever.

Principles of Quantarhei
------------------------

Subject of Quantarhei
~~~~~~~~~~~~~~~~~~~~~

Acknowledging the diversity of sub-fields characterising open quantum
systems, the present software package - Quantarhei - has the word **molecular**
inserted in its subtitle: *molecular open quantum systems simulator*. At 
least for the foreseable future, it will concentrate only on the molecular
branch of open quantum systems theory. That is the first rule of Quantarhei
philosophy:

    1. :ref:`Quantarhei simulates molecular systems <mol-sys-label>`
    
In order to provide the user the best possible access to molecular world,
we treat every molecule in the system as an object. Molecular objects can be
assigned properties, groupped into aggregates and interactions between 
molecules can be assigned or calculated. The similarity between microscopic
objects like molecules and the concept of object in programming leads to the
second principle of Quantarhei philosophy:

    2. :ref:`Quantarhei is object oriented. Basic objects
       are molecules and their aggregates <object-oriented-label>`
    
Best examples of application of the open quantum systems theory implemented
in Quantarhei come from the study of natural light-harvesting. Chlorophyll
and Bacteriochlorophyll light-harvesting antennae of plants, algae and 
bacteria fit best the idea of an open quantum system as simulated by 
Quantarhei. 

Physics and its represetation in Quantarhei
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The molecules of these light-harvesting antennae interact strongly with each
other and they are embedded in protein 
environment. Interplay between the mutual inter-molecular interaction and
the interaction between the molecules and the protein *bath* significantly 
shapes the electronic properties of the whole system. In partictular, it 
depends on the competition between these two interactions, if the electronic
eigenstates of this molecular systems will be localized on individual
molecules, or delocalized over several molecules. In other words, eigenstates
of the system have different form, and we arrived at the problem of choosing
the most convenient representation of our quantum mechanical problem. However,
there is another principle at the basis of quantum mechanics, which makes it
(with a litle extension) into the list of Quantarhei founding principles:

    3. :ref:`Physics is basis independent, theories are not 
    <basis-in-quantarhei-label>`
    
This means that user should know in what basis she defines quantities,
she should be given a freedom to choose a basis to work in,
and she should be provided some means to control that choice. Quantarhei does 
some automatic handling of the basis for the user, and it allows her to
specify in which basis she defines, reads, saves, displays and plots physical
quantities. 

The second part of the Principle 3, *theories are not (basis independent)*
applies, for instance, to approximations. There are certain approximative 
theories which only make sense when applied in certain basis representation, 
because the conditions of validity of those theories imply certain basis 
representation as preferrable. Do you know `Fermi's golden rule`_? It does
not apply to arbitrary initial and final states, but only to those for
which the coupling can be assumed small. No arbitrary rotation of basis is 
allowed, when you want apply it. A good example from the field
of light-harvesting is the well-known Redfield
theory of relaxation. It can be formulated in a basis independent fashion,
but at some point you have to evaluate its quantities numericaly. In
eigenstate representation one obtains closed formulation of all terms of 
the tensor in terms of bath spectral density and transition frequencies. One
has to strictly apply this equation in eigenstate basis. It is very easy to
overlook that the final expression for the Redfield tensor elements does
not apply in arbitrary basis, but only in the one in which the Hamiltonian
of the system is diagonal. Additional approximations, such as the so-called
secular approximation also properly apply only in this basis. Quantarhei
does everything it can to ensure that it calculates relaxation tensors
in the proper representation. When such calculation is done, the resulting
relaxation tensor can be applied in any basis of your wish. Quantarhei
keeps track of *current basis of representation* and always transforms all
quantities for you to the current basis.

Somewhat similar situation as with the basis representation occurs in physics
with the representation of numbers. Very often, physical quantities have
dimensions and have to be expressed in certain units. There are very different
attitudes towards units between theorists and experimentalists. For a theorist,
it is just an additional nuisnance having to express her beautiful results,
obtained in *natural* units of the problem in some strange units used by
experimentalists, let alone in units standardized by some international
body. It is clear that physics does not care about our units, but people need
to pass the information among themselves using some predefined standards. 
Quantarhei recoginzes this by its fourth principle: 

    4. :ref:`Physics does not depend on units, numbers do <units-label>`
    
Quantarhei provides simple mechanisms for specifying and convering values
in and out of various unit systems. Again, defining, reading, saving and
displaying physical quantities can be done easily in any of the predefined
units. Internally, Quantarhei uses the most suitable units so that numerics
does not suffer. But as a user, you can specify and display the results in
the units most comfortable for you.


    5. :ref:`There is Quantum Mechanics deep down there <quantum-mechanics>`
    
    
\... to be written
    
Transparency of Quantarhei's implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This far we have been concerned with how physics is represented in Quantarhei.
It is important that Quantarhei does things right, but who is going to judge
the quality of its results? The only way to ensure the results Quantarhei
produces are correct is through the community control. Quantarhei is an open
source project so that everybody call check its internals and see how it
works inside.

    5. :ref:`Quantarhei is open source and transparent <transparency-label>`
    
All Quantarhei's functionality is available in open source format, with 
important numerical sections well separated from the more administrative
parts of the code. The separation should be such that one can check the
numerical routines independent of Quantarhei. They should only accept
simple arrays (matrices, arrays of real, complex and integer numbers),
numbers and strings as arguments. They should be well documented so that one
can check the correctness of implementation. To ensure that Quantarhei is
completely transparent, we decided to have it written in Python using only
standard scientific libraries such as numpy and scipy. The principle is

    6. :ref:`Quantarhei is written in Python <why-python-label>`


On the face of it, both the 5th and 6th principles would be too restrictive
if applied strictly. Python, even with numpy, might not provide the best
performace for many numerical problems. Being completely open source might
also not always be the best practise. What if, in the future, a non open
source library will be available for free (e.g. for non-profit community)
with some of the desired functionality, and with huge performace benefits. It 
would be perhaps desirable to integrate it with Quantarhei. However, in order
to keep the promis of transparency, it is required that an alternative
(perhaps low performace) implementation of the functionality exists in 
Quantarhei for both reasons of testing the functionality and perhaps
pedagogical reasons. This requires too things, first non-Python extensions
must be possible, and a mechanism for keeping the transparency must be
provided. The seveth Quantarhei principle is

    7. :ref:`Quantarhei is extensible <extensions-label>`

Quantarhei provides a standardized mechanism for writting extensions.
Different implementations should be allowed to *compete* for the same 
functionality, i.e. the user should be able to select an implementation, if
more than one are available. At the same time, we want that a simple
Python implementation is always available so that people could install a 
Python-only source distribution and play with it and its source code. 
Therefore:
    
    8. :ref:`Non-Python extensions must be coverred by a Python implementation <non-python-label>`
    
Quantarhei does not want to prohibit its commercial use. If Quantarhei becomes
a platform which can be used for industrial simulations, and provides would
be motivated to sell high performace modules for Quantarhei, we (the 
originators of this project) would be extremely happy about it.    

\... to be continued

Reproducibility of Simulations with Quantarhei
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    9. Simulation by Quantarhei are reproducible
    
    
Teaching with Quantarhei
~~~~~~~~~~~~~~~~~~~~~~~~

    10. Quantarhei is well documented with examples and templates available
    
    
\... to be continued
    
.. 
    |Qrhei|_ is an object oriented library for writting Python scripts to simulate
    problems related to excited state dynamics and spectroscopy of molecular 
    systems. 

Some More Details of the Quantarhei Principles
----------------------------------------------

.. toctree::
    :maxdepth: 1
    
    details/molsys
    details/object
    details/basis
    details/units
    details/qm
    details/transparency
    details/python
    details/extensions
    details/non_python
    

.. |Qrhei| replace:: **Quantarhei**
.. _Qrhei: http://github.com/tmancal74/quantarhei

.. _`Fermi's golden rule`: https://en.wikipedia.org/wiki/Fermi%27s_golden_rule