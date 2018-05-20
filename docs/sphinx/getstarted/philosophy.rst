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

Principles of Quantarhei
------------------------

Acknowledging this diversity of sub-fields characterising open quantum
systems, the present software package - Quantarhei - has the word **molecular**
inserted in its subtitle: *molecular open quantum systems simulator*. At 
least for the foreseable future, it will concentrate only on the molecular
branch of open quantum systems theory. That is the first rule of Quantarhei
philosophy:

    1. **Quantarhei simulates molecular systems**
    
In order to provide the user the best possible access to molecular world,
we treat every molecule in the system as an object. Molecular objects can be
assigned properties, groupped into aggregates and interactions between 
molecules can be assigned or calculated. The similarity between microscopic
objects like molecules and the concept of object in programming leads to the
second principle of Quantarhei philosophy:

    2. **Quantarhei is object oriented. Basic objects
       are molecules and their aggregates**
    
Best examples of application of the open quantum systems theory implemented
in Quantarhei come from the study of natural light-harvesting. Chlorophyll
and Bacteriochlorophyll light-harvesting antennae of plants, algae and 
bacteria fit best the idea of an open quantum system as simulated by 
Quantarhei. 

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

    3. **Physics is basis independent, people are not**
    
This means that user should know in what basis she defines quantities,
she should be given a freedom to choose a basis to work with,
and provided some means to control that choice. Quantarhei does some automatic
handling of the basis for the user, and it allows her to specify in which 
basis she defines, displays and plots physical quantities.

\... to be continued
    
.. 
    Molecular and aggregate objects in Quantarhei form the layer closest to the
    user. Below this layer there is a layer describing quantum mechanics of the
    problem. In quantum mechanics of open system we deal with operators and 
    superoperators, which are very happy about being treated as objects. 
    Especially, the fact that an abstract object is independent of its mathematical
    representation (the basis of representation) is well represented by the 
    idea of object in programming languages. 


    |Qrhei|_ is an object oriented library for writting Python scripts to simulate
    problems related to excited state dynamics and spectroscopy of molecular 
    systems. 


.. |Qrhei| replace:: **Quantarhei**
.. _Qrhei: http://github.com/tmancal74/quantarhei

