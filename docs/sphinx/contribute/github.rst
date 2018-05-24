Using the Github Repository
===========================

Quantarhei is freely available in source code form from a Git repository at
github.com_. The source code is managed by the Git, a version control system
originated by Linus Torvalds. 

Whether you plan to actively develop Quantarhie, or just modify the source
code for yourself, the most convenient way to handle the source code
is through the `git` command from your terminal. If you don't like terminal,
there are plenty of Git based graphical managers, but they will not be described
here. To translate the commands discussed here to mouse clicks should be
way easier than the oposite.


Forking and cloning the repository
----------------------------------

While you can always just download the code from `github.com`_, it is way more
useful to **clone** the online repository, i.e. to create your local copy 
of the source code, which is (just as the original repository) aware of its
history and enables you to record and manage changes to the source code.
Ultimitely, such a cloned repository can be used to upload your changes and
additions to the original repository. If you mean contributing to the project
seriously, you should do one more step before you clone the repository. You
shoud **fork** it. Forking the repository on Github creates your own on-line
version of the repository to which you will push your changes and additions,
so that they do not stay localized on your computer only. The advantage on
having such a forked repository on Github is that you can do the so-called
**pull requests** using the Github website. Pull request is a process in which
your code is revised and eventually *pulled* into the original repository
from which you forked. This is the way you contribute to our project.

Consult `Github website`_ about getting an account. When you have an account
search for Quantarhei and **fork it** using a Fork button on top-right side
of the `project Github website`_. You will have your own, independent version
of Quantarhei on-line.

To proceed further you will need the **Git** command line tool installed on
your computer.
Consult `Git`_ website for installation instructions and download for your
platform. When the Git is installed on your computer, create a suitable
directory, where you will keep your Git repositories, and clone Quantarhei.
It is done by the following command

.. code:: bash

    $ git clone https://github.com/YOUR-GITHUB_NAME/quantarhei.git 

This will create your private copy of Quantarhei package in a directory called
`quantarhei`. 

It is your private copy, it is connected to your forked on-line repository on
Github, and you can do to it any changes you want, without disturbing anybody
and anything.

Building and installing from the cloned repository
--------------------------------------------------

One thing which you can do with your new repository is to build and install
Quantarhei. The build part is done using Python and a script which can be 
found the package root directory. First you enter the root directory of the 
package

.. code:: bash

    $ cd quantarhei/
    
and then you run the building process by typing the following:

.. code:: bash

    $ python setup.py sdist
    
This will create an installable package, stored in the `dist` directory inside
the package directory (where we are located now).

Installation of the package (let's assume it is called
`quantarhei-0.0.36.tar.gz`) can be done easily by using the `pip`
command

.. code:: bash

    $ pip install dist/quantarhei-0.0.36.tar.gz
    
The package is now install on your system.



Keeping track of remote changes
-------------------------------

If you plan to make changes to the sotware you just downloaded (by cloning the
repository), you would like to receive update as a new work is added to the 
original repository. This can be easily done.

\... to be continued


.. _`github.com`: https://github.com/tmancal74/quantarhei
.. _`Git`: https://git-scm.com
.. _`Github website`: http://github.com
.. _`project Github website`: https://github.com/tmancal74/quantarhei
