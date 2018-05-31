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

You can now test the package by calling the local `Makefile`. Before you do
that, install some packages on which the testing depends. The easiest way is
to you the `requirements_test.txt` file in the root directory of the package:

.. code:: bash

    $ pip install -r requirements_test.txt
    
This will install all required packages. 

The next step is to run the tests: 

.. code:: bash

    $ make test
    
It will take couple of minutes to finish, and it should end up with an `OK`.

Configuring Git
---------------

In order for your repository to recognize you when you push in your changes,
you need to configure the information Git has about you. Set your user name,
and set your email (the one that is associated with your GitHub account)
by the following commands

.. code:: bash

    $ git config --global user.name "Dobran Krasa"
    $ git config --global user.email dobrank@gmail.com

This should be enough to keep you well identified, when you push your changes
to the on-line repository.

\... to be continued

Keeping track of remote changes
-------------------------------

If you plan to make changes to the sotware you just downloaded (by cloning the
repository), you would like to receive updates as new work is added to the 
original repository. This can be easily done. In the `Makefile` which is
included in the root directory of the package, there is a predefined 
command which connects your local repository to the original on-line repository
and allows you to do regular updates.

First time you want to update, you need to type:

.. code:: bash

    $ make git_add_upstream
    
This will set the original on-line Quantarhei repository, where all the
development converges as an `remote upstream` repository for your project.
You can now do updates regularly using another task predefined in the Makefile:

.. code:: bash

    $ make git_update_master
    
This command merges all the changes from the master branch of the Quantarhei
project into your master branch. It leaves you in the *master* branch. If you
need to update the brach *my_branch* on which you work, do the following (it
is assumed you do it right after the previous command):

.. code:: bash

    $ git checkout my_branch
    $ git merge master
    
Now you are in your working branch *my_branch* and you are up-to-date. Do this
regularly so that you avoid large changes to be made to your code form outside,
which are more likely to lead to conflicts.

Contributing to the code
------------------------

Now you are free to make any changes to the code you want. To contribute to
Quantarhei for really, it is advisable to consult the :ref:`guidelines for 
contributors <how-to-contribute>`, which are part of this documentation. Here,
we concentrate only on the technical aspect of it.

At this point it is worth to mention that there is plent of excellent 
information about how to use Git on Internet and in books. If you mean your
open source involvement seriously, go on and read a lot about Git. Here we
only review basics.

When you made changes to the code, you can see which files were changed by the 
*status* subcommand of Git:

.. code:: bash

    $ git status
    
You may get the something similar to the following:

.. code:: bash

    On branch master
    Your branch is up to date with 'origin/master'.
    
    Changes not staged for commit:
      (use "git add <file>..." to update what will be committed)
      (use "git checkout -- <file>..." to discard changes in working directory)
    
    	modified:   docs/sphinx/contribute/github.rst
    
    Untracked files:
      (use "git add <file>..." to include in what will be committed)
    
    	docs/sphinx/contribute/howtocont.rst
    
    no changes added to commit (use "git add" and/or "git commit -a")

Git informs you that there is one changed file and one new untracked file. You
are advised to use `git add` command to include these files to the next commit,
or to use `git commit -a`.  

We will use the latter to record the changes you have made:

.. code:: bash

    $ git commit -a
    
Now the changes are commited, but in Git, everything is local. You changes 
still only live locally on your computer. To push them to the on-line
repository (to your private repository that you forked) you need to use the
`push` subcommand of the Git:

.. code:: bash

    $ git push

You might get an error message in some cases, which asks for `--set-upstream`
option. In this case you type:
    
.. code:: bash
 
    $ git push --set-upstream origin YOUR-BRANCH-NAME
    
and Git sets the correct target for your push. This should happen only first
time you use a given branch.
    
The push command will communicate with the remote repository,
and upload your changes.
Git will open your favourite editor - if you did not configure the editor,
be prepared to use `vi` (use ESC, : and wq! to save and leave after you
wrote a commit message describing your changes).

With your changes committed and pushed into repository you can perhaps run
tests (may you should have done this before the commit), or simply continue
working on your copy of Quantarhei. An important part of contributing is
also :ref:`writting tests <write-tests>`, which is also described in this
documentation. 


Pull Requests
-------------

The final stage of your contribution is to push the code into the main 
repository. This is something where the `Github website`_ can greatly 
help you. Consult Github's `help on pull requests`_. When you open a
pull request, you code will be reviewed and eventually (if it is good), 
it will be pulled
by the maintainer into the master branch of the Quantarhei repository.
Acceptance of a pull request may be a long process. It is greatly helped
by following advice of the :ref:`contribution guidelines <how-to-contribute>` 
and by remembering that there is no acceptance of the code with out
proper tests (which you are supposed to write with an advice of the
corresponding :ref:`guidelines <write-tests>`.)

When your pull request is accepted,
you can enjoy seeing your code appearing with the next release of Quantarhei.



.. _`github.com`: https://github.com/tmancal74/quantarhei
.. _`Git`: https://git-scm.com
.. _`Github website`: http://github.com
.. _`project Github website`: https://github.com/tmancal74/quantarhei
.. _`help on pull requests`: https://help.github.com/articles/about-pull-requests/

