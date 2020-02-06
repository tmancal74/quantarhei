# -*- coding: utf-8 -*-
#
#  Example for the parallel iterator and collector:
#      block_distributed_list
#      collect_distributed_list_data
#
#  The iterator `block_distributed_list` distributes a list among different
#  processes. The list is defined everywhere (in all processes), but only its
#  section is given to each process by the iterator. If one calculates some data
#  based on the section of the list (it is assumed that it is one set of data
#  per the list member), and stores them on a container defined in the process,
#  one can collect the data from the local containers into one container on
#  the master nod by the function `collect_distributed_list_data`. The user
#  supplies an empty container and the local container in a two membered list,
#  and provides functions to set and retrieve the data to and from
#  the container, based on the data and a tag. Also the mapping between
#  tags and the integer indices starting from zero can be submitted.
#
#
#
#  Usage:
#
#  > qrhei run -p -n NUMBER_OF_PROVESSES ex_300_ParallelIterator.py
#
#
#
import numpy
import quantarhei as qr


#
# The container is represented by a dictionary in this example. We define
# a setter and retriever function for a dictionary
#
def setter(cont, tag, data):
    """Setter function setting data to a primitive container (here a dictionary)
    """
    cont[tag] = data

def retriever(cont, tag):
    """Function retrieving the data from the container (here a dictionary)
    """
    return cont[tag]

#
# This list will be distributed
#
lst = [1,2,3,4,5,6,7]

#
# tags under which the data is saved to the container
#
tags = ["a","b","c","d","e","f","g"]

# container for local calculation
cont = dict()
# global container to hold the collected results
collected_cont = dict()

# two membered list of containers
containers = [collected_cont, cont]

# config holds information about parallel environment
config = qr.Manager().get_DistributedConfiguration()

# iteration over the list
for k, a in qr.block_distributed_list(lst, return_index=True):

    #print(config.rank, a)
    # calculation of the data
    b = numpy.zeros(2, dtype=qr.COMPLEX)
    b[0] = 2.0*a
    b[1] = 3.2*a + 1.3
    # storing data under an appropriate tag
    cont[tags[k]] = b

# here we collect the data to a single container available on the master nod
qr.collect_block_distributed_data(containers, setter, retriever, tags=tags)

# printing the result
#if config.rank == 0:
if True:
    print(config.rank, collected_cont)

