# -*- coding: utf-8 -*-
"""
    Steps of testing of SuperOperator class


"""
from behave import *

import quantarhei as qr
import numpy

@given('I have a general superoperator S and two operators A and B')
def step_given(context):
    
    S = qr.qm.TestSuperOperator("dim-3-AOA")
    A = qr.Hamiltonian(data=[[0.0, 0.1, 0.0], 
                             [0.1, 1.0, 0.2], 
                             [0.0, 0.2, 1.2]])
    B = qr.Hamiltonian(data=[[0.0, 0.3, 0.1],
                             [0.3, 1.0, 0.4],
                             [0.1, 0.4, 2.0]])
    
    context.S = S
    context.A = A
    context.B = B 

    
@when('I apply S to A to get operator C')
def step_when(context):

    S = context.S
    A = context.A
    
    # This happens in the site basis
    C = S.apply(A)

    context.C = C


@when('I transform S and A to the eigenbasis of B to get S_ and A_, respectively')
def step_and(context):

    S = context.S
    A = context.A
    B = context.B

    with qr.eigenbasis_of(B):
        S_matrix = S._data
        A_matrix = A._data
        
    context.S_ = S_matrix
    context.A_ = A_matrix
    
    
@when('I apply S_ to A_ to get operator D')
def step_and2(context):
    
    S_ = context.S_
    A_ = context.A_
    
    # This happens in the basis of B eigenstates
    D = numpy.tensordot(S_, A_)
    
    context.D = D

    
@when('I transform C into eigenbasis of B to get C_')
def step_and3(context):
    
    # here we retrieve in eigenbasis of B what was calculated in site basis 
    with qr.eigenbasis_of(context.B):
        C_ = context.C._data
    
    context.C_ = C_

    
@then('C_ equals D')
def step_then(context):
    
    # we compare results calculated in different bases
    numpy.testing.assert_allclose(context.C_, context.D)