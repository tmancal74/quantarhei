# -*- coding: utf-8 -*-
"""
    Steps of testing of EvolutionSuperOperator class


"""
from behave import *

import quantarhei as qr
import numpy

@given('I have a general evolution superoperator U, initial density matrix R'+
       ' and a Hamiltonian H')
def step_given(context):
    
    # create test aggregate
    agg = qr.TestAggregate("dimer-2-env")
    agg.build()
    
    # get the associated time axis and the relaxation tensor and Hamiltonian
    time = agg.get_SystemBathInteraction().TimeAxis
    RR, HH = agg.get_RelaxationTensor(time, relaxation_theory="stR")
    
    context.H = HH
    
    # define and calculate evolution superoperator
    U = qr.qm.EvolutionSuperOperator(time, ham=HH, relt=RR)
    U.calculate()
    
    context.U = U
    
    # initial density matrix
    R = qr.ReducedDensityMatrix(dim=HH.dim)
    R.data[2,2] = 1.0 
    
    context.R = R

    
@when('I apply U to R in in time {time} to get density matrix Rt')
def step_when(context, time):

    U = context.U
    R = context.R
    
    t = float(time)
    
    # This happens in the site basis
    Rt = U.apply(t, R)

    context.Rt = Rt
    context.time = t


@when('I transform U and R to the eigenbasis of H to get U_ and R_,'+
      ' respectively')
def step_and(context):

    H = context.H
    U = context.U
    R = context.R

    with qr.eigenbasis_of(H):
        Ut = U.at(context.time)
        U_matrix = Ut._data
        R_matrix = R._data
        
    context.U_ = U_matrix
    context.R_ = R_matrix
    
    
@when('I apply U_ to R_ to get operator Rt_trans')
def step_and2(context):
    
    U_ = context.U_
    R_ = context.R_
    
    # This happens in the basis of B eigenstates
    Rt_trans = numpy.tensordot(U_, R_)
    
    context.Rt_trans = Rt_trans

    
@when('I transform Rt into eigenbasis of H to get Rt_')
def step_and3(context):
    
    # here we retrieve in eigenbasis of B what was calculated in site basis 
    with qr.eigenbasis_of(context.H):
        Rt_ = context.Rt._data
    
    context.Rt_ = Rt_

    
@then('Rt_ equals Rt_trans')
def step_then(context):
    
    # we compare results calculated in different bases
    numpy.testing.assert_allclose(context.Rt_trans,
                                  context.Rt_)
