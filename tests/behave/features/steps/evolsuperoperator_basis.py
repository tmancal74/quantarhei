# -*- coding: utf-8 -*-
"""
    Steps of testing of EvolutionSuperOperator class


"""
from behave import *

import quantarhei as qr
import numpy

@given('I have a general evolution superoperator U, initial density matrix R'+
       ' and a Hamiltonian H')
def step_given_1(context):
    
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
def step_when_2(context, time):

    U = context.U
    R = context.R
    
    t = float(time)
    
    # This happens in the site basis
    Rt = U.apply(t, R)

    context.Rt = Rt
    context.time = t


@when('I transform U and R to the eigenbasis of H to get U_ and R_,'+
      ' respectively')
def step_when_3(context):

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
def step_when_4(context):
    
    U_ = context.U_
    R_ = context.R_
    
    # This happens in the basis of B eigenstates
    Rt_trans = numpy.tensordot(U_, R_)
    
    context.Rt_trans = Rt_trans

    
@when('I transform Rt into eigenbasis of H to get Rt_')
def step_when_5(context):
    
    # here we retrieve in eigenbasis of B what was calculated in site basis 
    with qr.eigenbasis_of(context.H):
        Rt_ = context.Rt._data
    
    context.Rt_ = Rt_

    
@then('Rt_ equals Rt_trans')
def step_then_6(context):
    
    # we compare results calculated in different bases
    numpy.testing.assert_allclose(context.Rt_trans,
                                  context.Rt_)


#
# Given ...
#
@given('I have a Hamiltonian H, Lidblad form L and initial density matrix R')
def step_given_7(context):
    """

        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R

    """
    # create test aggregatedimer
    agg = qr.TestAggregate("trimer-2")
    agg.build()
    
    # get the associated time axis and the relaxation tensor and Hamiltonian
    time = qr.TimeAxis(0, 320, 1.0)
    time2 = qr.TimeAxis(0, 32, 10.0)
    context.time = time
    context.time2 = time2
    
    HH = agg.get_Hamiltonian()
    context.H = HH
    
    SBI = qr.qm.TestSystemBathInteraction(name="trimer-2-lind")
    LL = qr.qm.LindbladForm(HH, SBI)
    
    context.L= LL
    
    # initial density matrix
    R = qr.ReducedDensityMatrix(dim=HH.dim)
    R.data[2,2] = 1.0 
    
    context.R = R
#
# When ...
#
@when('I calculate evolution superoperator using H and L in sitebasis')
def step_when_8(context):
    """

        When I calculate evolution superoperator using H and L in sitebasis

    """
    # define and calculate evolution superoperator
    time2 = context.time2
    LL = context.L
    HH = context.H
   
    U = qr.qm.EvolutionSuperOperator(time2, ham=HH, relt=LL)
    U.set_dense_dt(10)
    U.calculate()
    
    context.U = U


#
# And ...
#
@when('I calculate dynamics of R using H and L to get R1')
def step_when_9(context):
    """

        And I calculate dynamics of R using H and L to get R1

    """
    time = context.time
    HH = context.H
    LL = context.L
    
    prop = qr.ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                             RTensor=LL)

    R = context.R
    
    R1 = prop.propagate(R)

    context.R1 = R1
    
#
# And ...
#
@when('I apply the evolution superoperator to R to get R2 at times {t_prop}')
def step_when_10(context, t_prop):
    """

        And I apply the evolution superoperator to R to get R2

    """
    t = float(t_prop)

    R = context.R
    U = context.U
    
    R2 = U.apply(t, R)

    context.R2 = R2
    
#
# Then ...
#
@then('R1 equals R2 at times {t_prop}')
def step_then_11(context, t_prop):
    """

        Then R1 equals R2 at times {t_prop}

    """
    t = float(t_prop)
    R1_t = context.R1.at(t)
    R2 = context.R2
    
    numpy.testing.assert_allclose(R1_t.data, R2.data, rtol=1.0e-12, atol=1.0e-12)


#
# When ...
#
@when('I calculate evolution superoperator using H and L in exciton basis')
def step_when_12(context):
    """

        When I calculate evolution superoperator using H and L in exciton basis
    """
    # define and calculate evolution superoperator
    time2 = context.time2
    LL = context.L
    HH = context.H
   
    with qr.eigenbasis_of(HH):
        U = qr.qm.EvolutionSuperOperator(time2, ham=HH, relt=LL)
        U.set_dense_dt(10)
        U.calculate()
    
    context.U = U

#
# When ...
#
@when('I calculate evolution superoperator in site basis using H and L in exciton basis')
def step_when_12(context):
    """

        When I calculate evolution superoperator in site basis using H and L in exciton basis

    """
    # define and calculate evolution superoperator
    time2 = context.time2
    LL = context.L
    HH = context.H
   
    with qr.eigenbasis_of(HH):
        U = qr.qm.EvolutionSuperOperator(time2, ham=HH, relt=LL)
        U.set_dense_dt(10)
    
    U.calculate()
    
    context.U = U