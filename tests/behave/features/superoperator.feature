#
#  SuperOperator testing 
#

Feature: Super operator basis transformations
    Super operators need to transform correctly under basis transformation

    Scenario: Application of superoperator in different basis leads to the same result
        Given I have a general superoperator S and two operators A and B
         When I apply S to A to get operator C
          And I transform S and A to the eigenbasis of B to get S_ and A_, respectively
          And I apply S_ to A_ to get operator D
          And I transform C into eigenbasis of B to get C_
         Then C_ equals D
