Feature: Aggregate object can return PureDephasing object


    Scenario: Aggregate can return empty PureDephasing object
        Given that I have an aggregate of three molecules with no dephasing set
        When I ask the aggregate to return PureDephasing object
        Then I get an empty PureDephasing object
        

    Scenario: Aggregate can return PureDephasing object with pure dephasing rates
        Given that I have an aggregate of three molecules with dephasing rates set
        When I ask the aggregate to return PureDephasing object
        Then I get a PureDephasing object with some rates

    Scenario: Homo-dimer has electronic PureDephasing equal to zero
        Given that I have an aggregate of two identical molecules with dephasing rates set
        When I ask the aggregate to return PureDephasing object
        Then I get a PureDephasing object with electronic dephasing rates equal to zero

    Scenario: Dimer aggregate with nuclear mode can return PureDephasing object with pure dephasing rates
        Given that I have an aggregate of two molecules with nuclear mode each and with electronic dephasing rates set
        When I ask the aggregate to return PureDephasing object
        Then I get a PureDephasing object with some rates
