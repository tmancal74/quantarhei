Feature: TwoDSpectrum can store data in different level of details.
    There is a hierarchy of storage. The lowest detail is represented by
    individual Liouville pathways. When in the mode of storing pathways,
    each spectrum is store with pathways character "R1g", "R2f", etc. and
    with a unique tag. Individual pathways can be retrieved using tags, and
    spectra corresponding to group of pathways, such as "ESA", "GSB", etc. or
    "Tot", "Reph", "Nonr" can be calculated. 
    
    
    
    Scenario Outline: TwoDSpectrum can store data by liouville pathways 
        Given that I have data corresponding to individual liouville pathways in 2D spectrum
        When I create a new TwoDSpectrum object
        Then I can save 2D data using type and tag
        And I can retrieve spectra by tag
        And I can retrieve sum of spectra of a given type <type>
        And I can retrieve sum of spectra of a given process <process>
        And I can retrieve sum of spectra of a given signal <signal>
        
        
    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   Reph    |
        | R1g  |  GSB      |   Nonr    |
        | R3g  |  ESA      |   Reph    |
        | R4g  |  SE       |   Nonr    |
        | R1fs |  GSB      |   Reph    |
        | R2fs |  ESA      |   Nonr    |
        | R3f  |  SE       |   Reph    |
         