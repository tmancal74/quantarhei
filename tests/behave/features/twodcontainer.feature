#
#
#  TwoDSpectrumConteiner class functionality
#
#

Feature: TwoDSpectrumContainer should store TwoDSpectra index by various methods
    TwoDSpectrum container should be able to index its content by integer values,
    by arbitrary ValueAxis (including TimeAxis) and by a list of strings 
    (like a dictionary).
    
    
    Scenario Outline: TwoDSpectrumContainer indexes by an integer number
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        When I set the container to accept indexing by integers
        And I add the spectra to the container one by one
        Then TwoDSpectrum can be retrieved using the index <i>
        But the when index is out of bounds, I get an exception
        
    Examples:
    | N  |  i  |
    | 1  |  0  |
    | 20 |  8  |
    | 12 |  12 |
    | 9  |  3  |
