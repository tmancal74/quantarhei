#
#
#  TwoDSpectrumConteiner class functionality
#
#

Feature: TwoDSpectrumContainer should store TwoDSpectra index by various methods
    TwoDSpectrum container should be able to index its content by integer values,
    by arbitrary ValueAxis (including TimeAxis) and by a list of strings 
    (like a dictionary). It can perform Fourier transform in its content if
    it is indexed by ValueAxis, TimeAxis of FrequencyAxis.
    
    
    Scenario Outline: TwoDSpectrumContainer indexes by an integer number
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        When I set the container to accept indexing by integers
        And I add the spectra to the container one by one
        Then TwoDSpectrum can be retrieved using the index <i>
        But when index is out of bounds, I get an exception
        
    Examples:
    | N  |  i  |
    | 1  |  0  |
    | 20 |  8  |
    | 12 |  12 |
    | 9  |  3  |


    Scenario Outline: TwoDSpectrumContainer indexes by ValueAxis
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        And I have a ValueAxis of lenght <N> starting from zero with certain <step>
        When I set the container to accept index by ValueAxis
        And I add the spectra to the container using values from ValueAxis
        Then TwoDSpectrum can be retrieved using values <val> from ValueAxis
        But when values are out of bounds, I get an exception
        And TwoDSpectrum can be retrieved using the index <i>
        But when index is out of bounds, I get an exception
        
    Examples:
    |  N  |  step  |  val   |  i   |
    |  5  |   0.3  |  0.9   |  4   |
    | 12  |   1.0  |  10.0  |  4   |
    | 12  |   1.0  |  15.0  |  5   |
    | 15  |   2.0  |  16.0  |  10  |
    | 15  |   2.0  |  18.0  |  20  |
    
    
    Scenario Outline: TwoDSpectrumContainer indexes by TimeAxis
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        And I have a TimeAxis of lenght <N> starting from zero with certain <step>
        When I set the container to accept index by TimeAxis
        And I add the spectra to the container using values from TimeAxis
        Then TwoDSpectrum can be retrieved using values <val> from TimeAxis
        But when values are out of bounds, I get an exception
        And TwoDSpectrum can be retrieved using the index <i>
        But when index is out of bounds, I get an exception        
        
    Examples:
    |  N  |  step  |  val   |  i   |
    |  5  |   0.3  |  0.9   |  4   |
    | 12  |   1.0  |  10.0  |  4   |
    | 12  |   1.0  |  15.0  |  5   |
    | 15  |   2.0  |  16.0  |  10  |
    | 15  |   2.0  |  18.0  |  20  |


    Scenario Outline: TwoDSpectrumContainer indexes by FrequencyAxis
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        And I have a FrequencyAxis of lenght <N> starting from zero with certain <step>
        When I set the container to accept index by FrequencyAxis
        And I add the spectra to the container using values from FrequencyAxis
        Then TwoDSpectrum can be retrieved using values <val> from FrequencyAxis
        But when values are out of bounds, I get an exception
        And TwoDSpectrum can be retrieved using the index <i>
        But when index is out of bounds, I get an exception        
        
    Examples:
    |  N  |  step  |  val   |  i   |
    |  5  |   0.3  |  0.9   |  4   |
    | 12  |   1.0  |  10.0  |  4   |
    | 12  |   1.0  |  15.0  |  5   |
    | 15  |   2.0  |  16.0  |  10  |
    | 15  |   2.0  |  18.0  |  20  |
    

    Scenario Outline: TwoDSpectrumContainer indexes by strings
        Given that I have <N> TwoDSpectrum objects
        And I have an empty TwoDSpectrum container
        And I have a list of strings of lenght <N>
        When I set the container to accept index by strings
        And I add the spectra to the container using values from the list of strings
        Then TwoDSpectrum can be retrieved using values from the list of strings
        But when values are not in the list of strings, I get an exception
        And TwoDSpectrum can be retrieved using the index <i>
        But when index is out of bounds, I get an exception        
        
    Examples:
    |  N  |  i   |
    |  5  |  4   |
    | 12  |  4   |
    | 12  |  5   |
    | 15  |  10  |
    | 15  |  20  |


    Scenario Outline: TwoDSpectrumContainer can do Fourier transform when it is indexed by ValuesAxis
        Given that I have a TwoDSpectrumContainer containing <N> spectra indexed by ValueAxis
        When I calculate Fourier transform on the container
        Then I get correct pointwise Fourier transform of the spectra
        And the TwoDSpectrum container will be indexed by ValueAxis with frequencies corresponding to the original ValueAxis
        
    Examples:
    | N    |
    | 10   |
    | 20   |
    | 100  |


    Scenario Outline: TwoDSpectrumContainer can do Fourier transform when it is indexed by TimeAxis
        Given that I have a TwoDSpectrumContainer containing <N> spectra indexed by TimeAxis
        When I calculate Fourier transform on the container
        Then I get correct pointwise Fourier transform of the spectra
        And the TwoDSpectrum container will be indexed by FrequencyAxis which corresponds to the original TimeAxis
        
    Examples:
    | N    |
    | 10   |
    | 20   |
    | 100  |
    

    Scenario Outline: TwoDSpectrumContainer can do Fourier transform when it is indexed by FrequencyAxis
        Given that I have a TwoDSpectrumContainer containing <N> spectra indexed by FrequencyAxis
        When I calculate Fourier transform on the container
        Then I get correct pointwise Fourier transform of the spectra
        And the TwoDSpectrum container will be indexed by TimeAxis which corresponds to the original FrequencyAxis
        
    Examples:
    | N    |
    | 10   |
    | 20   |
    | 100  |      
        