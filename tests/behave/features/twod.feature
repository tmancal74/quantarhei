Feature: TwoDSpectrum can store data in different level of details.
    There is a hierarchy of storage. The lowest level detail (highest resolution) is represented by
    individual Liouville pathways. When in the mode of storing pathways,
    each spectrum is store with pathways character "R1g", "R2f", etc. and
    with a unique tag. Individual pathways can be retrieved using tags, and
    spectra corresponding to group of pathways, such as "ESA", "GSB", etc. or
    "total", "REPH", "NONR" can be calculated. 
    
    
    
    Scenario Outline: TwoDSpectrum can store data by liouville pathways 
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        When I create a new TwoDSpectrum object
        Then I can save 2D data using type and tag
        And I can retrieve spectra by type and tag
        And I can retrieve sum of spectra of a given type <type>
        And I can retrieve sum of spectra of a given process <process>
        And I can retrieve sum of spectra of a given signal <signal>
        And I can retrieve total spectrum
        
        
    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   REPH    |
        | R1g  |  GSB      |   NONR    |
        | R3g  |  ESA      |   REPH    |
        | R4g  |  SE       |   NONR    |
        | R1fs |  GSB      |   REPH    |
        | R2fs |  ESA      |   NONR    |
        | R3fs |  SE       |   REPH    |
         


    Scenario Outline: TwoDSpectrum storing data in form of pathways can convert into storage of types of pathways
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        When I convert the storage mode into the one storing only types of pathways
        Then I can retrieve sum of spectra of a given type <type>
        And I can retrieve sum of spectra of a given process <process>
        And I can retrieve sum of spectra of a given signal <signal>
        And I can retrieve total spectrum

    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   REPH    |
        | R1g  |  GSB      |   NONR    |
        | R3g  |  ESA      |   REPH    |
        | R4g  |  SE       |   NONR    |
        | R1fs |  GSB      |   REPH    |
        | R2fs |  ESA      |   NONR    |
        | R3fs |  SE       |   REPH    |

  
    
    Scenario Outline: TwoDSpectrum storing data in form of pathways can convert into storage of spectra related to processes
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        When I convert the storage mode into the one storing spectra of processes
        Then I can retrieve sum of spectra of a given process <process>
        And I can retrieve total spectrum
        But when I try to retrieve signal <signal> I get an exception

    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   REPH    |
        | R1g  |  GSB      |   NONR    |
        | R3g  |  ESA      |   REPH    |
        | R4g  |  SE       |   NONR    |
        | R1fs |  GSB      |   REPH    |
        | R2fs |  ESA      |   NONR    |
        | R3fs |  SE       |   REPH    |
        


    Scenario Outline: TwoDSpectrum storing data in form of pathways can convert into storage of rephasing and non-rephasing signals
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        When I convert the storage mode into the one storing rephasing and non-rephasing signals
        Then I can retrieve sum of spectra of a given signal <signal>
        And I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception

    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   REPH    |
        | R1g  |  GSB      |   NONR    |
        | R3g  |  ESA      |   REPH    |
        | R4g  |  SE       |   NONR    |
        | R1fs |  GSB      |   REPH    |
        | R2fs |  ESA      |   NONR    |
        | R3fs |  SE       |   REPH    |
        
        
        
    Scenario Outline: TwoDSpectrum storing data in form of pathways can convert into storage of total spectrum
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception


    Examples:
        | type |  process  |   signal  |
        | R2g  |  SE       |   REPH    |
        | R1g  |  GSB      |   NONR    |
        | R3g  |  ESA      |   REPH    |
        | R4g  |  SE       |   NONR    |
        | R1fs |  GSB      |   REPH    |
        | R2fs |  ESA      |   NONR    |
        | R3fs |  SE       |   REPH    |
        
        
    Scenario Outline: TwoDSpectra storing data accoring to pathway types can convert into storage of spectra related to processes
        Given that I have data corresponding to Liouville pathway types
        And I create a new TwoDSpectrum object
        And I save 2D data using type information
        When I convert the storage mode into the one storing spectra of processes
        Then I can retrieve sum of spectra of a given process <process>
        And I can retrieve total spectrum
        But when I try to retrieve signal <signal> I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |


    Scenario Outline: TwoDSpectra storing data accoring to pathway types can convert into storage of spectra related to signals
        Given that I have data corresponding to Liouville pathway types
        And I create a new TwoDSpectrum object
        And I save 2D data using type information
        When I convert the storage mode into the one storing spectra of rephasing and non-rephasing signals
        Then I can retrieve sum of spectra of a given signal <signal>
        And I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |


    Scenario Outline: TwoDSpectra storing data accoring to pathway types can convert into storage of total spectrum only
        Given that I have data corresponding to Liouville pathway types
        And I create a new TwoDSpectrum object
        And I save 2D data using type information
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception
        And when I try to retrieve pathways I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |


    Scenario Outline: TwoDSpectra storing data accoring to signal type can convert into storage of total spectrum only
        Given that I have data corresponding to signal types
        And I create a new TwoDSpectrum object
        And I save 2D data using signal types
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception
        And when I try to retrieve pathways I get an exception
        And when I try to retrieve types of spectra I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |


    Scenario Outline: TwoDSpectra storing data accoring to processes can convert into storage of total spectrum only
        Given that I have data corresponding to processes
        And I create a new TwoDSpectrum object
        And I save 2D data using processes
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception
        And when I try to retrieve pathways I get an exception
        And when I try to retrieve types of spectra I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |

        
    Scenario Outline: TwoDSpectra storing only total data 
        Given that I have data corresponding to total spectrum
        And I create a new TwoDSpectrum object
        When I save 2D data of the total spectrum
        Then I can retrieve total spectrum
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception
        And when I try to retrieve pathways I get an exception
        And when I try to retrieve types of spectra I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |


#
# Trimming
#

    Scenario Outline: TwoDSpectra storing data accoring to signal type can be trimmed
        Given that I have data corresponding to signal types
        And I create a new TwoDSpectrum object
        And I save 2D data using signal types
        And I trim the data to half the length of the axes
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum with half the length of the axes
        But when I try to retrieve process <process> I get an exception
        And when I try to retrieve signal <signal> I get an exception
        And when I try to retrieve pathways I get an exception
        And when I try to retrieve types of spectra I get an exception

    Examples:
        |  process  |   signal  |
        |  SE       |   REPH    |
        |  GSB      |   NONR    |
        |  ESA      |   REPH    |
        |  SE       |   NONR    |
        |  GSB      |   REPH    |
        |  ESA      |   NONR    |
        |  SE       |   REPH    |
        
        
    Scenario: TwoDSpectrum storing data by liouville pathways can be trimmed
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        And I trim the data to half the length of the axes
        When I convert the storage mode into the storing total spectrum
        Then I can retrieve total spectrum with half the length of the axes
        

#
# Operations on the spectrum
#

    Scenario Outline: TwoDSpectra in any storage mode can be devided by a number
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        And I convert the storage mode into <storage_mode>
        When I devide the spectrum by number <number>    
        Then I can retrieve total spectrum
        
    Examples:
        | storage_mode | number |
        |  pathways    |  2.0   |
        |  pathways    |  3     |
        |  pathways    |  10.1  |
        |  types       |  2.0   |
        |  types       |  3     |
        |  types       |  10.1  |
        |  signals     |  3     |
        |  signals     |  10.1  |
        |  processes   |  3     |
        |  processes   |  10.1  |
        |  off         |  3     |
        |  off         |  10.1  |


    Scenario Outline: From TwoDSpectra in any storage mode one can get value at a point
        Given that I have data corresponding to individual Liouville pathways in 2D spectrum
        And I create a new TwoDSpectrum object
        And I save 2D data using type and tag
        And I convert the storage mode into <storage_mode>
        When I devide the spectrum by number <number>    
        Then I can retrieve total spectrum at a point [<px>, <py>]
        
    Examples:
        | storage_mode | number | px    | py    |
        |  pathways    |  2.0   | 0.0   | 0.0   |
        |  pathways    |  3     | 0.0   | 0.0   |
        |  pathways    |  10.1  | 0.0   | 0.0   |
        |  types       |  2.0   | 0.0   | 0.0   |
        |  types       |  3     | 0.0   | 0.0   |
        |  types       |  10.1  | 0.0   | 0.0   |
        |  signals     |  3     | 0.0   | 0.0   |
        |  signals     |  10.1  | 0.0   | 0.0   |
        |  processes   |  3     | 0.0   | 0.0   |
        |  processes   |  10.1  | 0.0   | 0.0   |
        |  off         |  3     | 0.0   | 0.0   |
        |  off         |  10.1  | 0.0   | 0.0   |
 
                              