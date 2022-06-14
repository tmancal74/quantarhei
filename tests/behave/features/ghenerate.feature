Feature: Ghenerate script reads a feature file and converts it into step file in default directory


    Scenario: Calling ghenerate with no arguments

        Given that Quantarhei is installed
        When I run the ghenerate script
        Then I get a simple usage message


    Scenario: Calling ghenerate with a filename

        Given that Quantarhei is installed 
        And current directory contains a feature file 
        And the default destination directory exists
        When I run the ghenerate script with the name of the feature file
        Then feature file is converted into a Python step file
        And the step file is saved into default destination directory


    Scenario Outline: Calling ghenerate with a filename and destination option

        Given that Quantarhei is installed
        And current directory contains a feature file 
        And <destination_directory> exists
        When I run <ghenerate_command> with the option specifying destination directory
        Then feature file is converted into a Python step file
        And step file is saved into the destination directory

        Examples:
        | ghenerate_command                | destination_directory |
        | ghenerate                        | ghen                  |
        | ghenerate -d ghen                | ghen                  |
        | ghenerate -d steps               | steps                 |
        | ghenerate --destination steps    | steps                 |
        | ghenerate -d some_dir            | some_dir              |
        | ghenerate --destination some_dir | some_dir              |


    Scenario Outline: if destination directory is missing, it is created before file is saved

        Given that Quantarhei is installed
        And current directory contains a feature file
        And the <destination_directory> does not exist
        When I run <ghenerate_command> with the option specifying destination directory
        Then destination directory is created
        And feature file is converted into a Python step file
        And step file is saved into the destination directory

        Examples:
        | destination_directory   |   ghenerate_command                |
        | ghen                    |   ghenerate                        |
        | ghen                    |   ghenerate -d ghen                |
        | ghen                    |   ghenerate --destination ghen     |
        | steps                   |   ghenerate -d steps               |
        | steps                   |   ghenerate --destination steps    |
        | some_dir                |   ghenerate -d some_dir            |
        | some_dir                |   ghenerate --destination some_dir |
        