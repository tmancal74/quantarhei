Feature: qrhei fetch subcommand fetches examples to the working directory

Scenario: All examples are fetchable
Given that I have a list of examples from qrhei list
When I fetch all examples one by one
Then examples are all fetchable