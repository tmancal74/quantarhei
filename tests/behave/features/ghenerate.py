#! /Users/tomas/anaconda3/bin/python


from gherkin.token_scanner import TokenScanner
from gherkin.parser import Parser

print("")
print("Ghenerate, the Gherkin python step generator ...")


print("")
print("Found feature files:")

filen = "superoperator"
filename = filen+".feature"

print("   ", filename)

with open(filename, 'r') as myfile:
  data = myfile.read()

parser = Parser()
feature_file = parser.parse(TokenScanner(data))

children = feature_file["feature"]["children"]

print("")
print("Analyzing file: ", filename)
print("")
print("Number of features: ", len(children))

steps = children[0]["steps"]

print("Number of steps in the feature: ", len(steps))

print("Following steps found:")
for step in steps:
   print("    ", step["keyword"], "\t :\t", step["text"])

ofile = filen+".py"

print("")
print("Generating files: ")
print("  ghen/"+ofile)
print("")
print("... done")
print("")
 








