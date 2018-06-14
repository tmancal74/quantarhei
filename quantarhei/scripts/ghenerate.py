"""

    Steps template generator for `behave` acceptance tests

    Author: Tomas Mancal, Charles University, Prague, Czech Republic
    email: mancal@karlov.mff.cuni.cz


"""
# standard imports
import argparse
import datetime
import time
import os
import re

# third party imports
from gherkin.token_scanner import TokenScanner
from gherkin.parser import Parser

# quantarhei
import quantarhei as qr


def parsing():
    """This function handles parsing command line arguments


    """

    descr = 'Ghenerate, the Gherkin Python Step Generator from Quantarhei'
    parser = argparse.ArgumentParser(description=descr+' ...')


    parser.add_argument("file", metavar='file', type=str,
                        help='feature file to be processed', nargs='?')

    #
    # Generator options
    #
    parser.add_argument("-v", "--version", action="store_true",
                        help="shows Quantarhei package version")
    parser.add_argument("-i", "--info", action='store_true',
                        help="shows detailed information about Quantarhei"+
                        " installation")
    parser.add_argument("-d", "--destination", type=str,
                        help="specifies destination directory for the"+
                        " generated step file")
    parser.add_argument("-n", "--no-pass", action="store_true",
                        help="empty tests should not pass (default is"
                        +" passing empty tests)")
    parser.add_argument("-f", "--start-from", type=int,
                        help="step functions will be numberred starting"
                        +" from this value")

    #
    # Parsing all arguments
    #
    args = parser.parse_args()

    #
    # show longer info
    #
    if args.info:
        qr.printlog("\n"
                    +"ghenerate: Quantarhei Gherkin Python Step Generator\n",
                    verbose=True, loglevel=0)

        if not args.version:
            qr.printlog("Package version: ", qr.Manager().version, "\n",
                        verbose=True, loglevel=0)
        return 0

    #
    # show just Quantarhei version number
    #
    if args.version:
        qr.printlog("Quantarhei package version: ", qr.Manager().version, "\n",
                    verbose=True, loglevel=0)
        return 0

    if args.destination:
        ddir = args.destination
    else:
        ddir = "ghen"

    if args.file:

        print("")
        print(descr+" ...")

        filename = args.file

    else:
        print("No file specified: quiting")
        parser.print_help()
        return 0

    steps_pass = True
    if args.no_pass:
        steps_pass = False

    k_from = 0
    if args.start_from:
        k_from = args.start_from

    try:
        with open(filename, 'r') as myfile:
            data = myfile.read()
    except:
        raise Exception("Problems reading file: "+filename)

    parser = Parser()
    try:
        feature_file = parser.parse(TokenScanner(data))
    except:
        raise Exception("Problem parsing file: "+filename+
                        " - is it a feature file?")

    try:
        children = feature_file["feature"]["children"]
    except:
        raise Exception("No scenarii or scenario outlines")

    return dict(children=children, ddir=ddir,
                steps_pass=steps_pass, filename=filename, k_from=k_from)


def analyze_children(children):
    """Analyzes children of the feature and prints info

    """

    test_strings = []
    for scenario in children:

        print(scenario["keyword"], ":", scenario["name"])

        steps = scenario["steps"]

        print("Number of steps in the scenario: ", len(steps))

        print("Following steps found:")
        for step in steps:
            print("    ", step["keyword"], "\t :\t", step["text"])
            text = step["keyword"]+": "+step["text"]
            if text not in test_strings:
                test_strings.append(text)
            else:
                print("\t*Duplicate test string:", text)
                print("\tOnly first occurence will be used to generate step")
        print("")


def check_outputfile_exists(ddir, filename):
    """Checks the existance of destination directory and the target file

    """

    (filen, ext) = os.path.splitext(os.path.basename(filename))

    if ext != ".feature":
        raise Exception("Feature file has to have the .feature extension")

    if not os.path.exists(ddir):
        os.makedirs(ddir)
    ofile = os.path.join(ddir, filen+'.py')

    print("")
    print("Generating file: ", ofile)
    print("")

    if os.path.isfile(ofile):
        answr = input("file `"+ofile+"` exists. Overwrite? (y/n) [n]: ")
        if answr.strip() != "y":
            print("Your answer is `"+answr+"`")
            print("... aborting")
            return 1

    return ofile


def write_func_def(myfile, step, textrep, args, current, k_step):
    """Write step function implementation header

    """

    if step["keyword"].strip() == "Given":
        myfile.write("\n\n#\n# Given ...\n#\n")
        myfile.write("@given('"+textrep+"')\n")
        myfile.write("def step_given_"+str(k_step)+"("+args+"):\n")
        myfile.write('    """\n\n        Given '+textrep+"\n\n")
        myfile.write('    """\n')
        current = "given"
    elif step["keyword"].strip() == "When":
        myfile.write("\n\n#\n# When ...\n#\n")
        myfile.write("@when('"+textrep+"')\n")
        myfile.write("def step_when_"+str(k_step)+"("+args+"):\n")
        myfile.write('    """\n\n        When '+textrep+"\n\n")
        myfile.write('    """\n')
        current = "when"
    elif step["keyword"].strip() == "Then":
        myfile.write("\n\n#\n# Then ...\n#\n")
        myfile.write("@then('"+textrep+"')\n")
        myfile.write("def step_then_"+str(k_step)+"("+args+"):\n")
        myfile.write('    """\n\n        Then '+textrep+"\n\n")
        myfile.write('    """\n')
        current = "then"
    elif step["keyword"].strip() == "And":
        if current == "":
            raise Exception("`And` has to be preceeded by a "+
                            "line with `Given`, `When` or"+
                            "`Then`")
        myfile.write("\n\n#\n# And ...\n#\n")
        myfile.write("@"+current+"('"+textrep+"')\n")
        myfile.write("def step_"+current+"_"+str(k_step)+"("+args+"):\n")
        myfile.write('    """\n\n        And '+textrep+"\n\n")
        myfile.write('    """\n')
    elif step["keyword"].strip() == "But":
        if current == "":
            raise Exception("`But` has to be preceeded by a "+
                            "line with `Given`, `When` or"+
                            "`Then`")
        myfile.write("\n\n#\n# But ...\n#\n")
        myfile.write("@"+current+"('"+textrep+"')\n")
        myfile.write("def step_"+current+"_"+str(k_step)+"("+args+"):\n")
        myfile.write('    """\n\n        But '+textrep+"\n\n")
        myfile.write('    """\n')
    else:
        raise Exception("unknown keyword: "+step["keyword"])

    return current


def write_header(myfile):
    """Write the output file header

    """
    tstamp = datetime.datetime.fromtimestamp(time.time()).\
    strftime('%Y-%m-%d %H:%M:%S')
    myfile.write('''"""

    Autogenerated by ghenerate script, part of Quantarhei
    http://github.com/tmancal74/quantarhei
    Tomas Mancal, tmancal74@gmai.com

    Generated on: {}

    Edit the functions below to give them desired functionality.
    In present version of `ghenerate`, no edits or replacements
    are perfomed in the feature file text.

"""

from behave import given
from behave import when
from behave import then

'''.format(tstamp))


def main():
    """Script's main function


    """

    #(children, ddir, steps_pass, filename) = parsing()
    parse_data = parsing()
    if not isinstance(parse_data, dict):
        return 0

    print("")
    print("Analyzing file: ", parse_data["filename"])
    print("Number of scenarii/scenario outlines: ",
          len(parse_data["children"]))
    print("")

    analyze_children(parse_data["children"])

    ofile = check_outputfile_exists(parse_data["ddir"], parse_data["filename"])

    with open(ofile, 'w') as myfile:

        write_header(myfile)

        test_strings = []
        k_step = parse_data["k_from"]
        for scenario in parse_data["children"]:

            steps = scenario["steps"]
            
            current = ""
            for step in steps:
                text = step["text"]
                
                if step["keyword"] != "And ":
                    prepo = step["keyword"]
                    # if keyword is not And, prepo remains
                    
                check_text = prepo+": "+step["text"]

                # the step strings should not be duplicate
                # but the same text can follow different keywords with 
                # different code
                if check_text not in test_strings:
                    test_strings.append(check_text)

                    # the 'text' is searched for variables
                    args = "context"
                    iterator = re.compile('<([^<]*)>').finditer(text)
                    for mtch in iterator:
                        # when variable found store it to args
                        # and replace the brackets
                        var = mtch.group()
                        var = var[1:len(var)-1]
                        args += ", "+var
                        text = re.sub("<"+var+">", "{"+var+"}", text)

                    #
                    # Header
                    #
                    k_step += 1
                    current = write_func_def(myfile, step, text,
                                             args, current, k_step)

                    #
                    # Body
                    #
                    if parse_data["steps_pass"]:
                        myfile.write("    pass\n")
                    else:
                        myfile.write("    assert False is True")

                # in case of duplication we just keep going
                else:
                    if step["keyword"].strip() == "Given":
                        current = "given"
                    elif step["keyword"].strip() == "When":
                        current = "when"
                    elif step["keyword"].strip() == "Then":
                        current = "then"
                    elif step["keyword"].strip() == "And":
                        if current == "":
                            raise Exception("`And` has to be preceeded by a "+
                                            "line with `Given`, `When` or"+
                                            "`Then`")
                    elif step["keyword"].strip() == "But":
                        if current == "":
                            raise Exception("`But` has to be preceeded by a "+
                                            "line with `Given`, `When` or"+
                                            "`Then`")
                    else:
                        raise Exception("unknown keyword: "+step["keyword"])


        print("... done")
        print("")

        return 0
