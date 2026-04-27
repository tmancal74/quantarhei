from __future__ import annotations

import argparse
import traceback
from typing import Any

import quantarhei as qr


def code_hamiltonian(INP: Any, status: Any) -> str:

    code = "# Ahoj"

    return code



def tasks(tlist: Any, INP: Any) -> str:
    """Go through the list of tasks

    """
    status: list[Any] = []
    code = ""

    if isinstance(tlist, str):
        tlist = tlist.split()
    print("Identified ", len(tlist), "tasks")
    for task in tlist:
        print("Processing: ", task)

        if task == "hamiltonian":
            code += code_hamiltonian(INP, status)

    return code


def do_process(args: Any) -> None:
    """Process arguments of the script

    """
    if args.filename:
        fname = args.filename[0]
        print("Processing file: ", fname)

        try:
            INP = qr.Input(fname)
            code = tasks(INP.tasks, INP)
        except Exception:
            print(traceback.format_exc())
            return

    if args.output:

        print("Script saved as: ", args.output)

    else:

        print("\nResulting code:\n")
        print(code)



def main() -> None:

    global parser_list
    global parser_fetch

    parser = argparse.ArgumentParser(
            description='Quantarhei Task Compiler')


    #
    # Driver options
    #
    parser.add_argument("-o", "--output", metavar="output", type=str,
                        help="output file")
    parser.add_argument("filename", metavar='filename', type=str,
                        help='task file to be processed', nargs=1)
    parser.set_defaults(func=do_process)

    #
    # Parsing all arguments
    #
    args = parser.parse_args()

    print("\nQTask: Quantarhei Task Compilator\n")

    try:
        if args.func:
            args.func(args)
    except AttributeError:
        parser.error("No arguments provided")
