from __future__ import annotations

import ast
import json
import os
import re
from typing import Any

import yaml  # type: ignore[import-untyped]


#
# The code below is here to test against arbitrary code execution during
# evaluation of math expressions in yaml configuration files
#
def ahoj(a: float) -> None:
    import numpy

    print("A teď dělám, co mě napadne")
    print(numpy.cos(a))


a = 5.0


def expr(code: str, context: dict[str, Any] | None = None) -> Any:
    """Eval a math expression and return the result

    To be honest, I have no idea how this works
    I got this piece of code from:
        https://gist.github.com/mariocesar/4af926b01c5496563787


    Examples
    --------
    This example should ultimately fail

    >>> print(expr("ahoj(a)"))
    A teď dělám, co mě napadne
    0.283662185463
    None

    """
    if not context:
        context = {}
    code = code.format(**context)

    # An especial case for my own when using percents values.
    # TODO: Fail if you are not comparing same type value like "50 > 20%" haves to fail
    code = re.sub("%", "", code)

    # FIXME: is there any way to protect against arbitrary code execution?
    expr = ast.parse(code, mode="eval")
    code_object = compile(expr, "<string>", "eval")

    return eval(code_object)


class Input:
    """A class representing a json or yaml input file"""

    # YAML-injected attributes — set dynamically via setattr in __init__
    define_usecases: Any
    tasks: Any
    script: Any
    path: Any

    def __init__(
        self,
        file_or_dict: str | dict[str, Any],
        math_allowed_in: list[Any] | None = None,
        show_input: bool = False,
    ) -> None:

        if math_allowed_in is None:
            math_allowed_in = []
        conversions = math_allowed_in
        self._math_allowed_in = []

        if isinstance(file_or_dict, str):
            filename = file_or_dict

            # get the file extension
            extsplt = os.path.splitext(filename)
            ext = extsplt[1]

            #
            # Reading configuration file
            #

            # json
            if ext in [".json"]:
                with open(filename) as f:
                    self.data = json.load(f)

            # yaml
            elif ext in [".yaml", ".yml"]:
                with open(filename) as f:
                    self.data = yaml.safe_load(f)

            else:
                Exception("Unknown file type")

            self._from_file = True

        elif isinstance(file_or_dict, dict):
            self.data = file_or_dict
            self._from_file = False

        else:
            raise Exception(
                "Input file name (string) or a dictionary hs to be specified"
            )

        data = self.data
        for piece in data:
            if piece == "_math_allowed_in":
                self._math_allowed_in += data[piece]
            else:
                setattr(self, piece, data[piece])

        conversions += self._math_allowed_in

        #
        # Converts selected strings into floats by evaluating them
        # as mathematical expressions
        #
        for cnv in conversions:
            if isinstance(cnv, list) or isinstance(cnv, tuple):
                prop = cnv[0]
                keys = cnv[1]

                val = getattr(self, prop)

                if isinstance(val, dict):
                    for key in keys:
                        # if a variable is listed but not found, we ignor it
                        try:
                            valstr = val[key]
                            valstr = self.string_2_float_prop(valstr)
                            val[key] = valstr
                        except KeyError:
                            pass

                else:
                    raise Exception("... must be a dictionary")

            else:
                val = getattr(self, cnv)
                if isinstance(val, list):
                    self.strings_2_floats_lists(val)
                else:
                    val = self.string_2_float_prop(val)
                    setattr(self, cnv, val)

        self.find_usecases()

        if show_input:
            print("\nInput file summary:\n")
            print(self.data)

    def strings_2_floats_dictionary(
        self, dictionary: dict[str, Any], keys: list[str]
    ) -> dict[str, Any]:
        """Converts selected keys of a dictionary from string expression to float"""
        ndict = {}
        for key in dictionary:
            val = dictionary[key]
            if key in keys:
                if isinstance(val, str):
                    val = expr(val)

            ndict[key] = val
        return ndict

    def strings_2_floats_lists(self, inlist: list[Any]) -> None:
        """Converts every string in a list into float"""
        N = len(inlist)

        for k in range(N):
            val = inlist[k]
            val = self.string_2_float_prop(val)
            inlist[k] = val

    def string_2_float_prop(self, prop: Any) -> Any:
        """If the submitted object is a string it is converted to float"""
        if isinstance(prop, str):
            val = expr(prop)
        else:
            val = prop
        return val

    def dump_yaml(self, filename: str) -> None:

        yaml.dump(self.data, open(filename, "w"), default_flow_style=False)

    def find_usecases(self) -> None:
        """Finds replacements for values defined in the input file

        This function replaces values of some input parameters by the values
        predefined in the "define_usecases" construct.

        Parameters
        ----------
        inpt: input
            The input object constructed from the configuration file.

        """
        try:
            defined_usecases = self.define_usecases["usecases"]
            definitions = self.define_usecases["definitions"]
            du_Number = len(defined_usecases)
            print()
            print("Found", du_Number, "defined usecase(s)")
        except (AttributeError, KeyError):
            defined_usecases = []

        for duname in defined_usecases:
            # print("Usecase", "'"+duname+"'","has possible values:")
            defs = definitions[duname]["values"]
            vars = definitions[duname]["variables"]
            cass = definitions[duname]["cases"]
            for df in defs:
                # print(df,":")
                kk = 0
                for var in vars:
                    # print(var, "=", cass[df][kk])
                    kk += 1

            try:
                usecase_value = getattr(self, duname)
            except AttributeError:
                usecase_value = None

            if usecase_value is not None:
                print()
                print("Replacing for usecase:", duname)
                print("with the value of:", usecase_value)
                kk = 0
                for var in vars:
                    print(var, "=", cass[usecase_value][kk])
                    setattr(self, var, cass[usecase_value][kk])
                    kk += 1
