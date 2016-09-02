# -*- coding: utf-8 -*-


import aloe

@aloe.step(r"Given reorganization energy (\d+(?:\.\d+)?) and correlation time (\d+(?:\.\d+)?)")
def reorganization_energy_correlation_time(self, reorg, ctime):
    print("\n", reorg, ctime)

@aloe.step(r"""When I calculate the ([^"]*) correlation function""")
def correlation_function_of_type(self, ctype):
    print(ctype)

@aloe.step(r"""Then I get data from the file ([^"]*)""")
def compare_data_with_file(self, file):
    print(file)