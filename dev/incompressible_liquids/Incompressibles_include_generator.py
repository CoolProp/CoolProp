'''
Created on 8 Sep 2014

@author: jowr
'''

from __future__ import division, absolute_import, print_function
import os
from CPIncomp import getCoefficientFluids, getDigitalFluids, getMelinderFluids,\
    getPureFluids, getSecCoolFluids

# See http://stackoverflow.com/questions/11347505/what-are-some-approaches-to-outputting-a-python-data-structure-to-restructuredte
def make_table(grid):
    max_cols = [max(out) for out in map(list, zip(*[[len(item) for item in row] for row in grid]))]
    rst = table_div(max_cols, 1)

    for i, row in enumerate(grid):
        header_flag = False
        if i == 0 or i == len(grid)-1: header_flag = True
        rst += normalize_row(row,max_cols)
        rst += table_div(max_cols, header_flag )
    return rst

def table_div(max_cols, header_flag=1):
    out = ""
    if header_flag == 1:
        style = "="
    else:
        style = "-"

    for max_col in max_cols:
        out += max_col * style + " "

    out += "\n"
    return out


def writeTextToFile(path,text):
    f = open(path, "w")
    f.write(text)
    f.close()
    return True

def writeTableToFile(path,table):
    return writeTextToFile(path, make_table(table))


def normalize_row(row, max_cols):
    r = ""
    for i, max_col in enumerate(max_cols):
        r += row[i] + (max_col  - len(row[i]) + 1) * " "

    return r + "\n"



if __name__ == '__main__':

    FLUID_INFO_FOLDER=os.path.join("..","..","Web","fluid_properties")

    FLUID_INFO_INCOMP_LIST=os.path.join(FLUID_INFO_FOLDER,"Incompressibles_Fluids.txt")



    coeffObjs = getCoefficientFluids()
    digitObjs = getDigitalFluids()
    melinObjs = getMelinderFluids()
    pureFObjs = getPureFluids()
    secCoObjs = getSecCoolFluids()

    allObjs   = []
    allObjs += coeffObjs[:]
    allObjs += digitObjs[:]
    allObjs += melinObjs[:]
    allObjs += pureFObjs[:]
    allObjs += secCoObjs[:]

    testTable = []
    testTable.append(['Name', 'Description', 'Reference']) # Headline

    for fluid in allObjs:
        testTable.append([fluid.name, fluid.description, fluid.reference])


    writeTableToFile(FLUID_INFO_INCOMP_LIST, testTable)
