'''
Created on 8 Sep 2014

@author: jowr
'''

from __future__ import division, absolute_import, print_function
import os
from CPIncomp import getCoefficientFluids, getDigitalFluids, getMelinderFluids,\
    getPureFluids, getSecCoolFluids
from CPIncomp.DataObjects import SolutionData
from CPIncomp.BaseObjects import IncompressibleData

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
    print("Writing to file: {0}".format(path))
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

    FLUID_INFO_FOLDER=os.path.abspath(os.path.join("..","..","Web","fluid_properties"))

    FLUID_INFO_MASS_LIST=os.path.join(FLUID_INFO_FOLDER,"mass-based-fluids.txt")
    FLUID_INFO_MOLE_LIST=os.path.join(FLUID_INFO_FOLDER,"mole-based-fluids.txt")
    FLUID_INFO_VOLU_LIST=os.path.join(FLUID_INFO_FOLDER,"volume-based-fluids.txt")
    FLUID_INFO_PURE_LIST=os.path.join(FLUID_INFO_FOLDER,"pure-fluids.txt")



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

    # Parse the objects and make lists according to the data sources
    dataSourceObjs = []
    equaSourceObjs = []
    coefSourceObjs = []
    unknSourceObjs = []

    # Make other lists for pure solution fluids
    massBasedObjs = []
    moleBasedObjs = []
    voluBasedObjs = []
    pureBasedObjs = []
    unknBasedObjs = []

    for fluid in allObjs:

        if fluid.density.source==IncompressibleData.SOURCE_DATA:
            dataSourceObjs += [fluid]
        elif fluid.density.source==IncompressibleData.SOURCE_EQUATION:
            equaSourceObjs += [fluid]
        elif fluid.density.source==IncompressibleData.SOURCE_COEFFS:
            coefSourceObjs += [fluid]
        else:
            unknSourceObjs += [fluid]

        if fluid.xid==SolutionData.ifrac_mass:
            massBasedObjs += [fluid]
        elif fluid.xid==SolutionData.ifrac_mole:
            moleBasedObjs += [fluid]
        elif fluid.xid==SolutionData.ifrac_volume:
            voluBasedObjs += [fluid]
        elif fluid.xid==SolutionData.ifrac_pure:
            pureBasedObjs += [fluid]
        else:
            unknBasedObjs += [fluid]

    # After all the list got populated, we can process the entries
    # and generate some tables

    header = ['Name', 'Description', 'Reference']

    #
    objLists = [pureBasedObjs,massBasedObjs,moleBasedObjs,voluBasedObjs]
    filLists = [FLUID_INFO_PURE_LIST,FLUID_INFO_MASS_LIST]
    filLists +=[FLUID_INFO_MOLE_LIST,FLUID_INFO_VOLU_LIST]

    for i in range(len(objLists)):
        testTable = []
        testTable.append(header) # Headline

        for fluid in objLists[i]:
            testTable.append([fluid.name, fluid.description, fluid.reference])

        writeTableToFile(filLists[i], testTable)




