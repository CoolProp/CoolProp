#!/usr/bin/python
# -*- coding: ascii -*-
#
from __future__ import print_function, division
import os, sys
from os import path
import numpy as np
import CoolProp
import glob
from warnings import warn
from time import clock

import CoolProp.constants
from CoolProp.CoolProp import PropsSI, generate_update_pair, get_parameter_index, set_debug_level, get_phase_index
from CoolProp import AbstractState as State

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import copy
from itertools import cycle
from matplotlib import gridspec, ticker
#from jopy.dataPlotters import roundList, range_brace


def range_brace(x_min, x_max, mid=0.5,
                beta1=50.0, beta2=100.0, height=1,
                initial_divisions=11, resolution_factor=1.5):
    """
    http://stackoverflow.com/questions/1289681/drawing-braces-with-pyx
    x,y = range_brace(0, 100)
    ax.plot(x, y,'-')
    ax.plot(y, x,'-')
    """
    # determine x0 adaptively values using second derivative
    # could be replaced with less snazzy:
    #   x0 = NP.arange(0, 0.5, .001)
    x0 = np.array(())
    tmpx = np.linspace(0, 0.5, initial_divisions)
    tmp = beta1**2 * (np.exp(beta1 * tmpx)) * (1 - np.exp(beta1 * tmpx)) / np.power((1 + np.exp(beta1 * tmpx)), 3)
    tmp += beta2**2 * (np.exp(beta2 * (tmpx - 0.5))) * (1 - np.exp(beta2 * (tmpx - 0.5))) / np.power((1 + np.exp(beta2 * (tmpx - 0.5))), 3)
    for i in range(0, len(tmpx) - 1):
        t = int(np.ceil(resolution_factor * max(np.abs(tmp[i:i + 2])) / float(initial_divisions)))
        x0 = np.append(x0, np.linspace(tmpx[i], tmpx[i + 1], t))
    x0 = np.sort(np.unique(x0))  # sort and remove dups
    # half brace using sum of two logistic functions
    y0 = mid * 2 * ((1 / (1. + np.exp(-1 * beta1 * x0))) - 0.5)
    y0 += (1 - mid) * 2 * (1 / (1. + np.exp(-1 * beta2 * (x0 - 0.5))))
    # concat and scale x
    x = np.concatenate((x0, 1 - x0[::-1])) * float((x_max - x_min)) + x_min
    y = np.concatenate((y0, y0[::-1])) * float(height)
    return (x, y)

# try:
    #from jopy.dataPlotters import BasePlotter
    #bp = BasePlotter()
# except:
    #bp = None


bp = None


# The basic settings for he plots
xypoints = 1000
loops = 1
repeat = 1
runs = 0
maxruns = 5
plot = True
calc = True
check = True
folder = "dataTTSE"
figures = "figuresTTSE"
# np.random.seed(1984)


fluids = ["CO2", "Pentane", "R134a", "Water", "Air", "LiBr-0%"]
fluids = ["CO2", "Pentane", "R134a", "Water"]
fluids = ["Air"]

#glskeys = [r"\glsentryshort{co2}",r"\glsentryshort{pentane}",r"\glsentryshort{r134a}",r"\glsentryshort{water}",r"\glsentryshort{air}",r"\glsentryshort{libr} \SI{0}{\percent}"]
#glskeys = [r"\ce{CO2}",r"n-Ppentane",r"R134a",r"Water",r"Air",r"\glsentryshort{libr} \SI{0}{\percent}"]
repList = []
# for i in range(len(fluids)):
#    repList.append(fluids[i])
#    repList.append(glskeys[i])

#backends = ["INCOMP","HEOS","REFPROP"]
backends = ["HEOS", "REFPROP"]
backends = ["HEOS"]

# repList.append("HEOS")
# repList.append(r"\glsentryshort{cp}")
# repList.append("REFPROP")
# repList.append(r"\glsentryshort{rp}")


# CoolProp.CoolProp.set_debug_level(51)

pStr = path.dirname(path.abspath(__file__))
fStr = path.splitext(path.basename(__file__))[0]


def getFolderName():
    folderName = path.join(pStr, folder)
    if not path.isdir(folderName):
        print("Creating data directory " + folderName)
        os.makedirs(folderName)
    return folderName


def getFigureFolder():
    folderName = path.join(pStr, figures)
    if not path.isdir(folderName):
        print("Creating data directory " + folderName)
        os.makedirs(folderName)
    return folderName


repList.append("TimeComp-")
repList.append("chapters/FluidProperties/" + path.basename(getFigureFolder()) + "/TimeComp-")


def getFileName(qualifiers=[]):
    fileName = path.join(getFolderName(), "-".join(qualifiers))
    return fileName

# Some file handling


def loadNpzData(backend, fluid):
    dicts = {}
    globber = getFileName([backend, fluid]) + '_[0-9][0-9][0-9].npz'
    for fname in glob.glob(globber):
        dataDict = dict(np.load(fname))
        dicts[str(dataDict["name"])] = dataDict
    # if len(dicts)<1:
    #    #print("No readable file found for {0}".format(globber))
    #    dataDict = dict(name=str(0).zfill(3))
    #    dicts[str(dataDict["name"])] = dataDict
    return dicts


def saveNpzData(backend, fluid, dicts, start=0, stop=-1):
    keys = dicts.keys()
    keys.sort()
    for k in keys[start:stop]:
        data = dicts[k]
        fname = getFileName([backend, fluid]) + '_{0}.npz'.format(str(data['name']).zfill(3))
        np.savez(fname, **data)
    return True


def splitFluid(propsfluid):
    fld = propsfluid.split("::")
    if len(fld) == 2:
        backend = fld[0]
        fld = fld[1]
    else:
        backend = None
        fld = fld[0]
    fld = fld.split("-")
    if len(fld) == 2:
        conc = float(fld[1].strip('%')) / 100.0
        fld = fld[0]
    else:
        conc = None
        fld = fld[0]
    return backend, fld, conc


def getInpList(backend):
    if backend == "HEOS": return ["DT", "HP"]
    elif backend == "REFPROP": return ["DT", "HP"]
    elif backend == "INCOMP": return ["PT", "HP"]
    else: raise ValueError("Unknown backend.")


def getOutList(inp=None):
    if inp == "HP":
        return [["Tmax"], ["D"], ["S"], ["T"], ["D", "S", "T"]]
    elif inp == "DT":
        return [["Tmax"], ["H"], ["P"], ["S"], ["H", "P", "S"]]
    elif inp == "PT":
        return [["Tmax"], ["H"], ["D"], ["S"], ["H", "D", "S"]]
    else:
        raise ValueError("Unknown inputs.")

# def getPhaseString(iPhase):
#    for i in range(11):
#        if getPhaseString(i)==sPhase:
#            return i
#     if iPhase==1: return "liquid"
#     elif iPhase==2: return "supercritical"
#     elif iPhase==3: return "supercritical_gas"
#     elif iPhase==4: return "supercritical_liquid"
#     elif iPhase==5: return "critical_point"
#     elif iPhase==6: return "gas"
#     elif iPhase==7: return "twophase"
#     elif iPhase==8: return "unknown"
#     elif iPhase==9: return "not_imposed"
#     else: raise ValueError("Couldn't find phase.")


def getPhaseNum(sPhase):
    return get_phase_index(sPhase)
#     for i in range(11):
#         if getPhaseString(i)==sPhase:
#             return i


def getOutKey(out): return "".join(out)


def getOutLabel(out): return ",".join(out)


def getTimeKey(inp, out): return "_".join([inp, getOutKey(out)])


def getVectorKey(inp, out): return getTimeKey(inp, out) + "_V"


def getCriticalProps(propsfluid):
    backend, _, _ = splitFluid(propsfluid)
    if backend != "INCOMP":
        p_crit_m = PropsSI('pcrit', "T", 0, "D", 0, propsfluid) * 0.995
        T_crit_m = PropsSI('Tcrit', "T", 0, "D", 0, propsfluid) * 1.005
        d_crit_m = PropsSI('rhocrit', "T", 0, "D", 0, propsfluid) * 0.995
        h_crit_m = PropsSI('H', "T", T_crit_m, "D", d_crit_m, propsfluid)
        s_crit_m = PropsSI('H', "T", T_crit_m, "D", d_crit_m, propsfluid)
    else:
        p_crit_m = None
        T_crit_m = None
        d_crit_m = None
        h_crit_m = None
        s_crit_m = None
    return dict(P=p_crit_m, T=T_crit_m, D=d_crit_m, H=h_crit_m, S=s_crit_m)


def getPTRanges(propsfluid):
    backend, _, _ = splitFluid(propsfluid)
    # Setting the limits for enthalpy and pressure
    T_min = PropsSI('Tmin', "T", 0, "D", 0, propsfluid) + 1
    T_max = PropsSI('Tmax', "T", 0, "D", 0, propsfluid) - 1

    if backend == "REFPROP":
        T_min = max(T_min, PropsSI('Ttriple', "T", 0, "D", 0, propsfluid)) + 1
        p_min = PropsSI('P', "T", T_min, "Q", 0, propsfluid) + 1
        p_max = PropsSI('pmax', "T", 0, "D", 0, propsfluid) - 1
    elif backend == "INCOMP":
        p_min = 1.5 * 1e5
        p_max = 200.0 * 1e5
    else:
        T_min = max(T_min, PropsSI('Ttriple', "T", 0, "D", 0, propsfluid)) + 1
        p_min = PropsSI('ptriple', "T", 0, "D", 0, propsfluid)
        p_min = max(p_min, PropsSI('pmin', "T", 0, "D", 0, propsfluid)) + 1
        p_max = PropsSI('pmax', "T", 0, "D", 0, propsfluid) - 1

    # One more check to debug things:
    #p_min = max(p_min,0.01e5)
    #T_min = max(T_min,200)
    #p_max = min(p_max,200e5)
    #T_max = min(T_max,1750)

    p_range = np.logspace(np.log10(p_min), np.log10(p_max), xypoints)
    T_range = np.linspace(T_min, T_max, xypoints)
    return p_range, T_range

    #p_max = min(PropsSI('pcrit',"T",0,"D",0,fluid)*20, p_max)
    #T_max = min(PropsSI('Tcrit',"T",0,"D",0,fluid)* 3, T_max)


def getLists(propsfluid):
    backend, _, _ = splitFluid(propsfluid)
    """Returns randomised lists of all properties within the ranges"""
    p, T = getPTRanges(propsfluid)
    p_min = np.min(p)
    p_max = np.max(p)
    T_min = np.min(T)
    T_max = np.max(T)

    if backend == "INCOMP":
        h_min = PropsSI('H', 'T', T_min, 'P', p_min, propsfluid)
        h_max = PropsSI('H', 'T', T_max, 'P', p_max, propsfluid)
    else:
        critProps = getCriticalProps(propsfluid)
        h_min = PropsSI('H', 'T', T_min, 'Q', 0, propsfluid)
        h_max = PropsSI('H', 'T', T_min, 'Q', 1, propsfluid)
        h_max = max(PropsSI('H', 'T', critProps["T"], 'D', critProps["D"], propsfluid), h_max)
        h_max = (h_max - h_min) * 2.0 + h_min

    loop = True
    count = 0
    while loop:
        count += 1

        h_list = np.random.uniform(h_min, h_max, int(xypoints * 2.0))
        p_list = np.random.uniform(np.log10(p_min), np.log10(p_max), int(xypoints * 2.0))
        p_list = np.power(10, p_list)

        out = ["T", "D", "S"]
        res = PropsSI(out, "P", p_list, "H", h_list, propsfluid)
        T_list = res[:, 0]
        d_list = res[:, 1]
        s_list = res[:, 2]

        mask = np.isfinite(T_list) & np.isfinite(d_list) & np.isfinite(s_list)
        if np.sum(mask) < xypoints:
            if False:
                print(h_list); print(p_list); print(T_list); print(d_list); print(s_list)
            print("There were not enough valid entries in your result vector: {0:d} > {1:d} - rerunning".format(xypoints, np.sum(mask)))
            loop = True
        else:
            loop = False
            p_list = p_list[mask][0:xypoints]
            h_list = h_list[mask][0:xypoints]
            T_list = T_list[mask][0:xypoints]
            d_list = d_list[mask][0:xypoints]
            s_list = s_list[mask][0:xypoints]
            return dict(P=p_list, T=T_list, D=d_list, H=h_list, S=s_list)

        maxTries = 4
        if count > maxTries:
            loop = False
            raise ValueError("{0}: Could not fill the lists in {0} runs, aborting.".format(propsfluid, maxTries))


def getInpValues(inp, dataDict):
    in1 = inp[0]
    in2 = dataDict[in1]
    in3 = inp[1]
    in4 = dataDict[in3]
    return in1, in2, in3, in4


def getStateObj(propsfluid):
    backend, fld, conc = splitFluid(propsfluid)
    # fluidstr holds the full information and fluid is only the name
    # Initialise the state object
    if backend is not None:
        state = State(backend, fld)
    else:
        state = State(fld)
    # if backend=="INCOMP":
    #    state.set_mass_fractions([0.0])
    if conc is not None:
        try:
            state.set_mass_fractions([conc])
        except:
            pass
    return state


def getSpeedMeas(out, in1, in2, in3, in4, propsfluid, vector=False):
    pair, out1, _ = generate_update_pair(get_parameter_index(in1), in2[0], get_parameter_index(in3), in4[0])
    if out1 == in2[0]: swap = False
    else: swap = True
    if swap:
        input1 = in4
        input2 = in2
    else:
        input1 = in2
        input2 = in4
    state = getStateObj(propsfluid)
    outList = [get_parameter_index(i) for i in out]
    resLst = np.empty((repeat,))
    timLst = np.empty((repeat,))
    if vector:
        for j in range(repeat):
            timLst.fill(np.inf)
            lrange = range(len(input1))
            resTmp = np.inf
            if len(outList) == 1 and outList[0] == CoolProp.constants.iT_max:
                t1 = clock()
                for l in lrange:
                    for o in outList:
                        resTmp = state.keyed_output(o)
                t2 = clock()
                timLst[j] = (t2 - t1) * 1e6 / float(len(input1))
            else:  # We have to update before doing other things
                t1 = clock()
                for l in lrange:
                    state.update(pair, input1[l], input2[l])
                    for o in outList:
                        resTmp = state.keyed_output(o)
                t2 = clock()
                timLst[j] = (t2 - t1) * 1e6 / float(len(input1))
        res = None
        tim = np.min(timLst)  # Best of (repeat)
        return res, tim
    else:
        res = np.empty_like(input1)
        res.fill(np.inf)
        tim = np.empty_like(input1)
        tim.fill(np.inf)
        for i in range(len(input1)):
            resLst.fill(np.inf)
            timLst.fill(np.inf)
            for j in range(repeat):
                lrange = range(loops)
                resTmp = np.inf
                if len(outList) == 1 and outList[0] == CoolProp.constants.iT_max:
                    t1 = clock()
                    for _ in lrange:
                        for o in outList:
                            resTmp = state.keyed_output(o)
                    t2 = clock()
                    timLst[j] = (t2 - t1) * 1e6 / float(loops)
                    resLst[j] = resTmp
                else:  # We have to update before doing other things
                    inV1 = input1[i]
                    inV2 = input2[i]  # *(1.0+(l/1000.0)*pow(-1,l)) for l in lrange ]
                    t1 = clock()
                    for l in lrange:
                        state.update(pair, inV1, inV2)
                        for o in outList:
                            resTmp = state.keyed_output(o)
                    t2 = clock()
                    timLst[j] = (t2 - t1) * 1e6 / float(loops)
                    resLst[j] = resTmp
            if not np.all(resLst == resLst[0]):
                raise ValueError("Not all results were the same.")
            res[i] = resLst[0]
            tim[i] = np.min(timLst)  # Best of three (repeat)
        return res, tim


def checkDataSet(propsfluid, dataDict, fill=True, quiet=False):
    if not check: return
    backend, _, _ = splitFluid(propsfluid)
    if not quiet: print("\n\n-- {0:16s} --".format(propsfluid), end="")
    # Test for required inputs
    newInputs = False
    inLists = getLists(propsfluid)
    for inp in getInpList(backend):
        if not quiet: print("\n{0:2s}: ".format(inp), end="")
        for inVal in inp:
            if inVal not in dataDict:  # A problem
                if not fill:
                    raise ValueError("The input {0:1s} is missing or faulty, cannot continue.".format(inVal))
                #dataDict[inVal] = inLists[inVal]
                dataDict.update(inLists)
                newInputs = True
                if not quiet: print("{0:s}*({1:d}),".format(inVal, len(dataDict[inVal])), end="")
            else:
                if not quiet: print("{0:s} ({1:d}),".format(inVal, len(dataDict[inVal])), end="")
        # All inputs are there
        in1, in2, in3, in4 = getInpValues(inp, dataDict)
        #in2 = in2[:3]
        #in4 = in4[:3]
        if in2.shape != in4.shape:
            raise ValueError("The stored data for {0:s} and {1:s} do not have the same shape.".format(in1, in3))
        if in2.shape != inLists[in1].shape:
            raise ValueError("The stored data for {0:s} and its list do not have the same shape {1} vs {2}.".format(in1, in2.shape, inLists[in1].shape))
        # Check for time data
        for out in getOutList(inp):
            key = getTimeKey(inp, out)
            okey = getOutKey(out)
            if key not in dataDict or newInputs or not np.all(np.isfinite(dataDict[key])):
                if not fill:
                    raise ValueError("The time data for {0:s} is missing or faulty, cannot continue.".format(key))
                res, tim = getSpeedMeas(out, in1, in2, in3, in4, propsfluid)
                dataDict[key] = tim
                dataDict[okey] = res  # We calculated in, why not use it here...
                if not quiet: print("{0:s}*({1:d}),".format(key, len(dataDict[key])), end="")
            else:
                if not quiet: print("{0:s} ({1:d}),".format(key, len(dataDict[key])), end="")
            if dataDict[key].shape != in2.shape or not np.all(np.isfinite(dataDict[key])):
                raise ValueError("The stored time data for {0:s} does not have the same shape as the inputs.".format(key))
        # Check for vectors
        for out in getOutList(inp):
            key = getVectorKey(inp, out)
            if key not in dataDict or not np.all(np.isfinite(dataDict[key])):
                if not fill:
                    raise ValueError("The fluid data for {0:s} is missing or faulty, cannot continue.".format(key))
                res, tim = getSpeedMeas(out, in1, in2, in3, in4, propsfluid, vector=True)
                dataDict[key] = tim
                if not quiet: print("{0:s}*({1:d}),".format(key, dataDict[key].size), end="")
            else:
                if not quiet: print("{0:s} ({1:d}),".format(key, dataDict[key].size), end="")
            if dataDict[key].size != 1 or not np.all(np.isfinite(dataDict[key])):
                raise ValueError("The vector data for {0:s} does not have the correct size {1}..".format(key, dataDict[key].size))

#     inp = getInpList(backend)[0] # Explicit calls
#      # Check for properties
#      for out in getOutList(inp)[:-1]:
#          key = getOutKey(out)
#          if key not in dataDict or not np.all(np.isfinite(dataDict[key])):
#              if not fill:
#                  raise ValueError("The fluid data for {0:s} is missing or faulty, cannot continue.".format(key))
#              res = PropsSI(out,in1,in2,in3,in4,propsfluid)
#              dataDict[key] = res
#              if not quiet: print("{0:s}*({1:d}),".format(key,len(dataDict[key])),end="")
#          else:
#              if not quiet: print("{0:s} ({1:d}),".format(key,len(dataDict[key])),end="")
#          if dataDict[key].shape != in2.shape or not np.all(np.isfinite(dataDict[key])):
#              raise ValueError("The stored data for {0:s} does not have the same shape as the inputs {1} vs {2}..".format(key,dataDict[key].shape,in2.shape))
#     # Check for phase
#     for out in ["Phase"]:
#         if backend!="HEOS":
#             dataDict[key] = np.zeros_like(a, dtype, order, subok)
#         key = getOutKey(out)
#         if key not in dataDict or newInputs or not np.all(np.isfinite(dataDict[key])):
#             if not fill:
#                 raise ValueError("The phase data for {0:s} is missing or faulty, cannot continue.".format(key))
#             res = np.empty_like(in2)
#             res.fill(np.inf)
#             for i in range(len(in2)):
#                 res[i] = PropsSI(out,in1,in2[i],in3,in4[i],propsfluid)
#             dataDict[key] = res
#             if not quiet: print("{0:s}*({1:d}),".format(key,len(dataDict[key])),end="")
#         else:
#             if not quiet: print("{0:s} ({1:d}),".format(key,len(dataDict[key])),end="")
#         if dataDict[key].shape != in2.shape or not np.all(np.isfinite(dataDict[key])):
#             raise ValueError("The stored data for {0:s} does not have the same shape as the inputs {1} vs {2}..".format(key,dataDict[key].shape,in2.shape))
#
#             # Now we use the vector data
#             key = getVectorKey(inp, out)
#             if key not in dataDict or not np.all(np.isfinite(dataDict[key])):
#                 if not fill:
#                     raise ValueError("The vector data for {0:s} is missing or faulty, cannot continue.".format(key))
#                 dataDict[key] = np.empty_like(in2)
#                 dataDict[key].fill(np.inf)
#                 res = []
#                 for _ in range(repeat):
#                     t1=clock()
#                     PropsSI(out,in1,in2,in3,in4,propsfluid)
#                     t2=clock()
#                     res.append((t2-t1)/float(len(in2)))
#                 dataDict[key] = np.min(res)
#                 if not quiet: print("{0:s}*({1}),".format(key,dataDict[key]),end="")
#             else:
#                 if not quiet: print("{0:s} ({1}),".format(key,dataDict[key]),end="")
#             try:
#                 float(dataDict[key])
#             except:
#                 raise ValueError("The stored vector data for {0:s} cannot be casted to float.".format(key))
#        if not quiet: print("")

# All data is loaded and checked, we can calculate more now


def getEntryCount(dicts, backend, fld):
    return len(fluidData[fld][backend].keys())


def getUKey(fld, bck, inp, out):
    return "-".join([fld, bck, inp, "".join(out)])


def getData(fld, backend, inp, out, fluidData):
    inputs1 = []
    inputs2 = []
    values = []
    times = []
    i1key = inp[0]
    i2key = inp[1]
    vkey = getOutKey(out)
    tkey = getTimeKey(inp, out)
    for dkey in fluidData[fld][backend]:
        cData = fluidData[fld][backend][dkey]
        inputs1.append(cData[i1key])
        inputs2.append(cData[i2key])
        values.append(cData[vkey])
        times.append(cData[tkey])
    ret = {}
    if len(inputs1) > 0:
        ret[i1key] = np.concatenate(inputs1)
        ret[i2key] = np.concatenate(inputs2)
        ret[vkey] = np.concatenate(values)
        ret[tkey] = np.concatenate(times)
    return ret


def getSingleData(fld, backend, key, fluidData):
    #print("Getting: "+fld+", "+backend+", "+key)
    values = []
    for dkey in fluidData[fld][backend]:
        if key in fluidData[fld][backend][dkey]:
            if "P" in fluidData[fld][backend][dkey]:
                # TODO: Fix this, do we need the mask?
                #mask = fluidData[fld][backend][dkey]["P"]>0.3e5
                mask = fluidData[fld][backend][dkey]["P"] > 0.0e5
                try:
                    values.append(fluidData[fld][backend][dkey][key][mask])
                except Exception as e:
                    values.append(fluidData[fld][backend][dkey][key])
                    print(e)
                    pass
            else:
                values.append(fluidData[fld][backend][dkey][key])
    if len(values) > 0:
        if np.size(values[0]) > 1:
            return np.concatenate(values)
        else:
            return np.array(values)
    return None


def fillDict(fld, backend, fluidData, curDict, curKeys):
    if curDict is None: curDict = {}
    for key in curKeys:
        vals = getSingleData(fld, backend, key, fluidData)
        if vals is not None: curDict[key] = vals
    return curDict

################################################


fluidData = {}

for fluidstr in fluids:
    _, fld, _ = splitFluid(fluidstr)
    if fld not in fluidData: fluidData[fld] = {}
    for backend in backends:
        if backend not in fluidData[fld]:  # Try to add it
            propsfluid = "::".join([backend, fluidstr])
            try:
                PropsSI('Tmax', "T", 0, "D", 0, propsfluid)
                fluidData[fld][backend] = loadNpzData(backend, fld)
            except:
                pass
            if backend in fluidData[fld]:
                for k in fluidData[fld][backend]:
                    checkDataSet(propsfluid, fluidData[fld][backend][k], fill=False, quiet=True)
        else:  # Already added backend for fluid
            pass

lastSave = 0
while runs < maxruns and calc:
    check = True  # force checking for new records
    runs += 1
    # Now we have the data structure with the precalculated data
    for fluidstr in fluids:
        _, fld, _ = splitFluid(fluidstr)
        for backend in fluidData[fld]:
            propsfluid = "::".join([backend, fluidstr])
            dicts = fluidData[fld][backend]
            keys = list(dicts.keys())
            keys.sort()
            while len(keys) < runs:
                if len(keys) < 1: newKey = 0
                else: newKey = int(keys[-1]) + 1
                if newKey > 999:
                    raise ValueError("Too many dicts: {0}>999".format(newKey))
                k = str(newKey).zfill(3)
                dicts[k] = {}
                dicts[k]['name'] = k
                try:
                    checkDataSet(propsfluid, dicts[k], fill=True, quiet=False)
                except Exception as e:
                    print("There was an error, dataset {0} from {1}.might be faulty:\n{2}".format(k, propsfluid, str(e)))
                    pass
                todel = []
                # for k in dicts:
                try:
                    checkDataSet(propsfluid, dicts[k], fill=False, quiet=True)
                except Exception as e:
                    print("There was an error, removing dataset {0} from {1}.\n{2}".format(k, propsfluid, str(e)))
                    todel.append(k)
                for t in todel: del dicts[t]
                keys = list(dicts.keys())
                keys.sort()

    # Updated all dicts for this backend, saving data
    if runs >= maxruns or (lastSave + 4) < runs:
        print("\n\nDone processing fluids, saving data: ")
        for fld in fluidData:
            for backend in fluidData[fld]:
                saveNpzData(backend, fld, fluidData[fld][backend], start=lastSave, stop=runs)
                print("{0} ({1})".format(backend + "::" + fld, len(fluidData[fld][backend].keys()[lastSave:runs])), end=", ")
        print("")
        lastSave = runs


if not plot: sys.exit(0)

# exclLst = [["Tmax"],["H"],["D"],["S"],["H","D","S"],["D"],["S"],["T"],["D","S","T"]]
#
# Start with a temporary dictionary that holds all the data we need
# for fld in fluidData:
#     cstData = {} # Data from calls to constants (overhead)
#     expData = {} # Data from explicit EOS calls
#     impData = {} # Data from EOS calls that require iterations
#     for backend in fluidData[fld]:
#
#         curKeys = []
#         for inp in getInpList(backend):
#             for out in getOutList(inp)[:1]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         cstData[backend] = curDict
#
#         curKeys = []
#         for inp in getInpList(backend)[:1]:
#             for out in getOutList(inp)[1:]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         expData[backend] = curDict
#
#         curKeys = []
#         for inp in getInpList(backend)[1:]:
#             for out in getOutList(inp)[1:]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         impData[backend] = curDict


#         curDict[backend] = {}
#         for key in :
#             vals = getSingleData(fld, backend, key, fluidData)
#             if vals is not None: curDict[backend][key] = vals
#             if curDict
#
#         cstData[backend] = {}
#         for out in getOutList(inp[0])[0]:
#             res = getData(fld,backend,inp[0],out,fluidData)
#             cstData[backend].update(res)
#         for out in getOutList(inp[1])[0]:
#             res = getData(fld,backend,inp[1],out,fluidData)
#             cstData[backend].update
#
#         expData[backend] = {}
#         for out in getOutList(inp[0])[1:]:
#             res = getData(fld,backend,inp[0],out,fluidData)
#             expData[backend].update(res)
#
#         impData[backend] = {}
#         for out in getOutList(inp[1])[1:]:
#             res = getData(fld,backend,inp[1],out,fluidData)
#             impData[backend].update(res)

    #############################################################
    # All data is available in the dicts now.
    #############################################################
    # The first thuing to do is to print some statistical
    # measures to give you an idea about the data.
#     try:
#         #dataOHCP = [cstData["HEOS"]["DT_Tmax"]   , cstData["HEOS"]["HP_Tmax"]   ]
#         #dataOHRP = [cstData["REFPROP"]["DT_Tmax"], cstData["REFPROP"]["HP_Tmax"]]
#         print("\n{0} - {1} points ".format(fld,np.size(cstData["HEOS"]["DT_Tmax"])))
#         print("Overhead CoolProp: {0:5.3f} us".format(np.mean(cstData["HEOS"]["DT_Tmax"])))#,   np.mean(cstData["HEOS"]["HP_Tmax"]))
#         print("Overhead REFPROP : {0:5.3f} us".format(np.mean(cstData["REFPROP"]["DT_Tmax"])))#,np.mean(cstData["REFPROP"]["HP_Tmax"]))
#         print("Mean EOS in  CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.mean(expData["HEOS"]["DT_HPS"]),np.mean(impData["HEOS"]["HP_DST"])))
#         print("Std. dev. in CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.std(expData["HEOS"]["DT_HPS"]) ,np.std(impData["HEOS"]["HP_DST"])))
#         print("Minimum in   CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.min(expData["HEOS"]["DT_HPS"]) ,np.min(impData["HEOS"]["HP_DST"])))
#         print("Maximum in   CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.max(expData["HEOS"]["DT_HPS"]) ,np.max(impData["HEOS"]["HP_DST"])))
#         print("Mean EOS in   REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.mean(expData["REFPROP"]["DT_HPS"]),np.mean(impData["REFPROP"]["HP_DST"])))
#         print("Std. dev. in  REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.std(expData["REFPROP"]["DT_HPS"]) ,np.std(impData["REFPROP"]["HP_DST"])))
#         print("Minimum in    REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.min(expData["REFPROP"]["DT_HPS"]) ,np.min(impData["REFPROP"]["HP_DST"])))
#         print("Maximum in    REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.max(expData["REFPROP"]["DT_HPS"]) ,np.max(impData["REFPROP"]["HP_DST"])))
#         print("")
#
#         print("\n{0} - {1} points ".format(fld,np.size(cstData["HEOS"]["DT_Tmax_V"])))
#         print("Overhead CoolProp: {0:5.3f} us".format(np.mean(cstData["HEOS"]["DT_Tmax_V"])))#,   np.mean(cstData["HEOS"]["HP_Tmax"]))
#         print("Overhead REFPROP : {0:5.3f} us".format(np.mean(cstData["REFPROP"]["DT_Tmax_V"])))#,np.mean(cstData["REFPROP"]["HP_Tmax"]))
#         print("Mean EOS in  CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.mean(expData["HEOS"]["DT_HPS_V"]),np.mean(impData["HEOS"]["HP_DST_V"])))
#         print("Std. dev. in CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.std(expData["HEOS"]["DT_HPS_V"]) ,np.std(impData["HEOS"]["HP_DST_V"])))
#         print("Minimum in   CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.min(expData["HEOS"]["DT_HPS_V"]) ,np.min(impData["HEOS"]["HP_DST_V"])))
#         print("Maximum in   CoolProp: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.max(expData["HEOS"]["DT_HPS_V"]) ,np.max(impData["HEOS"]["HP_DST_V"])))
#         print("Mean EOS in   REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.mean(expData["REFPROP"]["DT_HPS_V"]),np.mean(impData["REFPROP"]["HP_DST_V"])))
#         print("Std. dev. in  REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.std(expData["REFPROP"]["DT_HPS_V"]) ,np.std(impData["REFPROP"]["HP_DST_V"])))
#         print("Minimum in    REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.min(expData["REFPROP"]["DT_HPS_V"]) ,np.min(impData["REFPROP"]["HP_DST_V"])))
#         print("Maximum in    REFPROP: f(rho,T): {0:9.3f} us, f(h,p): {1:9.3f} us".format(np.max(expData["REFPROP"]["DT_HPS_V"]) ,np.max(impData["REFPROP"]["HP_DST_V"])))
#         print("")
#
#     except Exception as e:
#         print(str(e.message))
#         pass
#
# try:
#     fld = 'Water'
#     cstData = {} # Data from calls to constants (overhead)
#     expData = {} # Data from explicit EOS calls
#     impData = {} # Data from EOS calls that require iterations
#     for backend in fluidData[fld]:
#         curKeys = []
#         for inp in getInpList(backend):
#             for out in getOutList(inp)[:1]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         cstData[backend] = curDict
#
#         curKeys = []
#         for inp in getInpList(backend)[:1]:
#             for out in getOutList(inp)[1:]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         expData[backend] = curDict
#
#         curKeys = []
#         for inp in getInpList(backend)[1:]:
#             for out in getOutList(inp)[1:]:
#                 curKeys.append(getTimeKey(  inp, out))
#                 curKeys.append(getVectorKey(inp, out))
#         curDict = {}
#         fillDict(fld,backend,fluidData,curDict,curKeys)
#         impData[backend] = curDict
#
#
#     print("Done")

def autolabel(ax, rects, means, stds, lens=0, fac=1):
    return
    # attach some text labels
    yerr = (stds * 100.0) / means
    ypos = np.max(means)  # + np.max(stds)
    for i in range(len(rects)):
        xpos = rects[i].get_x() + rects[i].get_width() / 2.
        #ax.text(xpos, 1.05*ypos, '{0:s}{1:4.2f}{2:s}{3:3.1f}{4:s}'.format(r'',means[i]*fac,r'us +- ',yerr[i],r'%'), rotation=90, ha='center', va='bottom', fontsize='smaller')
        ax.text(xpos, 1.05 * ypos, '{0:s}{1:4.2f}{2:s}{3:3.1f}{4:s}'.format(r'\SI{', means[i] * fac, r'}{\us} (\SI{\pm', yerr[i], r'}{\percent})'), rotation=90, ha='center', va='bottom', fontsize='smaller')
        #ax.text(xpos, 0.25*ypos, str(lens[i]), rotation=90, ha='center', va='bottom', fontsize='smaller')


def axislabel(txt, ax, xPos, yPos=-1):
    ax.text(xPos, yPos, txt, rotation=45, ha='center', va='top', fontsize='xx-small')


#############################################################
# All data is available in the dicts now.
#############################################################
# The first plot contains the time data, this is averaged and
# plotted as a bar graph with the standard deviation.
hatchLst = ["", "///", "\\\\\\"]
backendsLst = ["INCOMP", "HEOS", "REFPROP"]
for fluidstr in fluids[:-1]:
    _, fld, _ = splitFluid(fluidstr)
    outCst = []
    labCst = []
    outExp = []
    labExp = []
    outImp = []
    labImp = []
    DEBUG = True
    for backend in backendsLst:
        # Backend exists in fluid data?
        try:
            for i, inp in enumerate(getInpList(backend)):
                for j, out in enumerate(getOutList(inp)):
                    if j == 0:  # First output is Tmax
                        if backend not in fluidData[fld]:
                            outCst.append([0])
                            labCst.append("Dummy")
                            if DEBUG: print("Added a dummy for {0} and {1},{2},{3}".format(fld, backend, inp, out))
                        else:
                            outCst.append(getSingleData(fld, backend, getTimeKey(inp, out), fluidData))
                            labCst.append(getTimeKey(inp, out))
                        continue
                    elif i == 0:  # First input is explicit
                        if backend not in fluidData[fld]:
                            outExp.append([0])
                            labExp.append("Dummy")
                            if DEBUG: print("Added a dummy for {0} and {1},{2},{3}".format(fld, backend, inp, out))
                        else:
                            outExp.append(getSingleData(fld, backend, getTimeKey(inp, out), fluidData))
                            labExp.append(getTimeKey(inp, out))
                        continue
                    elif i == 1:
                        if backend not in fluidData[fld]:
                            outImp.append([0])
                            labImp.append("Dummy")
                            if DEBUG: print("Added a dummy for {0} and {1},{2},{3}".format(fld, backend, inp, out))
                        else:
                            outImp.append(getSingleData(fld, backend, getTimeKey(inp, out), fluidData))
                            labImp.append(getTimeKey(inp, out))
                        continue
                    else:
                        raise ValueError("Wrong number.")
        except Exception as e:
            print(e)
            sys.exit(1)

    # Do the plotting
    if bp is not None:
        bp.figure = None
        fg = bp.getFigure()
        ccycle = bp.getColorCycle(length=3)
    else:
        fg = plt.figure()
        ccycle = cycle(["b", "g", "r"])

    #fg = plt.figure()
    ax1 = fg.add_subplot(111)
    ax2 = ax1.twinx()
    rects1 = []
    labels1 = []
    rects2 = []
    labels2 = []
    rects3 = []
    labels3 = []
    col1 = ccycle.next()
    col2 = ccycle.next()
    col3 = ccycle.next()

    numBackends = len(backendsLst)

    step = 1
    width = step / (numBackends + 1)  # the width of the bars
    offset = -0.5 * numBackends * width

    entries = int((len(outCst) + len(outExp) + len(outImp)) / numBackends)

#     for o in range(entries):
#         ids = np.empty((numBackends,))
#         for b in range(numBackends):
#             i = o*numBackends+b
#             ids[b] = offset + o*step + b*width
#             j = i - 0
#             if j < len(outCst):
#                 rects1.extend(ax1.bar(ids[b], np.mean(outCst[j]), width, color=col1, hatch=hatchLst[b]))#, yerr=np.std(curList[i]), ecolor='k'))
#
#
#             j = i - 0
#             if j < len(outCst):
#                 rects1.extend(ax1.bar(ids[b], np.mean(outCst[j]), width, color=col1, hatch=hatchLst[b]))#, yerr=np.std(curList[i]), ecolor='k'))
#             else:
#                 j = i-len(outCst)
#                 if j < len(outExp):
#                     rects2.extend(ax1.bar(ids[b], np.mean(outExp[j]), width, color=col2, hatch=hatchLst[b]))#, yerr=np.std(curList[i]), ecolor='k'))
#                 else:
#                     j = i-len(outCst)-len(outExp)
#                     if j < len(outImp):
#                         rects3.extend(ax2.bar(ids[b], np.mean(outImp[j]), width, color=col3, hatch=hatchLst[b]))#, yerr=np.std(curList[i]), ecolor='k'))
#                     else:
#                         raise ValueError("Do not go here!")
    DEBUG = True

    entries = 2
    for o in range(entries):
        ids = np.empty((numBackends,))
        for b in range(numBackends):
            i = b * entries + o
            try:
                ids[b] = offset + o * step + b * width
                rects1.extend(ax1.bar(ids[b], np.mean(outCst[i]), width, color=col1, hatch=hatchLst[b], rasterized=False))  # , yerr=np.std(curList[i]), ecolor='k'))
                if DEBUG:
                    print("Plotting {0}: {1:7.1f} - {2} - {3}".format(fld, np.mean(outCst[i]), "cst.", b))
                    # print(ids[b],labCst[i])
                    #axislabel(labCst[i], ax1, ids[b]+0.5*width)
            except:
                pass

    offset += entries * step
    #entries = int(len(outExp)/numBackends)
    entries = 4
    for o in range(entries):
        ids = np.empty((numBackends,))
        for b in range(numBackends):
            i = b * entries + o
            try:
                ids[b] = offset + o * step + b * width
                rects2.extend(ax1.bar(ids[b], np.mean(outExp[i]), width, color=col2, hatch=hatchLst[b], rasterized=False))  # , yerr=np.std(curList[i]), ecolor='k'))
                if DEBUG:
                    print("Plotting {0}: {1:7.1f} - {2} - {3}".format(fld, np.mean(outExp[i]), "exp.", b))
                    # print(ids[b],labExp[i])
                    #axislabel(labExp[i], ax1, ids[b]+0.5*width)
            except:
                pass

    x_newaxis = np.max(ids) + 1.5 * width
    plt.axvline(x_newaxis, color='k', linestyle='dashed')

    offset += entries * step
    entries = 4
    for o in range(entries):
        ids = np.empty((numBackends,))
        for b in range(numBackends):
            i = b * entries + o
            try:
                ids[b] = offset + o * step + b * width
                rects3.extend(ax2.bar(ids[b], np.mean(outImp[i]), width, color=col3, hatch=hatchLst[b], rasterized=False))  # , yerr=np.std(curList[i]), ecolor='k'))
                if DEBUG:
                    print("Plotting {0}: {1:7.1f} - {2} - {3}".format(fld, np.mean(outImp[i]), "imp.", b))
                    # print(ids[b],labImp[i])
                    #axislabel(labImp[i], ax1, ids[b]+0.5*width)
            except:
                pass

    # ax1.set_xlim([ids.min()-2.5*width,ids.max()+2.5*width])
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_ticks_position('bottom')

    labels = [r"ex.", r"im.", r"$h$", r"$\rho|p$", r"$s$", r"all", r"$\rho$", r"$s$", r"$T$", r"all"]
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels)
    # ax1.yaxis.get_label().set_verticalalignment("baseline")

    x_min = rects1[0].get_x()
    dx = rects1[0].get_width()
    x_max = rects3[-1].get_x()

    x_min = x_min - 1 * dx
    x_max = x_max + 2 * dx

    ax1.set_xlim([x_min, x_max])

    y_min = 0

    y_max_c = np.nanmax([a.get_height() for a in rects1])
    y_max_e = np.nanmax([a.get_height() for a in rects2])
    y_max_i = np.nanmax([a.get_height() for a in rects3])

    y_max = np.max([y_max_c, y_max_e, y_max_i / 10.0])

    y_max = np.ceil(1.3 * y_max / 10.0) * 10.0

    ax1.set_ylim([y_min, y_max])
    ax2.set_ylim([y_min, y_max * 10.0])

    ratio = 10.0 / 4.0 * y_max / 250.0  # height of 10 for 4 points if y_max==250

    x_min = rects1[0].get_x()
    x_max = rects1[-1].get_x() + dx
    x, y = range_brace(x_min, x_max)
    dy = np.ceil(y_max_c / 10.0) * 10.0
    y = dy + y * ratio * (x[-1] - x[0])
    ax1.plot(x, y, ls='-', color='k')
    ax1.text(np.mean(x), np.max(y), "const.", rotation=0, ha='center', va='bottom', fontsize='medium')

    x_min = rects2[0].get_x()
    x_max = rects2[-1].get_x() + dx
    x, y = range_brace(x_min, x_max)
    dy = np.ceil(y_max_e / 10.0) * 10.0
    y = dy + y * ratio * (x[-1] - x[0])
    ax1.plot(x, y, ls='-', color='k')
    ax1.text(np.mean(x), np.max(y), "explicit", rotation=0, ha='center', va='bottom', fontsize='medium')

    x_min = rects3[0].get_x()
    x_max = rects3[-1].get_x() + dx
    x, y = range_brace(x_min, x_max)
    dy = np.ceil(y_max_i / 100.0) * 10
    y = dy + y * ratio * (x[-1] - x[0])
    ax1.plot(x, y, ls='-', color='k')
    ax1.text(np.mean(x), np.max(y), "implicit", rotation=0, ha='center', va='bottom', fontsize='medium')

    #ax1.text(x_newaxis*0.9, y_max*0.9, "<- left axis", rotation=0, ha='right', va='bottom', fontsize='medium')
    #ax1.text(x_newaxis*1.1, y_max*0.9, "right axis ->", rotation=0, ha='left', va='bottom', fontsize='medium')

    handles = []
    for h in (rects1[0], rects2[1], rects3[2]):
        handles.append(copy.copy(h))
        handles[-1].set_facecolor('white')
        handles.append(copy.copy(h))
        handles[-1].set_hatch('')
    labels = (r'$p,T$-fit', r'constant', r'CoolProp', r'explicit, $f(p|\rho,T)$', r'REFPROP', r'implicit, $f(h,p)$')

    if bp is not None:
        bp.drawLegend(ax=ax1,
          loc='lower center',
          bbox_to_anchor=(0.5, 1.05),
          ncol=3,
          handles=handles,
          labels=labels)
    else:
        ax1.legend(handles, labels,
          loc='lower center',
          bbox_to_anchor=(0.5, 1.05),
          ncol=3)

    ax1.set_ylabel(r'Time per explicit call (us)')
    ax2.set_ylabel(r'Time per implicit call (us)')
    fg.savefig(path.join(getFigureFolder(), "TimeComp-" + fld.lower() + ".pdf"))
    if bp is not None:
        ax1.set_ylabel(r'Time per explicit call (\si{\us})')
        ax2.set_ylabel(r'Time per implicit call (\si{\us})')
        mpl.rcParams['text.usetex'] = True
        fg.savefig(path.join(getFigureFolder(), "TimeComp-" + fld.lower() + "-tex.pdf"))
        mpl.rcParams['text.usetex'] = False
        # Fix the wrong baseline
        for tick in ax1.get_xaxis().get_major_ticks():
            tick.set_pad(2 * tick.get_pad())
            tick.label1 = tick._get_text1()
        for lab in ax1.xaxis.get_ticklabels():
            lab.set_verticalalignment("baseline")
            # lab.set_pad(1.5*lab.get_pad())

        # ax1.set_xticklabels(labels)
        #
        # for tick in ax1.xaxis.get_major_ticks():
        #    tick.label1.set_horizontalalignment('center')
        bp.savepgf(path.join(getFigureFolder(), "TimeComp-" + fld.lower() + ".pgf"), fg, repList)

    plt.close()

# for fld in fluids:
#     try:
#         if bp is not None:
#             bp.figure = None
#             fg = bp.getFigure()
#             ccycle = bp.getColorCycle(length=3)
#         else:
#             fg = plt.figure()
#             ccycle = cycle(["b","g","r"])
#
#         #fg = plt.figure()
#         ax1 = fg.add_subplot(111)
#         ax2 = ax1.twinx()
#
#         if "INCOMP" in fluidData[fld]:
#             el = 3
#             hatch = ["","//","x"]
#         else:
#             el = 2
#             hatch = ["//","x"]
#
#         #one standard deviation above and below the mean of the data
#         width  = 0.25        # the width of the bars
#         step   = 1
#         offset = -step-0.5*el*width
#
#         lab   = []
#         rects1 = []
#         rects2 = []
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "PT_Tmax", fluidData)*1e2)
#         curList.append(getSingleData(fld, "HEOS"   , "DT_Tmax", fluidData)*1e2)
#         curList.append(getSingleData(fld, "REFPROP", "DT_Tmax", fluidData)*1e2)
#
#         lab.extend(["Tmax"])
#
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         curCol = ccycle.next()
#         for i in range(len(hatch)):
#             rects1.extend(ax1.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax1,rects1[0:],np.mean(curList,axis=1),np.std(curList,axis=1),fac=1e-2)
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "PT_H", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "DT_H", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "DT_H", fluidData))
#         lab.extend(["H"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         curCol = ccycle.next()
#         for i in range(len(hatch)):
#             rects1.extend(ax1.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax1,rects1[el:],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "PT_D", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "DT_P", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "DT_P", fluidData))
#         if el==3: lab.extend(["D/P"])
#         else: lab.extend(["P"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "g"
#         for i in range(len(hatch)):
#             rects1.extend(ax1.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax1,rects1[int(2*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "PT_S", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "DT_S", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "DT_S", fluidData))
#         lab.extend(["S"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "g"
#         for i in range(len(hatch)):
#             rects1.extend(ax1.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax1,rects1[int(3*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "PT_HDS", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "DT_HPS", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "DT_HPS", fluidData))
#         lab.extend(["all"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "g"
#         for i in range(len(hatch)):
#             rects1.extend(ax1.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax1,rects1[int(4*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "HP_D", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "HP_D", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "HP_D", fluidData))
#         lab.extend(["D"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         curCol = ccycle.next()
#         for i in range(len(hatch)):
#             rects2.extend(ax2.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax2,rects2[0:],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "HP_T", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "HP_T", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "HP_T", fluidData))
#         lab.extend(["T"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "r"
#         for i in range(len(hatch)):
#             rects2.extend(ax2.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax2,rects2[int(1*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "HP_S", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "HP_S", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "HP_S", fluidData))
#         lab.extend(["S"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "r"
#         for i in range(len(hatch)):
#             rects2.extend(ax2.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax2,rects2[int(2*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#         curList = []
#         if el==3: curList.append(getSingleData(fld, "INCOMP" , "HP_DST", fluidData))
#         curList.append(getSingleData(fld, "HEOS"   , "HP_DST", fluidData))
#         curList.append(getSingleData(fld, "REFPROP", "HP_DST", fluidData))
#         lab.extend(["all"])
#         offset += step
#         curN   = len(curList)
#         curInd = [ offset + i * width for i in range(curN)]
#         #curCol = "r"
#         for i in range(len(hatch)):
#             rects2.extend(ax2.bar(curInd[i], np.mean(curList[i]), width, color=curCol, hatch=hatch[i]))#, yerr=np.std(curList[i]), ecolor='k'))
#         autolabel(ax2,rects2[int(3*el):],np.mean(curList,axis=1),np.std(curList,axis=1))
#
#
#         ids = np.arange(len(lab))
#         # add some text for labels, title and axes ticks
#         #if backend=="INCOMP": ax1.set_ylabel(r'Time per $f(p,T)$ call (\si{\us})')
#         if el==3: ax1.set_ylabel(r'Time per \texttt{Tmax} call (\SI{0.01}{\us}) and'+"\n"+r'per $f(p,T)$ and $f(\rho,T)$ call (\si{\us})')
#         else: ax1.set_ylabel(r'Time per \texttt{Tmax} call (\SI{0.01}{\us})'+"\n"+r'and per $f(\rho,T)$ call (\si{\us})')
#         ax2.set_ylabel(r'Time per $f(h,p)$ call (\si{\us})')
#
#         ax1.set_xticks(ids)
#         ax1.set_xticklabels([r"\texttt{"+i+r"}" for i in lab], rotation=0)
#         ax1.set_xlim([ids.min()-2.5*width,ids.max()+2.5*width])
#
#         ax1.spines['top'].set_visible(False)
#         ax2.spines['top'].set_visible(False)
#         ax1.xaxis.set_ticks_position('bottom')
#         ax2.xaxis.set_ticks_position('bottom')
#
#         handles = []
#         if el==3:
#             for h in (rects1[0], rects1[4], rects2[2]):
#                 handles.append(copy.copy(h))
#                 handles[-1].set_facecolor('white')
#                 handles.append(copy.copy(h))
#                 handles[-1].set_hatch('')
#             labels = (r'$p,T$-fit', r'\texttt{Tmax}', r'CoolProp', r'explicit, $f(p|\rho,T)$', r'REFPROP', r'implicit, $f(h,p)$')
#         else:
#             for h in (rects1[0], rects1[2], rects2[1]):
#                 handles.append(copy.copy(h))
#                 handles[-1].set_facecolor('white')
#                 handles.append(copy.copy(h))
#                 handles[-1].set_hatch('')
#             labels = (r'', r'\texttt{Tmax}', r'CoolProp', r'explicit, $f(\rho,T)$', r'REFPROP', r'implicit, $f(h,p)$')
#             handles[0] = mpatches.Patch(visible=False)
#
#         if bp is not None:
#             bp.drawLegend(ax=ax1,
#               loc='upper center',
#               bbox_to_anchor=(0.5, 1.4),
#               ncol=3,
#               handles=handles,
#               labels=labels)
#         else:
#             ax1.legend(handles,labels,
#               loc='upper center',
#               bbox_to_anchor=(0.5, 1.4),
#               ncol=3)
#
#         fg.savefig(path.join(getFigureFolder(),"TimeComp-"+fld+".pdf"))
#         if bp is not None: bp.savepgf(path.join(getFigureFolder(),"TimeComp-"+fld+".pgf"),fg,repList)
#         plt.close()
#
#     except Exception as e:
#         print(e)
#         pass

#############################################################
# The second figure compares the backend for the full calculation
#############################################################
backendExp = []
backendImp = []
fluidLabel = []
hatchLst = ["///", "\\\\\\"]
backendsLst = ["HEOS", "REFPROP"]
for fluidstr in fluids[:-1]:
    _, fld, _ = splitFluid(fluidstr)
    #if fld=="CO2": fluidLabel.append("\ce{CO2}")
    # else:
    fluidLabel.append(fld)
    for backend in backendsLst:
        if backend not in fluidData[fld]:
            backendExp.append([0])
            backendImp.append([0])
            continue
        # Backend exists in fluid data
        try:
            inp = getInpList(backend)
            outExp = getTimeKey(inp[0], getOutList(inp[0])[-1])
            outImp = getTimeKey(inp[1], getOutList(inp[1])[-1])
            backendExp.append(getSingleData(fld, backend, outExp, fluidData))
            backendImp.append(getSingleData(fld, backend, outImp, fluidData))
        except Exception as e:
            backendExp.append([0])
            backendImp.append([0])
            print(e)
            pass

# Data is prepared, we can plot now.
if bp is not None:
    bp.figure = None
    fg1 = bp.getFigure()
    bp2 = BasePlotter()
    fg2 = bp2.getFigure()
    ccycle = bp.getColorCycle(length=3)
else:
    fg1 = plt.figure()
    fg2 = plt.figure()
    ccycle = cycle(["b", "g", "r"])

fg1.set_size_inches((fg1.get_size_inches()[0] * 1, fg1.get_size_inches()[1] * 0.75))
fg2.set_size_inches((fg2.get_size_inches()[0] * 1, fg2.get_size_inches()[1] * 0.75))
ccycle.next()  # No incomp
#
#ax1 = fg.add_subplot(111)
#ax2 = ax1.twinx()
ax1 = fg1.add_subplot(111)
ax2 = fg2.add_subplot(111)

#entries = int(len(backendExp)/len(fluidLabel))
# one standard deviation above and below the mean of the data
rects1 = []
rects2 = []
col1 = ccycle.next()
col2 = ccycle.next()

numFluids = len(fluidLabel)
numBackends = len(backendsLst)

step = 1
width = step / (numBackends + 1)  # the width of the bars
offset = -0.5 * numBackends * width

for f in range(numFluids):
    ids = np.empty((numBackends,))
    for b in range(numBackends):
        i = f * numBackends + b
        ids[b] = offset + f * step + b * width
        rects1.extend(ax1.bar(ids[b], np.mean(backendExp[i]), width, color=col1, hatch=hatchLst[b], rasterized=False))  # , yerr=np.std(curList[i]), ecolor='k'))
        rects2.extend(ax2.bar(ids[b], np.mean(backendImp[i]), width, color=col2, hatch=hatchLst[b], rasterized=False))  # , yerr=np.std(curList[i]), ecolor='k'))

y_max = np.max(np.concatenate((np.ravel(ax1.get_ylim()), np.ravel(ax2.get_ylim()) / 10.0)))

ax1.set_ylim([0, y_max])
ax2.set_ylim([0, y_max * 10.0])

for ax in [ax1, ax2]:
    ax.set_xticks(range(numFluids))
    ax.set_xticklabels(fluidLabel, rotation=25)
    ax.set_xlim([0.0 - 0.5 * step, numFluids - 1 + 0.5 * step])
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')

    if ax == ax1: rects = rects1
    elif ax == ax2: rects = rects2
    handles = (rects[0], rects[1])
    labels = (r'CoolProp', r'REFPROP')
    anchor = (0.5, 1.2)

    if bp is not None:
        bp.drawLegend(ax=ax,
          loc='upper center',
          bbox_to_anchor=anchor,
          ncol=numBackends,
          handles=handles,
          labels=labels)
    else:
        ax.legend(handles, labels,
          loc='upper center',
          bbox_to_anchor=anchor,
          ncol=numBackends)

ax1.set_ylabel(r'Time per $f(\rho,T)$ call (us)')
ax2.set_ylabel(r'Time per $f(h,p)$ call (us)')
fg1.savefig(path.join(getFigureFolder(), "TimeComp-backends-exp.pdf"))
fg2.savefig(path.join(getFigureFolder(), "TimeComp-backends-imp.pdf"))
if bp is not None:
    ax1.set_ylabel(r'Time per $f(\rho,T)$ call (\si{\us})')
    ax2.set_ylabel(r'Time per $f(h,p)$ call (\si{\us})')
    mpl.rcParams['text.usetex'] = True
    fg1.savefig(path.join(getFigureFolder(), "TimeComp-backends-exp-tex.pdf"))
    fg2.savefig(path.join(getFigureFolder(), "TimeComp-backends-imp-tex.pdf"))
    mpl.rcParams['text.usetex'] = False
    bp.savepgf(path.join(getFigureFolder(), "TimeComp-backends-exp.pgf"), fg1, repList)
    bp.savepgf(path.join(getFigureFolder(), "TimeComp-backends-imp.pgf"), fg2, repList)
plt.close('all')

#############################################################
# The third figure is a heat map of the execution times in
# log p h diagram
#############################################################
for fluidstr in fluids:
    try:
        _, fld, _ = splitFluid(fluidstr)
        for backend in fluidData[fld]:
            propsfluid = "::".join([backend, fluidstr])

            if backend != "INCOMP":
                TP = {}
                points = max(int(xypoints / 2), 250)
                T_range_TP = np.linspace(PropsSI('Ttriple', "T", 0, "D", 0, propsfluid) + 1, PropsSI('Tcrit', "T", 0, "D", 0, propsfluid) - 0.1, points)
                T_TP = np.append(T_range_TP, T_range_TP[::-1])
                Q_TP = np.zeros_like(T_TP)
                Q_TP[points:] = 1
                points *= 2

                out = ["D", "H", "P", "S"]
                res = PropsSI(out, "T", T_TP, "Q", Q_TP, propsfluid)
                D_TP = res[:, 0]
                H_TP = res[:, 1]
                P_TP = res[:, 2]
                S_TP = res[:, 3]

                mask = np.isfinite(D_TP)
                if np.sum(mask) < points:
                    warn("There were not enough valid entries in your result vector. Reducing the number of points from {0:d} to {1:d}.".format(points, np.sum(mask)))
                    points = np.sum(mask)

                TP["T"] = T_TP[mask]
                TP["D"] = D_TP[mask]
                TP["H"] = H_TP[mask]
                TP["P"] = P_TP[mask]
                TP["S"] = S_TP[mask]
                # saveNpzData(TP)
            else:
                TP = None

            state = getStateObj(propsfluid)
            if backend == "HEOS" and state.has_melting_line():
                p_melt = np.logspace(np.log10(state.melting_line(CoolProp.constants.iP_min, CoolProp.constants.iT, 0)), np.log10(state.melting_line(CoolProp.constants.iP_max, CoolProp.constants.iT, 0)), xypoints)
                #p_melt = p_range
                ML = dict(T=[], D=[], H=[], S=[], P=p_melt)
                for p in p_melt:
                    try:
                        ML["T"].append(state.melting_line(CoolProp.constants.iT, CoolProp.constants.iP, p))
                    except Exception as ve:
                        ML["T"].append(np.inf)
                res = PropsSI(["D", "H", "P", "S", "T"], "T", ML["T"], "P", ML["P"], propsfluid)
                ML["D"] = res[:, 0]
                ML["H"] = res[:, 1]
                ML["P"] = res[:, 2]
                ML["S"] = res[:, 3]
                ML["T"] = res[:, 4]
                mask = np.isfinite(ML["T"])
                ML["P"] = ML["P"][mask]
                ML["T"] = ML["T"][mask]
                ML["D"] = ML["D"][mask]
                ML["H"] = ML["H"][mask]
                ML["S"] = ML["S"][mask]
            else:
                ML = None
            #ML = {}

            IP = {}
            p_range, T_range = getPTRanges(propsfluid)
            critProps = getCriticalProps(propsfluid)
            try:
                IP["T"] = T_range
                IP["P"] = np.zeros_like(T_range) + critProps["P"]
                res = PropsSI(["D", "H"], "T", IP["T"], "P", IP["P"], propsfluid)
                IP["D"] = res[:, 0]
                IP["H"] = res[:, 1]
            except Exception as ve:
                IP = None

            IT = {}
            try:
                IT["P"] = p_range
                IT["T"] = np.zeros_like(p_range) + critProps["T"]
                res = PropsSI(["D", "H"], "T", IT["T"], "P", IT["P"], propsfluid)
                IT["D"] = res[:, 0]
                IT["H"] = res[:, 1]
            except Exception as ve:
                IT = None

            ID = {}
            try:
                ID["T"] = T_range
                ID["D"] = np.zeros_like(p_range) + critProps["D"]
                res = PropsSI(["P", "H"], "T", ID["T"], "D", ID["D"], propsfluid)
                ID["P"] = res[:, 0]
                ID["H"] = res[:, 1]
            except Exception as ve:
                ID = None

            IH = {}
            try:
                IH["P"] = p_range
                IH["H"] = np.zeros_like(p_range) + critProps["H"]
                res = PropsSI(["D", "T"], "P", IH["P"], "H", IH["H"], propsfluid)
                IH["D"] = res[:, 0]
                IH["T"] = res[:, 1]
            except Exception as ve:
                IH = None

            IS = {}
            try:
                IS["P"] = p_range
                IS["S"] = np.zeros_like(p_range) + critProps["S"]
                res = PropsSI(["D", "H", "T"], "P", IS["P"], "S", IS["S"], propsfluid)
                IS["D"] = res[:, 0]
                IS["H"] = res[:, 1]
                IS["T"] = res[:, 2]
            except Exception as ve:
                IS = None

            for I in [IP, IT, ID, IH, IS]:
                if I is not None:
                    mask = np.isfinite(I["D"]) & np.isfinite(I["H"])
                    if np.sum(mask) < 20: I = None
                    else:
                        for k in I:
                            I[k] = I[k][mask]

            for inp in getInpList(backend):
                if bp is not None:
                    bp.figure = None
                    fg = bp.getFigure()
                else:
                    fg = plt.figure()

                kind = getTimeKey(inp, getOutList(inp)[-1])
                t_data = getSingleData(fld, backend, kind, fluidData)
                x_data = getSingleData(fld, backend, "H", fluidData)
                y_data = getSingleData(fld, backend, "P", fluidData)

                gs = gridspec.GridSpec(1, 2, wspace=None, hspace=None, width_ratios=[10, 1])
                ax1 = fg.add_subplot(gs[0, 0], axisbg='Tan')
                ax1.set_yscale('log')
                #ax2 = ax1.twinx()

                minHP = np.min(t_data)
                maxHP = np.max(t_data)
                minIT = np.percentile(t_data, 10)
                maxIT = np.percentile(t_data, 90)
                difIT = np.log10(maxIT / minIT) * 0.25

                print(kind, ": {0:7.2f} to {1:7.2f}".format(minHP, maxHP))
                if kind == "DT":
                    if bp is not None:
                        cx1 = bp.getColourMap(reverse=True)
                    else:
                        cx1 = mpl.cm.get_cmap('cubehelix_r')
                    minHP = minIT
                    maxHP = np.power(10, np.log10(maxIT) + difIT)
                    #minHP = np.power(10,np.log10(np.percentile(t_data,10)*1e6))
                    #maxHP = np.power(10,np.log10(np.percentile(t_data,90)*1e6)*1.10)
                    #maxHP = np.power(10,1.10*np.log10(maxHP))
                    #minHP = np.percentile(t_data,10)*1e6
                    #maxHP = np.percentile(t_data,99)*1e6
                    #print(kind,": {0:7.2f} to {1:7.2f}".format(minHP,maxHP))
                    #minHP = 100
                    #maxHP = 20000
                else:
                    if bp is not None:
                        cx1 = bp.getColourMap()
                    else:
                        cx1 = mpl.cm.get_cmap('cubehelix')
                    minHP = np.power(10, np.log10(minIT) - difIT)
                    maxHP = maxIT
                    #minHP = np.power(10,np.log10(np.percentile(t_data,10)*1e6)*0.90)
                    #maxHP = np.power(10,np.log10(np.percentile(t_data,90)*1e6))
                    # minHP = np.percentile(t_data,01)*1e6
                    #maxHP = np.percentile(t_data,90)*1e6
                    #print(kind,": {0:7.2f} to {1:7.2f}".format(minHP,maxHP))
                    #minHP = 100
                    #maxHP = 20000

                #cx1_r = reverse_colourmap(cx1)
                cNorm = mpl.colors.LogNorm(vmin=minHP, vmax=maxHP)
                #cNorm  = mpl.colors.LogNorm(vmin=ceil(minHP/1e1)*1e1, vmax=floor(maxHP/1e2)*1e2)
                #cNorm  = mpl.colors.Normalize(vmin=round(minHP,-2), vmax=round(maxHP,-2))
                colourSettings = dict(c=t_data, edgecolors='none', cmap=cx1, norm=cNorm)
                pointSettings = dict(s=6)
                scatterSettings = dict(rasterized=True, alpha=0.5)
                #scatterSettings = dict(rasterized=False, alpha=0.5)
                scatterSettings.update(colourSettings)
                scatterSettings.update(pointSettings)
                SC = ax1.scatter(x_data / 1e6, y_data / 1e5, **scatterSettings)

                for I in [TP, ML]:
                    if I is not None:
                        ax1.plot(I["H"] / 1e6, I["P"] / 1e5, lw=1.5, c='k')

                for I in [IP, IT, ID, IS, IH]:
                    if I is not None:
                        ax1.plot(I["H"] / 1e6, I["P"] / 1e5, lw=1.0, c='k', alpha=1)

                # ax1.set_xlim([0e+0,6e1])
                # ax1.set_ylim([5e-1,2e4])
                ax1.set_xlim([np.percentile(x_data / 1e6, 0.1), np.percentile(x_data / 1e6, 99.9)])
                ax1.set_ylim([np.percentile(y_data / 1e5, 0.1), np.percentile(y_data / 1e5, 99.9)])

                formatter = ticker.LogFormatter(base=10.0, labelOnlyBase=False)
                #formatter = ticker.ScalarFormatter()
                #ticks     = roundList(np.logspace(np.log10(ax1.get_ylim()[0]), np.log10(ax1.get_ylim()[1]), 5))
                #locator   = ticker.FixedLocator(ticks)
                # ax1.yaxis.set_major_locator(locator)
                ax1.yaxis.set_major_formatter(formatter)

                cax = fg.add_subplot(gs[0, 1])
                formatter = ticker.ScalarFormatter()
                CB = fg.colorbar(SC, cax=cax, format=formatter)
                CB.set_alpha(1)
                CB.locator = ticker.MaxNLocator(nbins=7)
                #ticks = roundList(np.logspace(np.log10(minHP), np.log10(maxHP), 5))
                #CB.locator = ticker.FixedLocator(ticks)
                CB.update_ticks()
                CB.draw_all()

                # fg.suptitle("f("+inp+")-"+backend.lower()+"-"+fld.lower())
                CB.set_label(backend.upper() + "::" + fld + ', execution time per f(' + inp[0] + "," + inp[1] + ') call (us)')
                ax1.set_xlabel(r'Specific enthalpy (MJ/kg)')
                ax1.set_ylabel(r'Pressure (bar)')
                fg.tight_layout()
                fg.savefig(path.join(getFigureFolder(), "TimeComp-" + inp + "-" + backend.lower() + "-" + fld.lower() + ".pdf"))

                #CB.set_label(r'Execution time per call (\si{\us})')
                #ax1.set_xlabel(r'Specific enthalpy (\si{\mega\J\per\kg})')
                #ax1.set_ylabel(r'Pressure (\si{\bar})')
                # fg.tight_layout()
                # bp.savepgf(path.join(getFigureFolder(),"TimeComp-"+inp+"-"+backend.lower()+"-"+fld.lower()+".pgf"),fg,repList)
                plt.close()

    except Exception as e:
        print(e)
        pass
    plt.close('all')
