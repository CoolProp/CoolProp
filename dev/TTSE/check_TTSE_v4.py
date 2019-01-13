import CoolProp.CoolProp as CP
import matplotlib
matplotlib.rc('font', family='serif', serif='Times New Roman')
#from matplotlib2tikz import save as tikz_save

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker
from matplotlib.patches import Ellipse
from matplotlib.transforms import ScaledTranslation
import numpy as np
import random
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from matplotlib.gridspec import GridSpec

# Create the colourmap
#import numpy as np
#import matplotlib.pyplot as plt
import matplotlib._cm, matplotlib.cm
specs = matplotlib._cm.cubehelix(gamma=1.4, s=0.4, r=-0.8, h=2.0)
specs_r = matplotlib.cm._reverse_cmap_spec(specs)
matplotlib.cm.register_cmap(name="jorrithelix", data=specs)
matplotlib.cm.register_cmap(name="jorrithelix" + "_r", data=specs_r)


def makeGrid(x, y, z, resX=200, resY=200):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z


def getErrors(p, h, out='D', Ref=''):
    "Get the relative errors from table-based interpolation"
    errorTTSE = 1e3
    errorBICUBIC = 1e3
    try:
        # Using the EOS
        CP.disable_TTSE_LUT(Ref)
        EOS = CP.PropsSI(out, 'P', p, 'H', h, Ref)
        # Using the TTSE method
        CP.enable_TTSE_LUT(Ref)
        CP.set_TTSE_mode(Ref, "TTSE")
        TTSE = CP.PropsSI(out, 'P', p, 'H', h, Ref)
        # Using the Bicubic method
        CP.enable_TTSE_LUT(Ref)
        CP.set_TTSE_mode(Ref, "BICUBIC")
        BICUBIC = CP.PropsSI(out, 'P', p, 'H', h, Ref)
        errorTTSE = abs(TTSE / EOS - 1.0) * 100.0
        errorBICUBIC = abs(BICUBIC / EOS - 1.0) * 100.0
    except ValueError as VE:
        print(VE)
        pass

    return errorTTSE, errorBICUBIC


# ['YlOrRd', 'PuBuGn', 'hot', 'cubehelix', 'gnuplot', 'gnuplot2']:
for colourmap in ['jorrithelix']:

    for out in ['D']:
        # landscape figure
        #fig = plt.figure(figsize=(10,5))
        #ax1 = fig.add_axes((0.08,0.1,0.32,0.83))
        #ax2 = fig.add_axes((0.50,0.1,0.32,0.83))
        #cbar_ax = fig.add_axes([0.80, 0.075, 0.05, 0.875])

        # portrait figure
        #fig = plt.figure(figsize=(5,8))
        #ax1     = plt.subplot2grid((2,8), (0,0), colspan=7)
        #ax2     = plt.subplot2grid((2,8), (1,0), colspan=7)
        #cbar_ax = plt.subplot2grid((2,8), (0,7), colspan=1, rowspan=2)

        #fig = plt.figure(figsize=(8,4))
        #ax1     = plt.subplot2grid((1,7), (0,0), colspan=3)
        #ax2     = plt.subplot2grid((1,7), (0,3), colspan=3)
        #cbar_ax = plt.subplot2grid((1,7), (0,6), colspan=1, rowspan=1)
        # plt.tight_layout()
        fig = plt.figure(figsize=(8, 4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        #cbar_ax = plt.subplot2grid((1,7), (0,6), colspan=1, rowspan=1)
        # plt.tight_layout()

        #Ref = 'R245fa'
        #Ref = 'Isopentane'
        Ref = 'Air'

        T = np.linspace(CP.PropsSI(Ref, 'Tmin') + 0.1, CP.PropsSI(Ref, 'Tcrit') - 0.01, 300)
        pV = CP.PropsSI('P', 'T', T, 'Q', 1, Ref)
        hL = CP.PropsSI('H', 'T', T, 'Q', 0, Ref)
        hV = CP.PropsSI('H', 'T', T, 'Q', 1, Ref)
        hTP = np.append(hL, [hV[::-1]])
        pTP = np.append(pV, [pV[::-1]])

        HHH1, PPP1, EEE1 = [], [], []
        HHH2, PPP2, EEE2 = [], [], []

        cNorm = colors.LogNorm(vmin=1e-10, vmax=1e-1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap(colourmap))

        # Setting the limits for enthalpy and pressure
        p_min = CP.PropsSI(Ref, 'ptriple')
        p_max = 60e5
        h_min = CP.PropsSI('H', 'T', CP.PropsSI(Ref, 'Ttriple') + 0.5, 'Q', 0, Ref)
        h_max = CP.PropsSI('H', 'T', 500 + 273.15, 'P', p_max, Ref)

        # Creating some isotherms for better illustration of the cycle
        isoT = np.array([0, 100, 200, 300, 400]) + 273.15
        isoP = np.logspace(np.log10(p_min), np.log10(p_max), base=10)
        ones = np.ones(isoP.shape)
        isoH = [CP.PropsSI('H', 'T', T * ones, 'P', isoP, Ref) for T in isoT]

        print("Lower left and upper right coordinates: ({0},{1}), ({2},{3})".format(h_min, p_min, h_max, p_max))

        CP.set_TTSESinglePhase_LUT_range(Ref, h_min, h_max * 1.05, p_min, p_max * 1.05)

        for a_useless_counter in range(40000):

            h = random.uniform(h_min, h_max)
            p = 10**random.uniform(np.log10(p_min), np.log10(p_max))

            try:
                # Using the EOS
                CP.disable_TTSE_LUT(Ref)
                rhoEOS = CP.PropsSI('D', 'P', p, 'H', h, Ref)
                TEOS = CP.PropsSI('T', 'P', p, 'H', h, Ref)
                if out == 'C': cpEOS = CP.PropsSI('C', 'P', p, 'H', h, Ref)

                # Using the TTSE method
                CP.enable_TTSE_LUT(Ref)
                CP.set_TTSE_mode(Ref, "TTSE")
                rhoTTSE = CP.PropsSI('D', 'P', p, 'H', h, Ref)
                TTTSE = CP.PropsSI('T', 'P', p, 'H', h, Ref)
                if out == 'C': cpTTSE = CP.PropsSI('C', 'P', p, 'H', h, Ref)

                # Using the Bicubic method
                CP.enable_TTSE_LUT(Ref)
                CP.set_TTSE_mode(Ref, "BICUBIC")
                rhoBICUBIC = CP.PropsSI('D', 'P', p, 'H', h, Ref)
                TBICUBIC = CP.PropsSI('T', 'P', p, 'H', h, Ref)
                if out == 'C': cpBICUBIC = CP.PropsSI('C', 'P', p, 'H', h, Ref)

                if out == 'D':
                    errorTTSE = abs(rhoTTSE / rhoEOS - 1) * 100
                    errorBICUBIC = abs(rhoBICUBIC / rhoEOS - 1) * 100
                elif out == 'T':
                    errorTTSE = abs(TTTSE / TEOS - 1) * 100
                    errorBICUBIC = abs(TBICUBIC / TEOS - 1) * 100
                elif out == 'C':
                    errorTTSE = abs(cpTTSE / cpEOS - 1) * 100
                    errorBICUBIC = abs(cpBICUBIC / cpEOS - 1) * 100

                HHH1.append(h)
                PPP1.append(p)
                EEE1.append(errorTTSE)

                HHH2.append(h)
                PPP2.append(p)
                EEE2.append(errorBICUBIC)

            except ValueError as VE:
                # print VE
                pass

        HHH1 = np.array(HHH1)
        PPP1 = np.array(PPP1)
        SC1 = ax1.scatter(HHH1 / 1e3, PPP1 / 1e5, s=8, c=EEE1, edgecolors='none', cmap=plt.get_cmap(colourmap), norm=cNorm, rasterized=True)

        #X, Y, Z = makeGrid(HHH1, np.log10(PPP1), EEE1)
        # SC1 = matplotlib.pyplot.contourf(X, Y, Z,
        #                    alpha=0.75,
        #                    norm=cNorm,
        #                    cmap=matplotlib.pyplot.get_cmap(colourmap))#,
        #                    #rasterized=True)
        HHH2 = np.array(HHH2)
        PPP2 = np.array(PPP2)
        SC2 = ax2.scatter(HHH2 / 1e3, PPP2 / 1e5, s=8, c=EEE2, edgecolors='none', cmap=plt.get_cmap(colourmap), norm=cNorm, rasterized=True)

        if out == 'D':
            ax1.set_title('rel. density error, TTSE')
            ax2.set_title('rel. density error, bicubic')
        elif out == 'T':
            ax1.set_title('rel. temperature error, TTSE')
            ax2.set_title('rel. temperature error, bicubic')
        elif out == 'C':
            ax1.set_title('rel. heat capacity error, TTSE')
            ax2.set_title('rel. heat capacity error, bicubic')

        for ax in [ax1, ax2]:
            #h_min = np.ceil(h_min)
            delta = 0.1
            delta_min = 1.0 + delta
            delta_max = 1.0 - delta

            #ax.set_xlim(delta_min*h_min/1e3, delta_max*h_max/1e3)
            #ax.set_ylim(delta_min*p_min/1e5, delta_max*p_max/1e5)
            ax.set_xlim(-155, 800)
            ax.set_ylim(0.025, 58)

            ax.set_yscale('log')

            #ticks = np.array([0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50])
            ticks = np.array([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50])
            labels = [str(tick) for tick in ticks]
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

            #ticks = [150,250,350,450,550]
            #labels = [str(tick) for tick in ticks]
            # ax.set_xticks(ticks)
            # ax.set_xticklabels(labels)

            # ax.tick_params(axis='y',which='minor',left='off')

            #ax.set_xlabel('Enthalpy [kJ \cdot kg^{-1}]')
            ax.set_xlabel('Specific Enthalpy [kJ$\cdot$kg$\mathdefault{^{-1}\!}$]')
            ax.set_ylabel('Pressure [bar]')

            #ax.plot(hL/1e3,pV/1e5,'k',lw = 4)
            #ax.plot(hV/1e3,pV/1e5,'k',lw = 4)

            ax.plot(hTP / 1e3, pTP / 1e5, 'k', lw=3)

            for i, T in enumerate(isoT):
                ax.plot(isoH[i] / 1e3, isoP / 1e5, 'k', lw=1)

        #CB = fig.colorbar(SC1)
        #cbar_ax = fig.add_axes([0.80, 0.075, 0.05, 0.875])
        #CB = fig.colorbar(SC1, cax=cbar_ax)

        #CB = matplotlib.pyplot.colorbar(SC2)
        # CB.solids.set_rasterized(True)
        # ax2.yaxis.set_visible(False)

        # [x0,y0,width,height]
        #cbar_ax = fig.add_axes([0.95, 0.00, 0.05, 1.00])
        #CB = fig.colorbar(SC2, ax=[ax1,ax2], cax=cbar_ax)
        # CB.solids.set_rasterized(True)

        #from mpl_toolkits.axes_grid1 import make_axes_locatable
        #divider = make_axes_locatable(ax2)
        #cbar_ax = divider.append_axes("right", "5%", pad="0%")
        #CB = plt.colorbar(SC2, cax=cbar_ax)
        # CB.solids.set_rasterized(True)

        #CB = fig.colorbar(SC2)
        # CB.solids.set_rasterized(True)

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax2)
        ax_cb = divider.new_horizontal(size="5%", pad=0.05)
        #fig1 = ax.get_figure()
        fig.add_axes(ax_cb)
        CB = fig.colorbar(SC2, cax=ax_cb)

        #aspect = 5./2.
        # ax1.set_aspect(aspect)
        # ax2.set_aspect(aspect)

        CB.solids.set_rasterized(True)

        if out == 'D':
            CB.set_label(r'$\|\rho/\rho\mathdefault{_{EOS}-1\|\times 100}$ [%]')
        elif out == 'T':
            CB.set_label(r'$\|T/T\mathdefault{_{EOS}-1\|\times 100}$ [%]')
        elif out == 'C':
            CB.set_label(r'$\|c\mathdefault{_p}/c\mathdefault{_{p,EOS}-1\|\times 100}$ [%]')

        # The plot is finished, now we add an ellipse
        # circle=plt.Circle((5,5),.5,color='b',fill=False)
        # A scale-free ellipse.
        # xy - center of ellipse
        # width - total length (diameter) of horizontal axis
        # height - total length (diameter) of vertical axis
        #angle - rotation in degrees (anti-clockwise)
        p_op_min = 1e5
        p_op_max = 3e5
        h_op_min = CP.PropsSI('H', 'T', 400 + 273.15, 'P', p_op_max, Ref)
        h_op_max = CP.PropsSI('H', 'T', 25 + 273.15, 'P', p_op_max, Ref)

        p_op_cen = (p_op_min + p_op_max) / 2.0
        h_op_cen = (h_op_min + h_op_max) / 2.0

        p_op_hei = p_op_max - p_op_min
        h_op_wid = h_op_max - h_op_min

        # for ax in [ax1, ax2]:
        ##x,y = 10,0
        # use the axis scale tform to figure out how far to translate
        ##circ_offset = ScaledTranslation(x,y,ax.transScale)
        # construct the composite tform
        ##circ_tform = circ_offset + ax.transLimits + ax.transAxes
        # ellipse = Ellipse(xy=(h_op_cen,p_op_cen), width=h_op_wid, height=p_op_hei, angle=15, color='black')#, transform=circ_tform)
        # ax.add_artist(ellipse)

#        font_def = font_manager.FontProperties(family='Helvetica', style='normal',
#    size=sizeOfFont, weight='normal', stretch='normal')
#
#        for a in fig.axes:
#            for label in [a.get_xticklabels(), a.get_yticklabels()]:
#                label.set_fontproperties(ticks_font

        #plt.savefig(out+'_'+colourmap+'_TTSE_BICUBIC.png', dpi = 300, transparent = True)
        # plt.savefig(out+'_'+colourmap+'_TTSE_BICUBIC.eps')
        #    plt.savefig(out+'_'+colourmap+'_TTSE_BICUBIC.pdf')
        plt.tight_layout()
        plt.savefig('check_TTSE_' + colourmap + '.pdf')
        #tikz_save(  'check_TTSE.tikz')
        #plt.savefig(out+'_'+colourmap+'_TTSE_BICUBIC.jpg', dpi = 1200)
        plt.close()
