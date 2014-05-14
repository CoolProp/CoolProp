# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 18:39:22 2013

@author: logan
"""

from Plots import PropsPlot #TODO: Change to absolute import


def main():
    fluid_ref = 'n-Pentane'
    for plot_type in ['Ts']: #['pt', 'ph', 'ps', 'ts', 'pt', 'prho', 'trho']:
        plt = PropsPlot(fluid_ref, plot_type)
        plt.set_axis_limits([-0.5, 1.5, 300, 530])
        plt.draw_isolines('Q', [0.1, 0.9])
        plt.draw_isolines('P', [100, 2000])
        plt.draw_isolines('D', [2, 600])
        plt.show()

if __name__ == "__main__":
    main()
