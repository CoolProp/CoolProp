# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

from CoolProp.Plots import PropertyPlot  # TODO: Change to absolute import


def main():
    fluid_ref = 'n-Pentane'
    for plot_type in ['Ts']:  # ['pt', 'ph', 'ps', 'ts', 'pt', 'prho', 'trho']:
        plt = PropertyPlot(fluid_ref, plot_type)
        plt.set_axis_limits([-0.5 * 1e3, 1.5 * 1e3, 300, 530])
        plt.draw_isolines('Q', [0.1, 0.9])
        plt.draw_isolines('P', [100 * 1e3, 2000 * 1e3])
        plt.draw_isolines('D', [2, 600])
        plt.show()


if __name__ == "__main__":
    main()
