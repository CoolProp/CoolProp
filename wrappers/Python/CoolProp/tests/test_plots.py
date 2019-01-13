import numpy as np
import matplotlib.pyplot as plt


def test_back_compatibility():
    fluid_ref = 'R290'

    def Ts_plot_tests():
        from CoolProp.Plots import Ts
        Ts(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ts(fluid_ref, show=False, axis=ax)
        plt.close()
        Ts(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def Ph_plot_tests():
        from CoolProp.Plots import Ph
        Ph(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ph(fluid_ref, show=False, axis=ax)
        plt.close()
        Ph(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def PT_plot_tests():
        from CoolProp.Plots import PT
        PT(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        PT(fluid_ref, show=False, axis=ax)
        plt.close()
        PT(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def Ps_plot_tests():
        from CoolProp.Plots import Ps
        Ps(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ps(fluid_ref, show=False, axis=ax)
        plt.close()
        Ps(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def Prho_plot_tests():
        from CoolProp.Plots import Prho
        Prho(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Prho(fluid_ref, show=False, axis=ax)
        plt.close()
        Prho(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def Trho_plot_tests():
        from CoolProp.Plots import Trho
        Trho(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Trho(fluid_ref, show=False, axis=ax)
        plt.close()
        Trho(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def hs_plot_tests():
        from CoolProp.Plots import hs
        hs(fluid_ref, show=False)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        hs(fluid_ref, show=False, axis=ax)
        plt.close()
        hs(fluid_ref, show=False, Tmin=200, Tmax=300)
        plt.close()

    def Isolines_plot_tests():
        from matplotlib import pyplot
        from CoolProp.Plots import Ts, drawIsoLines
        ax = Ts(fluid_ref)
        #ax.set_xlim([-0.5, 1.5])
        #ax.set_ylim([300, 530])
        quality = drawIsoLines(fluid_ref, 'Ts', 'Q', [0.3, 0.5, 0.7, 0.8], axis=ax)
        isobars = drawIsoLines(fluid_ref, 'Ts', 'P', [100, 2000], num=5, axis=ax)
        isochores = drawIsoLines(fluid_ref, 'Ts', 'D', [2, 600], num=7, axis=ax)
        pyplot.close()

    Ts_plot_tests()
    Ph_plot_tests()
    Ps_plot_tests()
    PT_plot_tests()
    Prho_plot_tests()
    Trho_plot_tests()
    hs_plot_tests()
    Isolines_plot_tests()


def test_new_code():
    fluid_ref = 'Water'

    def Ts_plot_tests():
        from CoolProp.Plots import PropsPlot
        PP = PropsPlot(fluid_ref, 'Ts')
        plt.close()

    def Ph_plot_tests():
        from CoolProp.Plots import PropsPlot
        PP = PropsPlot(fluid_ref, 'Ph')
        plt.close()

    def Isolines_plot_tests():
        from CoolProp.Plots import PropsPlot
        PP = PropsPlot(fluid_ref, 'Ts')
        #plt.set_axis_limits([-0.5, 1.5, 300, 530])
        PP.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
        PP.draw_isolines('P', [100, 2000], num=5)
        PP.draw_isolines('D', [2, 600], num=7)
        plt.close()

    def Graph_annotations():
        from CoolProp.Plots import PropsPlot, IsoLines
        PP = PropsPlot(fluid_ref, 'Ts')
        PP.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
        PP.draw_isolines('P', [100, 2000], num=5)
        PP.draw_isolines('D', [2, 600], num=7)
        plt.title('New Title')
        PP.xlabel('New x label')
        PP.ylabel('New y label')
        PP = IsoLines(fluid_ref, 'Ts', 'P')
        PP.draw_isolines([100, 2000], num=5)
        plt.close()

    def Mixture():
        from CoolProp.Plots import PropsPlot
        PP = PropsPlot('REFPROP-MIX:R32[0.47319469]&R125[0.2051091]&R134a[0.32169621]', 'TD')
        PP._plot_default_annotations()
        plt.close()

    Ts_plot_tests()
    Ph_plot_tests()
    Isolines_plot_tests()
    Graph_annotations()
    Mixture()


if __name__ == '__main__':
    import nose
    nose.runmodule()
