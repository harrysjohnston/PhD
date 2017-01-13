from cosmosis.postprocessing.cosmology_theory_plots import Plot, plot_list
from cosmosis.postprocessing import lazy_pylab as pylab

class PmIplot(MatterPowerPlot):
    filename='matter_intrinsic_power'
    def plot(self):
        super(PmIplot,self).plot()
        done_any=False
        if os.path.exists("{0}/matter_intrinsic_power".format(self.dirname)):
            self.plot_section("matter_intrinsic_power", "matter-intrinsic")
            done_any=True
        if not done_any:
            raise IOError("Not making plot: %s (no data in this sample)"% self.__class__.__name__[:-4])
        pylab.xlabel("k / (Mpc/h)")
        pylab.ylabel("P(k) / (h^-1 Mpc)^3")
        pylab.grid()
        pylab.legend()