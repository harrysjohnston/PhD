from cosmosis.postprocessing.cosmology_theory_plots import Plot, plot_list
from cosmosis.postprocessing import lazy_pylab as pylab

class PmIplot(Plot):
    filename='matter_intrinsic_power'
    def plot(self):
        super(PmIPlot,self).plot()
        k = self.load_file("matter_intrinsic_power", "k_h")
        p_k = self.load_file("matter_intrinsic_power", "p_k")
        pylab.loglog(k_h, p_k)
        pylab.xlabel("k / (Mpc/h)")
        pylab.ylabel("P(k) / (h^-1 Mpc)^3")