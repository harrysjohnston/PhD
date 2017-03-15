from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class wgg_likelihood(GaussianLikelihood):

    x_section = self.options.get_string("theory_section")
    x_name    = "theta"
    y_section = x_section
    y_name    = "wgg_r"
    like_name = "wgg"

    def build_data(self):
        #use self.options to find the data_file and load ell, tt from it
        data_file = self.options.get_string("data_file")
        print("loading wgg measurement, NOT skipping rows...!")
        rp,wgg = np.loadtxt(data_file).T[:2]
        return rp,wgg

    def build_covariance(self):
        cov_file = self.options.get_string("cov_file")
        print("loading covariance, skipping row=1 (column heads)...!")
        covmat = np.loadtxt(cov_file,skiprows=1)
        return covmat

setup,execute,cleanup = wgg_likelihood.build_module()