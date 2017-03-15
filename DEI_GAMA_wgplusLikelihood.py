from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class DEI_GAMA_wgplusLikelihood(GaussianLikelihood):

    def build_data(self):
        #use self.options to find the data_file and load ell, tt from it
        data_file = self.options.get_string("data_file")
        print("loading wgp measurement, skipping row=1 (column heads)...!")
        rp,wgp = np.loadtxt(data_file,skiprows=1).T[:2]
        return rp,wgp

    def build_covariance(self):
        cov_file = self.options.get_string("cov_file")
        print("loading covariance, skipping row=1 (column heads)...!")
        covmat = np.loadtxt(cov_file,skiprows=1)
        return covmat

setup,execute,cleanup = DEI_GAMA_wgplusLikelihood.build_module()