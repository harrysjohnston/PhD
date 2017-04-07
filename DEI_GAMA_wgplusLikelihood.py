from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class DEI_GAMA_wgplusLikelihood(GaussianLikelihood):

    def build_data(self):
        #use self.options to find the data_file and load ell, tt from it
        data_file = self.options.get_string("data_file")
        try:
            rp,wgp = np.loadtxt(data_file).T[:2]
        except ValueError:
            rp,wgp = np.loadtxt(data_file,skiprows=1).T[:2]
            print('ValueError thrown for %s; skipping row 1..'%data_file)
        return rp,wgp

    def build_covariance(self):
        cov_file = self.options.get_string("cov_file")
        try:
            covmat = np.loadtxt(cov_file)
        except ValueError:
            covmat = np.loadtxt(cov_file,skiprows=1)
            print('ValueError thrown for %s; skipping row 1..'%cov_file)
        return covmat

setup,execute,cleanup = DEI_GAMA_wgplusLikelihood.build_module()