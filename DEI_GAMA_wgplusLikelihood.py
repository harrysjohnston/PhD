from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class DEI_GAMA_wgplusLikelihood(GaussianLikelihood):

    def build_data(self):
        #use self.options to find the data_file and load ell, tt from it
        data_file = self.options.get_string("data_file")
        self.NLA = self.options.get_bool("NLA")
	self.drop_large_bins = self.options.get_int("drop_large_bins", default=0)
        try:
            rp,wgp = np.loadtxt(data_file).T[:2]
        except ValueError:
            rp,wgp = np.loadtxt(data_file,skiprows=1).T[:2]
            print('ValueError thrown for %s; skipping row 1..'%data_file)
	if self.NLA:
		print('fitting NLA model; discard r_p < 6 Mpc/h...!')
		self.nla_cov_cut = sum(rp>6)
		rp, wgp = rp[rp>6], wgp[rp>6]
	print('DISCARDING %s largest r_p bins'%self.drop_large_bins)
	if self.drop_large_bins!=0:
		rp, wgp = rp[:-self.drop_large_bins], wgp[:-self.drop_large_bins]
        self.data_x_range = (rp.min()*0.8, rp.max()*1.2)
	print('RP SHAPE: %s'%rp.shape, self.drop_large_bins)
        return rp, wgp

    def build_covariance(self):
        cov_file = self.options.get_string("cov_file")
        try:
            covmat = np.loadtxt(cov_file)
        except ValueError:
            covmat = np.loadtxt(cov_file,skiprows=1)
            print('ValueError thrown for %s; skipping row 1..'%cov_file)
	if self.NLA:
		covmat = covmat[-self.nla_cov_cut:,-self.nla_cov_cut:]
	if self.drop_large_bins!=0:
                covmat = covmat[:-self.drop_large_bins,:-self.drop_large_bins]
        return covmat

    def extract_theory_points(self, block):
        "Extract relevant theory from block and get theory at data x values"
        theory_x = block[self.x_section, self.x_name]
        theory_y = block[self.y_section, self.y_name]
	range_cut = (theory_x > self.data_x_range[0]) & (theory_x < self.data_x_range[1])
	theory_x, theory_y = theory_x[range_cut], theory_y[range_cut]
	
        return self.generate_theory_points(theory_x, theory_y)

setup,execute,cleanup = DEI_GAMA_wgplusLikelihood.build_module()
