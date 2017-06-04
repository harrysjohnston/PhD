from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class DEI_GAMA_wggLikelihood(GaussianLikelihood):

    def __init__(self, options):
	if options.has_value("Nsims"):
            self.Nsims = options["Nsims"]
        else:
            print('No number of mock sims in options for Hartlap factor!!\n
                        better be doing TEST SAMPLER!!')
            self.testsample = 1

    def build_data(self):
        #use self.options to find the data_file and load ell, tt from it
        data_file = self.options.get_string("data_file")
	self.NLA = self.options.get_bool("NLA")
        try:
            rp,wgp = np.loadtxt(data_file).T[:2]
        except ValueError:
            rp,wgp = np.loadtxt(data_file,skiprows=1).T[:2]
            print('ValueError thrown for %s; skipping row 1..'%data_file)
	if self.NLA:
		rplim = 2 # Mpc/h
                print('fitting NLA model; discard r_p < %s Mpc/h...!'%rplim)
		self.nla_cov_cut = sum(rp>rplim)
                rp, wgp = rp[rp>rplim], wgp[rp>rplim]
        self.data_x_range = (rp.min()*0.8, rp.max()*1.2)
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
        return covmat

    def build_inverse_covariance(self):
        # apply Hartlap factor
        invcov = np.linalg.inv(self.cov)
	if not self.testsample:
        	Hf = (self.Nsims-1)/(self.Nsims-self.Nrpbins-2)
        	invcov *= Hf
        return invcov

    def extract_theory_points(self, block):
        "Extract relevant theory from block and get theory at data x values"
        theory_x = block[self.x_section, self.x_name]
        theory_y = block[self.y_section, self.y_name]
        range_cut = (theory_x > self.data_x_range[0]) & (theory_x < self.data_x_range[1])
        theory_x, theory_y = theory_x[range_cut], theory_y[range_cut]

        return self.generate_theory_points(theory_x, theory_y)

setup,execute,cleanup = DEI_GAMA_wggLikelihood.build_module()
