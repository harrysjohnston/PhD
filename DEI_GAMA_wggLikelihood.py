from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np

class DEI_GAMA_wggLikelihood(GaussianLikelihood):

    def __init__(self, options):
	self.options=options
        self.data_x, self.data_y = self.build_data()
	if options.has_value("Nsims"):
            self.Nsims = options["Nsims"]
	    self.testsample = 0
        else:
            print('No number of mock sims in options for Hartlap factor!!'
                        '\nbetter be doing TEST SAMPLER!!')
            self.testsample = 1
        if self.constant_covariance:
            self.cov = self.build_covariance()
            self.inv_cov = self.build_inverse_covariance()
        self.kind = self.options.get_string("kind", "cubic")

        #Allow over-riding where the inputs come from in 
        #the options section
        if options.has_value("x_section"):
            self.x_section = options['x_section']
        if options.has_value("y_section"):
            self.y_section = options['y_section']
        if options.has_value("x_name"):
            self.x_name = options['x_name']
        if options.has_value("y_name"):
            self.y_name = options['y_name']
        if options.has_value("like_name"):
            self.like_name = options['like_name']

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
        print('plot points: r_p = %s'%(rp))
        self.Nrpbins = len(rp)
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
