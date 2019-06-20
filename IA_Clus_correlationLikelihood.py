from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names, option_section
import numpy as np

class DEI_GAMA_wgplusLikelihood(GaussianLikelihood): # called 'wgplus' but works on wgg as well

	def __init__(self, options):
		self.options = options
		self.data_x, self.data_y = self.build_data()
		self.likelihood_only = options.get_bool('likelihood_only', False)
		if options.has_value("Njkregions"):
			self.Njkregions = options["Njkregions"]
			self.testsample = 0
		elif options.has_value("Nsims"):
			self.Njkregions = options["Nsims"]
			self.testsample = 0
		else:
			print('No number of realisations in options for Hartlap factor!!',
				'\nbetter be doing TEST SAMPLER!!')
			self.testsample = 1
################################################################
# PULLED FROM GAUSSIAN_LIKELIHOOD.__INIT__
		if self.constant_covariance:
			self.cov = self.build_covariance()
			self.inv_cov = self.build_inverse_covariance()
			if not self.likelihood_only:
				self.chol = np.linalg.cholesky(self.cov)
			include_norm = self.options.get_bool("include_norm", False)
			if include_norm:
				self.log_det_constant = GaussianLikelihood.extract_covariance_log_determinant(self,None)
				print("Including -0.5*|C| normalization in {} likelihood where |C| = {}".format(self.like_name, self.log_det_constant))
			else:
				self.log_det_constant = 0.0
		self.kind = self.options.get_string("kind", "cubic")
		# Allow over-riding where the inputs come from in 
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
################################################################

	def build_data(self):
		# use self.options to find the data_file and load r_p and wg+
		data_file = self.options.get_string("data_file")
		self.NLA = self.options.get_bool("NLA")
		self.drop_large_bins = self.options.get_int("drop_large_bins", default=0)
		try:
			rp,wgp = np.loadtxt(data_file).T[:2]
		except ValueError:
			rp,wgp = np.loadtxt(data_file,skiprows=1).T[:2]
			#print('ValueError thrown for %s; skipping row 1..'%data_file)

		# limit x-range of data for fitting
		if self.NLA:
			self.rplim = self.options.get_double('rplim', default=6.)
			rplim = self.rplim

			print('fitting NLA model; discard r_p < %s Mpc/h...!'%rplim)
			self.nla_cov_cut = sum(rp>rplim)
			rp, wgp = rp[rp>rplim], wgp[rp>rplim]
		if self.drop_large_bins!=0:
			print('DISCARDING %s largest r_p bins'%self.drop_large_bins)
			rp, wgp = rp[:-self.drop_large_bins], wgp[:-self.drop_large_bins]

		# define relevant x-range for theory curve
		self.data_x_range = (rp.min()*0.8, rp.max()*1.2)
		print('plot points: r_p = %s'%(rp))
		self.Nrpbins = len(rp)
		return rp, wgp

	def build_covariance(self):
		cov_file = self.options.get_string("cov_file")
		self.diag_only = self.options.get_bool("diag_only", default=False)
		try:
			covmat = np.loadtxt(cov_file)
		except ValueError:
			covmat = np.loadtxt(cov_file,skiprows=1)
			#print('ValueError thrown for %s; skipping row 1..'%cov_file)

		# clip covariance into range for fitting
		if self.NLA:
			covmat = covmat[-self.nla_cov_cut:,-self.nla_cov_cut:]
		if self.drop_large_bins!=0:
			covmat = covmat[:-self.drop_large_bins,:-self.drop_large_bins]

		if self.diag_only:
			print('==================== ===================== ===================== DIAGONAL COVARIANCE ONLY!!')
			covmat = np.diag(covmat) * np.identity(len(covmat))

		return covmat

	def build_inverse_covariance(self):
		# apply Hartlap factor
		invcov = np.linalg.inv(self.cov)
		if not self.testsample:
			Hf = (self.Njkregions - self.Nrpbins - 2.) / (self.Njkregions - 1.)
			invcov *= Hf
		return invcov

	def extract_theory_points(self, block):
		"Extract relevant theory from block and get theory at data x values"
		theory_x = block[self.x_section, self.x_name]
		theory_y = block[self.y_section, self.y_name]
		range_cut = (theory_x > self.data_x_range[0]) & (theory_x < self.data_x_range[1])
		theory_x, theory_y = theory_x[range_cut], theory_y[range_cut]

		return self.generate_theory_points(theory_x, theory_y)

setup,execute,cleanup = DEI_GAMA_wgplusLikelihood.build_module()
