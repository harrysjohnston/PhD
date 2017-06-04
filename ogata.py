# coding: utf-8
tanh, sinh, pi = np.tanh, np.sinh, np.pi
pi
import scipy.special as ss
Jnu = ss.jv
xi = ss.jn_zeros
w_nul = 2 / (pi**2 * xi(2,5) * Jnu(3, pi*xi(2,5)))
w_nul
k = np.loadtxt('./cosmosis/modules/euclid_ias/HM_HOD_TEST/matter_intrinsic_power_1_1/ell.txt')
k.shape
pk = np.loadtxt('./cosmosis/modules/euclid_ias/HM_HOD_TEST/matter_intrinsic_power_1_1/p_k.txt')
pk.shape
k = np.loadtxt('./cosmosis/modules/euclid_ias/HM_HOD_TEST/matter_intrinsic_power_1_1/k_h.txt')
get_ipython().magic(u'rerun 9')
pk = pk[0]
psi = lambda t: t*tanh((pi/2)*sinh(t))
def fivept_stencil(func,x,h):
    # returns f'(x), via 5pt stencil, for grid-spacing h
    return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)
np.diff(xi)
np.diff(xi(2,5))
w_nul
spline = scipy.interpolate.UnivariateSpline
import scipy
spline = scipy.interpolate.UnivariateSpline
import scipy.interpolate
spline = scipy.interpolate.UnivariateSpline
kpk_spline = spline(k,k*pk)
get_ipython().magic(u'matplotlib')
plt.plot(k, k*pk)
