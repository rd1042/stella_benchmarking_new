""" """

import scipy.integrate as integrate
import scipy.special as special
import numpy as np

result = integrate.quad(lambda x: special.jv(1,x)/x, 0, np.inf)
print("result = ", result)
