import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox


def sphericalh1(n, z, p=0):
    if p:
        return np.sqrt(np.pi / (2*z)) * ((n / z) * sp.hankel1(n+0.5, z) - sp.hankel1(n+1.5, z))
    else:
        return np.sqrt(np.pi / (2*z)) * sp.hankel1(n+0.5, z)


def an_(n, ka, rho1, beta1):
    gamma = np.sqrt(beta1/rho1)
    k1a = ka * np.sqrt(beta1*rho1)
    jn1 = sp.spherical_jn(n, k1a)
    jn = sp.spherical_jn(n, ka)
    jn1p = sp.spherical_jn(n, k1a, 1)
    jnp = sp.spherical_jn(n, ka, 1)
    hn = sphericalh1(n, ka)
    hnp = sphericalh1(n, ka, 1)
    up = (gamma * jn1p * jn - jn1 * jnp)
    down = (jn1 * hnp - gamma * jn1p * hn)
    if np.isnan(up/down):
        return 0
    else:
        return up/down


def sigmaSc_(k, a, rho0, beta0, rho1, alpha1, nmin=0, nmax=50):
    ka = a * k
    sum = np.zeros(ka.size, dtype=np.float64)

    multipoles = [[] for i in range(nmin, nmax+1)]

    i = 0
    for ka_ in ka:
        beta1 = (1j * alpha1 / (ka_ / a) - 1) ** 2 / (c1 ** 2 * rho1)
        sum_ = 0
        for n in range(nmin, nmax+1):
            an = an_(n, ka_, rho1/rho0, beta1/beta0)
            multipoles[n].append(4 * np.pi / ((ka_/a) ** 2) * (2*n+1)*abs(an**2))
            sum_ = sum_ + (2*n+1)*abs(an**2)
        sum[i] = (4 * np.pi / ((ka_/a) ** 2) * sum_)
        i += 1
    return sum, multipoles


def update_graph():
    f = np.linspace(f_start, f_stop, 1000, dtype=np.float64)
    k = 2 * np.pi * np.sqrt(beta0 * rho0) * f
    sca, multipoles = sigmaSc_(k, a, rho0, beta0, rho1, alpha1)
    plt.grid()
    plt.xlabel('f, Гц')
    plt.ylabel('σ')
    plt.plot(f, sca)
    legend = ['σ']
    for i in range(MULTIPOLES):
        plt.plot(f, multipoles[i])
        legend.append('n=%i'%i)
    plt.legend(legend)
    plt.title('$c_0$ = %.2f, $\\rho_0$ = %.2f, $c_1$ = %.2f, $\\rho_1$ = %.2f, $\\alpha_1$ = %.2f, $e^{-\\alpha_1*a}$ = %.2f' % (c0, rho0, c1, rho1, alpha1, np.exp(-alpha1 * a)))


if __name__ == '__main__':
    a = 0.02
    rho0 = 1.39
    c0 = 331.45
    beta0 = 1 / (c0 ** 2 * rho0)
    rho1 = 97
    c1 = 123
    alpha1 = 10
    f_start = 1000
    f_stop = 20000

    MULTIPOLES = 5

    fig = plt.figure()
    update_graph()
    plt.show()







