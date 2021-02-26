import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def sphericalh1(n, z, p=0):
    if p:
        return np.sqrt(np.pi / (2 * z)) * ((n / z) * sp.hankel1(n + 0.5, z) - sp.hankel1(n + 1.5, z))
    else:
        return np.sqrt(np.pi / (2 * z)) * sp.hankel1(n + 0.5, z)


def an_(n, ka, rho1, beta1):
    gamma = np.sqrt(beta1 / rho1)
    k1a = ka * np.sqrt(beta1 * rho1)
    jn1 = sp.spherical_jn(n, k1a)
    jn = sp.spherical_jn(n, ka)
    jn1p = sp.spherical_jn(n, k1a, 1)
    jnp = sp.spherical_jn(n, ka, 1)
    hn = sphericalh1(n, ka)
    hnp = sphericalh1(n, ka, 1)
    up = (gamma * jn1p * jn - jn1 * jnp)
    down = (jn1 * hnp - gamma * jn1p * hn)
    if np.isnan(up / down):
        return 0
    else:
        return up / down


def sigmaSc_(k, a, rho0, beta0, rho1, alpha, nmin=0, nmax=50):
    ka = a * k
    sum = np.zeros(ka.size, dtype=np.float64)

    i = 0
    for ka_ in ka:
        sum_ = 0
        beta1 = (1j * alpha / (ka_ / a) - 1)**2 / (c1 ** 2 * rho1)
        for n in range(nmin, nmax + 1):
            an = an_(n, ka_, rho1 / rho0, beta1 / beta0)
            sum_ = sum_ + (2 * n + 1) * abs(an ** 2)
        sum[i] = (4 * np.pi / ((ka_ / a) ** 2) * sum_)
        i += 1
    return sum


def update_graph():
    ax.clear()
    ax.grid()

    f = np.linspace(f_start, f_stop, X_TICKS, dtype=np.float64)
    k = 2 * np.pi * np.sqrt(beta0 * rho0) * f
    alpha_ = np.linspace(alpha_start, alpha_stop, Y_TICKS)

    sca = []

    for alpha in alpha_:
        sca.append(sigmaSc_(k, a, rho0, beta0, rho1, alpha))
    plt.subplot(aspect='equal')
    plt.pcolor(sca, norm=NORM, cmap=COLOR_MAP)

    f_ticks = np.linspace(f_start, f_stop, round(X_TICKS / GAP_IN_TICKS))
    c_ticks = np.linspace(alpha_start, alpha_stop, round(Y_TICKS / GAP_IN_TICKS))
    plt.xticks([GAP_IN_TICKS * i + 0.5 for i in range(round(X_TICKS / GAP_IN_TICKS))], np.round(f_ticks, 3),
               rotation='vertical')
    plt.yticks([GAP_IN_TICKS * i + 0.5 for i in range(round(Y_TICKS / GAP_IN_TICKS))], np.round(c_ticks, 3))
    plt.title('c0 = %.2f, rho0 = %.2f, c1 = %.2f, rho1 = %.2f' % (c0, rho0, c1, rho1))
    plt.xlabel('f, Гц')
    plt.ylabel('alpha')


if __name__ == '__main__':
    a = 0.01
    rho0 = 1.39
    c0 = 331
    beta0 = 1 / (c0 ** 2 * rho0)
    c1 = 123
    rho1 = 97

    X_TICKS = 200  # кол-во пикселей по Х
    Y_TICKS = 50  # кол-во пикселей по У
    GAP_IN_TICKS = 10  # расстояние между подписями на графике

    alpha_start = 0
    alpha_stop = 5

    f_start = 1000
    f_stop = 20000

    fig, ax = plt.subplots()
    fig.set_size_inches(16, 9)

    COLOR_MAP = cm.inferno
    NORM = colors.Normalize(vmax=0.01)

    update_graph()
    plt.show()
