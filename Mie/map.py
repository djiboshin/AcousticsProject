import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib import cm, colors


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


def sigmaSc_(k, a, rho0, beta0, rho1, beta1, nmin=0, nmax=50):
    ka = a * k
    sum = np.zeros(ka.size, dtype=np.float64)

    i = 0
    for ka_ in ka:
        sum_ = 0
        for n in range(nmin, nmax+1):
            an = an_(n, ka_, rho1/rho0, beta1/beta0)
            sum_ = sum_ + (2*n+1)*abs(an**2)
        sum[i] = (4 * np.pi / ((ka_/a) ** 2) * sum_)
        i+=1
    return sum


def update_graph():
    ax.clear()
    ax.grid()
    ax.set_xlabel('f, Гц')
    ax.set_ylabel('c')

    f = np.linspace(f_start, f_stop, X_TICKS, dtype=np.float64)
    k = 2 * np.pi * np.sqrt(beta0 * rho0) * f
    c_ = np.linspace(c1_start, c1_stop, Y_TICKS)

    sca = []

    for c in c_:
        beta1 = 1 / (c**2 * rho1)
        sca.append(sigmaSc_(k, a, rho0, beta0, rho1, beta1))
    plt.subplot(aspect='equal')
    plt.pcolor(sca, norm=NORM, cmap=COLOR_MAP)

    f_ticks = np.linspace(f_start, f_stop, round(X_TICKS/GAP_IN_TICKS))
    c_ticks = np.linspace(c1_start, c1_stop, round(Y_TICKS / GAP_IN_TICKS))
    plt.xticks([GAP_IN_TICKS*i+0.5 for i in range(round(X_TICKS/GAP_IN_TICKS))], np.round(f_ticks, 3), rotation='vertical')
    plt.yticks([GAP_IN_TICKS*i+0.5 for i in range(round(Y_TICKS/GAP_IN_TICKS))], np.round(c_ticks, 3))
    plt.xlabel('freq')
    plt.ylabel('c1')
    plt.title('c0 = %.2f, rho0 = %.2f, rho1 = %.2f'%(c0, rho0, rho1))


if __name__ == '__main__':
    a = 1
    rho0 = 1.39
    c0 = 331
    beta0 = 1 / (c0**2 * rho0)
    rho1 = 2

    X_TICKS = 100   # кол-во пикселей по Х
    Y_TICKS = 100   # кол-во пикселей по У
    GAP_IN_TICKS = 10   # расстояние между подписями на графике

    c1_start = 100
    c1_stop = 500

    f_start = 30
    f_stop = 2000

    fig, ax = plt.subplots()
    fig.set_size_inches(16,9)

    COLOR_MAP = cm.inferno
    NORM = colors.Normalize(vmax=40)

    update_graph()
    plt.show()








