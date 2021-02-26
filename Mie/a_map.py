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

    f = np.linspace(c0/lambda_start, c0/lambda_stop, X_TICKS, dtype=np.float64)
    k = 2 * np.pi * np.sqrt(beta0 * rho0) * f

    a_ = np.linspace(a_start, a_stop, Y_TICKS)

    sca = []

    for a in a_:
        sca.append(sigmaSc_(k, a, rho0, beta0, rho1, beta1))
    plt.subplot(aspect='equal')
    plt.pcolor(sca, norm=NORM, cmap=COLOR_MAP)

    lambda_ticks= np.linspace(lambda_start, lambda_stop, round(X_TICKS/GAP_IN_TICKS))
    a_ticks = np.linspace(a_start, a_stop, round(Y_TICKS / GAP_IN_TICKS))
    plt.xticks([GAP_IN_TICKS*i+0.5 for i in range(round(X_TICKS/GAP_IN_TICKS))], np.round(lambda_ticks, 3), rotation='vertical')
    plt.yticks([GAP_IN_TICKS*i+0.5 for i in range(round(Y_TICKS/GAP_IN_TICKS))], np.round(a_ticks, 3))
    plt.xlabel('freq')
    plt.ylabel('a')
    plt.title('c0 = %.2f, rho0 = %.2f, c1 = %.2f, rho1 = %.2f'%(c0, rho0, c1, rho1))


if __name__ == '__main__':
    a = 0.01
    rho0 = 1.39
    c0 = 331.45
    beta0 = 1 / (c0**2 * rho0)
    rho1 = 97
    c1 = 123
    beta1 = 1 / (c1**2 * rho1)

    X_TICKS = 400   # кол-во пикселей по Х
    Y_TICKS = 400   # кол-во пикселей по У
    GAP_IN_TICKS = 20   # расстояние между подписями на графике

    a_start = 1
    a_stop = 2

    lambda_start = 0.4
    lambda_stop = 2

    fig, ax = plt.subplots()
    fig.set_size_inches(16,9)

    COLOR_MAP = cm.jet
    NORM = colors.Normalize()

    update_graph()
    plt.show()








