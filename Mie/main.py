import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox


class SliderClass:
    def __init__(self, index, name, min, max, start):
        self.slider_ax = fig.add_axes([0.05, 0.84 - 0.04*index, 0.24, 0.02])
        self.slider = Slider(self.slider_ax, name, np.real(min), np.real(max), valinit=np.real(start))
        self.slider.on_changed(on_change_graph)


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

    multipoles = [[] for i in range(nmin, nmax+1)]

    i = 0
    for ka_ in ka:
        sum_ = 0
        for n in range(nmin, nmax+1):
            an = an_(n, ka_, rho1/rho0, beta1/beta0)
            multipoles[n].append(4 * np.pi / ((ka_/a) ** 2) * (2*n+1)*abs(an**2))
            sum_ = sum_ + (2*n+1)*abs(an**2)
        sum[i] = (4 * np.pi / ((ka_/a) ** 2) * sum_)
        i+=1
    return sum, multipoles


def update_graph():
    f = np.linspace(f_start, f_stop, 200, dtype=np.float64)
    k = 2 * np.pi * np.sqrt(beta0 * rho0) * f
    sca, multipoles = sigmaSc_(k, a, rho0, beta0, rho1, beta1)
    ax.clear()
    ax.grid()
    ax.set_xlabel('f, Гц')
    ax.set_ylabel('σ')
    ax.plot(f, sca)
    legend = ['σ']
    for i in range(MULTIPOLES):
        ax.plot(f, multipoles[i])
        legend.append('n=%i'%i)
    ax.legend(legend)


def on_change_c_particle(value):
    global beta1, beta1_slider
    beta1 = 1/((float(value))**2 * rho1)
    beta1_slider.slider.set_val(beta1)


def on_change_c_host(value):
    global beta0, beta0_slider
    beta0 = 1/((float(value))**2 * rho0)
    beta0_slider.slider.set_val(beta0)


def on_change_rho_particle(value):
    global rho1_slider
    rho1_slider.slider.set_val(float(value))


def on_change_rho_host(value):
    global rho0_slider
    rho0_slider.slider.set_val(float(value))


def on_change_graph(value):
    global a, rho0, beta0, f_stop, rho1, beta1, f_start
    a = a_slider.slider.val
    f_start = f_start_slider.slider.val
    f_stop = f_stop_slider.slider.val
    beta0 = beta0_slider.slider.val
    rho0 = rho0_slider.slider.val
    beta1 = beta1_slider.slider.val
    rho1 = rho1_slider.slider.val
    # с_particle.set_val(round(1/np.sqrt(rho1*beta1), 5))
    # c_host.set_val(round(1 / np.sqrt(rho0 * beta0), 5))
    rho_particle.set_val(round(rho1, 5))
    rho_host.set_val(round(rho0, 5))
    refractive_index.set_val(round(np.sqrt((rho1*beta1)/(rho0*beta0)),7))
    n_particle.set_val(round(np.sqrt(rho1*beta1),5))
    n_host.set_val(round(np.sqrt(rho0*beta0),5))
    update_graph()


if __name__ == '__main__':
    a = 0.01
    rho0 = 1.39
    c0 = 331.45
    beta0 = 1 / (c0 ** 2 * rho0)
    rho1 = 97
    c1 = 123
    beta1 = 1 / (c1**2 * rho1)
    f_start = 1
    f_stop = 15000

    MULTIPOLES = 5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.4)
    update_graph()

    c_particle = TextBox(fig.add_axes([0.15, 0.5, 0.14, 0.03]), 'с1', initial=round(1 / np.sqrt(rho1 * beta1), 6))
    c_host = TextBox(fig.add_axes([0.15, 0.46, 0.14, 0.03]), 'с0', initial=round(1/np.sqrt(rho0*beta0),6))
    c_particle.on_submit(on_change_c_particle)
    c_host.on_submit(on_change_c_host)

    rho_particle = TextBox(fig.add_axes([0.15, 0.42, 0.14, 0.03]), 'rho1', initial=round(rho1, 6))
    rho_host = TextBox(fig.add_axes([0.15, 0.38, 0.14, 0.03]), 'rho0', initial=round(rho0, 6))
    rho_particle.on_submit(on_change_rho_particle)
    rho_host.on_submit(on_change_rho_host)

    refractive_index = TextBox(fig.add_axes([0.15, 0.32, 0.14, 0.03]), 'refractive_index', initial=round(np.sqrt((rho1*beta1)/(rho0*beta0)),7))
    n_particle = TextBox(fig.add_axes([0.15, 0.28, 0.14, 0.03]), 'n_particle', initial=round(np.sqrt(rho1*beta1),5))
    n_host = TextBox(fig.add_axes([0.15, 0.24, 0.14, 0.03]), 'n_host', initial=round(np.sqrt(rho0*beta0),5))

    f_start_slider = SliderClass(0, 'f_start', min=0, max=2 * f_start, start=f_start)
    f_stop_slider = SliderClass(1, 'f_max', min=0, max=2 * f_stop, start=f_stop)
    rho0_slider = SliderClass(2, 'rho0', min=0, max=10 * rho0, start=rho0)
    beta0_slider = SliderClass(3, 'beta0', min=0, max=10 * beta0, start=beta0)
    rho1_slider = SliderClass(4, 'rho1', min=0, max=10 * rho1, start=rho1)
    beta1_slider = SliderClass(5, 'beta1', min=0, max=10 * beta1, start=beta1)
    a_slider = SliderClass(6, 'a', min=0, max=5 * a, start=a)

    plt.show()







