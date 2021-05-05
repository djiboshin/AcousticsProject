import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import csv
from matplotlib import gridspec
from matplotlib import artist
import math
import scipy.special as spf
from matplotlib.widgets import Slider, TextBox


class SliderClass:
    def __init__(self, index, name, min, max, start):
        self.slider_ax = fig.add_axes([0.05, 0.84 - 0.04*index, 0.24, 0.02])
        self.slider = Slider(self.slider_ax, name, np.real(min), np.real(max), valinit=np.real(start))
        self.slider.on_changed(on_change_graph)


def f(x):
    return 0.5*(x-1j*(g1+g2)+np.sqrt((x-1j*(g2-g1))**2+4*(v+1j*v_im)**2))


def g(x):
    return 0.5*(x-1j*(g1+g2)-np.sqrt((x-1j*(g2-g1))**2+4*(v+1j*v_im)**2))


def on_change_graph(value):
    global g1, g2, v, E0, v_im
    g1 = g1_slider.slider.val
    g2 = g2_slider.slider.val
    v = v_slider.slider.val
    E0 = E0_slider.slider.val
    v_im = v_im_slider.slider.val
    update_graph()


def update_graph():
    ax1.clear()
    ax1.grid()
    ax1.set_xlabel('$\Delta E$')
    ax1.set_ylabel('Re(freq)')
    ax1.plot(x, np.real(f(x)))
    ax1.plot(x, np.real(g(x)))

    ax2.clear()
    ax2.grid()
    ax2.set_xlabel('$\Delta E$')
    ax2.set_ylabel('|Im(freq)|')
    ax2.plot(x, np.abs(np.imag(f(x))))
    ax2.plot(x, np.abs(np.imag(g(x))))

    ax3.clear()
    ax3.grid()
    ax3.set_xlabel('$\Delta E$')
    ax3.set_ylabel('Q')
    ax3.plot(x, np.divide(np.abs(np.real(f(x)) + E0),np.abs(np.imag(f(x)))))
    ax3.plot(x, np.divide(np.abs(np.real(g(x)) + E0),np.abs(np.imag(g(x)))))


x = np.linspace(-20,20,3000)

fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
fig.subplots_adjust(left=0.3)
g1 = 0.1
g2 = 2
v = 1
v_im = 2
E0 = 500

g1_slider = SliderClass(0, 'g1', min=0, max=2 * g1, start=g1)
g2_slider = SliderClass(1, 'g2', min=0, max=2 * g2, start=g2)
v_slider = SliderClass(2, 'Re(v)', min=0, max=10 * v, start=v)
v_im_slider = SliderClass(3, 'Im(v)', min=-2 * v_im, max=2 * v_im, start=v_im)
E0_slider = SliderClass(4, 'E0', min=0, max=10 * E0, start=E0)

update_graph()
plt.show()
