#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
#
#
#
#

from numpy import zeros, exp, arange, sqrt
from bisect import bisect_left

def numerov(dx, kn, psin, kn_p, psin_p, kn_n):
    g  = dx**2/12.
    c0 = 1. + g*kn_p
    c1 = 2. - 10.*g*kn
    c2 = 1. + g*kn_n
    return (c1*psin - c0*psin_p) / c2;

class OneDQuantumScatteringNumerov:
    def __init__(self, N=256, L=10):
        self._N = N
        self._psi = zeros(N, complex)

        self._dx = 2.*L/N
        self._x = (arange(N) - N/2)*self._dx

        self._V = zeros(N)

        self._m = .5
        self._h_bar = 1.

    def solve(self, E=0.5, a=1.5):

        l = bisect_left(self._x, 0.)
        r = bisect_left(self._x, a)
        print "a = ", a, " shifted to ", self._x[r]

        self._V[l:r] = 1.

        self._k = 2.*self._m*(E-self._V)/self._h_bar**2

        self._psi[r:] = map(lambda x, k: exp(-1j*x*sqrt(k)), self._x[r:], self._k[r:])

        for i in range(r, l-1, -1):
            self._psi[i-1] = numerov(self._dx, self._k[i], self._psi[i], self._k[i+1], self._psi[i+1], self._k[i-1])
        
        # calculate the coefficients A, B for the pre-barrier solution
        p1 = self._psi[l]
        p2 = self._psi[l-1]
        x = self._x[l]
        y = self._x[l-1]
        c1 = exp(1j*sqrt(abs(self._k[l]))*x)
        c2 = exp(1j*sqrt(abs(self._k[l-1]))*y)
        A = (c1*c2*(c2*p1-c1*p2))/(-c1*c1 + c2*c2);
        B = (c1*p1 - c2*p2)/(c1*c1 - c2*c2);

        # the transmission and reflection probabilities
        T = 1./abs(A)**2
        R = T*abs(B)**2

        # plot the solution
        self._psi[:l] = map(lambda x, k: A*exp(-1j*x*sqrt(k)) + B*exp(1j*x*sqrt(k)), self._x[:l], self._k[:l])

        return (T, R)

    def psi(self):
        return self._psi

    def x(self):
        return self._x

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="solve the 1D quantum scattering problem using the Numerov algorithm and calculate the transmission/reflection probabilities")
    args = parser.parse_args()

    from matplotlib.pyplot import show, plot, legend, figure

    s = OneDQuantumScatteringNumerov(512)

    fig = figure()

    ax = fig.add_subplot(211)
    for E in arange(.1, 1., .1):
        s.solve(E)
        plot(s.x(), s.psi().real, label="E = %f" % E)

    legend(loc=1, borderaxespad=0.)

    ax = fig.add_subplot(212)

    a = arange(.0, 6., .05)
    T = map(lambda a: s.solve(a=a)[0], a)
    print T
    plot(a, T)

    show()
