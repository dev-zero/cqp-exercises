#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
#
#
#
#

from numpy import zeros, exp, arange, sqrt, array, arctan2, sin, arctan, pi
from scipy.special import sph_jn, sph_yn
from bisect import bisect_left
from numerov import numerov

def LennardJonesPotential(e, s, r):
    return e*((s/r)**12 - 2.*(s/r)**6)

class HKrScattering:
    def __init__(self, a, b, N=512):
        self._N = N
        self._u = zeros(N)

        self._dr = (b-a)/float(N)
        self._r = a + arange(N)*self._dr

        self._V = zeros(N)

    def solve(self, E=0.5):
        e = 5.9
        s = 1. # 3.57 A, but we normalized to s
        r0 = .5
        mh = 6.12

        # find the position of the initial value for the numerov algorithm
        r0_pos = bisect_left(self._r, r0)
        # .... and initialize the solution up to that point using the solution for small r
        self._u[:r0_pos+1] = map(lambda r: exp(-sqrt(6.12*e/25.)/r**5), self._r[:r0_pos+1])

        def solve4l(l):
            def CombinedPotential(r):
                return LennardJonesPotential(e, s, r) + l*(l+1)/(mh*r**2)

            # initialize the (combined) potential
            self._V = array(map(lambda r: CombinedPotential(r), self._r))

            # initialize the k's, based on the potential
            self._k = mh*(E-self._V)

            for i in xrange(r0_pos, len(self._r)-1):
                self._u[i+1] = numerov(self._dr, self._k[i], self._u[i], self._k[i-1], self._u[i-1], self._k[i+1])


        l_max = 10
        s_tot = 0.
        for l in range(0, l_max+1):
            solve4l(l)

            r1 = self._r[-4]
            u1 = self._u[-4]
            r2 = self._r[-1]
            u2 = self._u[-1]

            K = r1*u2/(r2*u1)

            k = sqrt(E*mh)

            delta = arctan2( K*sph_jn(l, k*r1)[0][l] - sph_jn(l, k*r2)[0][l] , K*sph_yn(l, k*r1)[0][l] - sph_yn(l, k*r2)[0][l] )
            s_tot += (2*l + 1)*sin(delta)**2 * 4*pi/k**2

        print s_tot
        return s_tot

    def u(self):
        return self._u

    def r(self):
        return self._r

    def V(self):
        return self._V

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="H-Kr Scattering")
    args = parser.parse_args()

    from matplotlib.pyplot import show, plot, legend, figure

    s = HKrScattering(0.3, 6., 1024)

    s.solve()

    E_range = arange(.1, 3.6, .05)
    s_tot = map(lambda E: s.solve(E), E_range) 
    plot(E_range, abs(array(s_tot)))
    show()
