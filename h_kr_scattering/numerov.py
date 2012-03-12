#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2012 Tiziano Müller <tm@dev-zero.ch>
#
#
#
#

def numerov(dx, kn, psin, kn_p, psin_p, kn_n):
    g  = dx**2/12.
    c0 = 1. + g*kn_p
    c1 = 2. - 10.*g*kn
    c2 = 1. + g*kn_n
    return (c1*psin - c0*psin_p) / c2;
