#!/usr/bin/env python

# Copyright 2017 Perttu Luukko

# This file is part of libeemd.
# 
# libeemd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# libeemd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with libeemd.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt

def plot_complex_signal(s, axis):
    """Common helper function for plotting complex signals"""
    axis.plot(np.real(s), 'b-')
    axis.plot(np.imag(s), 'k--')

# Load output data
data = np.genfromtxt("bemd_example.out", dtype=complex)
x = np.real(data[0,:])
y = np.imag(data[0,:])

plt.figure()
plt.plot(x, y)

# Plot imfs
f, axes = plt.subplots(5, sharex=True)

plot_complex_signal(data[0,:], axes[0])
for i in range(4):
    plot_complex_signal(data[i+1,:], axes[i+1])

plt.show()
