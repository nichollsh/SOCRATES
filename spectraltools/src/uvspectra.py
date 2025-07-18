from scipy.interpolate import interp1d
import numpy as np
import os

import src.phys as phys
import src.utils as utils

def interpolation(wavelength, flux, step):

    ius = interp1d(wavelength, flux, kind='linear')
    grid_new = np.arange(min(wavelength), max(wavelength), step)
    flux_new = ius(grid_new)

    return grid_new, flux_new

def readUV(formula:str):
    path = os.path.join(utils.dirs["moleculesUV"], formula+'.txt')

    data = np.loadtxt(path)
    wvl = data[:, 0]
    xsec = data[:, 1]

    return wvl, xsec

def plot(wvl, xsec, show=True, saveout="plot_xsec_wavenumber", label="xsec"):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.plot(wvl, xsec, label=label)
    ax.set_yscale('log')
    ax.set_xlabel("Wavenumber [cm$^{-1}$]")

    if len(saveout) > 0:
        save_path = os.path.join(utils.dirs["moleculesUV"], saveout)
        print("Saving plot to '%s'"%save_path)
        utils.rmsafe(save_path)
        fig.savefig(save_path, bbox_inches="tight")
        show=True
        plt.close()

    if show:
        plt.show()
        return

#wvl, xsec = readUV(formula)
#plot(wvl, xsec, show=True, saveout="plot_xsec_wavenumber", label="xsec")

from netCDF4 import Dataset

nc = Dataset(os.path.join(utils.dirs["tools"], "output/Original_H2O.nc"), "r")
print(nc.variables.keys())

with nc as ds:
    nu = ds.variables["nu"][:]
    xsec = ds.variables["kabs"][:]

    idx_band1 = (nu >= 93170) & (nu <= 100000)
    print(f"Points in Band 1: {np.sum(idx_band1)}")

    sub_xsec = xsec[0, :]
    print(f"Min: {np.min(sub_xsec)}, Max: {np.max(sub_xsec)}")





#xsec_var = nc.variables["kabs"]
#print(xsec_var.dimensions)  # names of dimensions (e.g., ('nu', 'p', 't'))
#print(xsec_var.shape)            # actual size of each dimension (e.g., (3000, 10, 5))
