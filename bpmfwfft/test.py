# !/usr/bin/env python
# coding: utf-8

# In[1]:


# jupyter nbconvert --to python extract_info_from_YANK.ipynb

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align

import os
import numpy as np
import netCDF4 as netcdf

import time

# In[2]:


# # Extracts information from YANK output
#
# From args.YANK_nc,
# * trajectories of the bound and unbound state as AMBER nc files
# * reduced energies of the bound and unbound state as npy files

# The atom selection that will be aligned to the first frame
sel_str_align = 'protein and name CA'
# Default location of data
base_path = '/home/jim/BPMFwFFT/bpmfwfft/'

import argparse

parser = argparse.ArgumentParser(description='Extracts information from an end state of a YANK simulation')
parser.add_argument('--thermo_state', choices=['bound', 'unbound'], default='bound',
                    help='Thermodynamic state to extract')
parser.add_argument('--prmtop', default=base_path + 'complex.prmtop', help='AMBER parameters and topology')
parser.add_argument('--inpcrd', default=base_path + 'complex.inpcrd', help='AMBER coordinates')
parser.add_argument('--YANK_nc', default=base_path + 'trajectory.nc', help='YANK output')
parser.add_argument('-f', help='Dummy argument for jupyter-notebook')
args = parser.parse_args()


# In[3]:


def init_AMBER_nc(FN, natoms):
    """Initializes an AMBER-format netcdf file
    """
    # Adapted from
    # https://www.mdanalysis.org/mdanalysis/_modules/MDAnalysis/coordinates/TRJ.html#NCDFWriter

    # Set global attributes
    import netCDF4
    nc = netCDF4.Dataset(FN, 'w', format='NETCDF3_64BIT')

    for attribute in ['application', 'program', 'programVersion']:
        setattr(nc, attribute, getattr(nc_in, attribute))
    setattr(nc, 'Conventions', 'AMBER')
    setattr(nc, 'ConventionVersion', '1.0')

    # Create dimensions
    nc.createDimension('frame',
                       None)  # unlimited number of steps (can append)
    nc.createDimension('atom', natoms)  # number of atoms in system
    nc.createDimension('spatial', 3)  # number of spatial dimensions
    nc.createDimension('cell_spatial', 3)  # unitcell lengths
    nc.createDimension('cell_angular', 3)  # unitcell angles
    nc.createDimension('label', 5)  # needed for cell_angular

    # Create variables.
    coords = nc.createVariable('coordinates', 'f4',
                               ('frame', 'atom', 'spatial'))
    setattr(coords, 'units', 'angstrom')

    spatial = nc.createVariable('spatial', 'c', ('spatial',))
    spatial[:] = np.asarray(list('xyz'))

    time = nc.createVariable('time', 'f4', ('frame',))
    setattr(time, 'units', 'picosecond')

    nc.sync()
    return nc


# In[4]:


# Load netCDF4 file from YANK output
#nc_in = netcdf.Dataset(args.YANK_nc)
#nits = nc_in.variables['positions'].shape[0]
#nstates = nc_in.variables['states'].shape[1]
#natoms = nc_in.variables['positions'].shape[2]
"""
# Determine state index and replica indices
state_ind = 0 if args.thermo_state == 'bound' else nstates - 1
(its, rep_inds) = np.nonzero(nc_in.variables['states'][:, :] == state_ind)
"""
# Extract energies
energy_FN = args.thermo_state + '_energies.npy'
if not os.path.isfile(energy_FN):
    print(f'Extracting energies from {args.thermo_state} state, just kidding doing nothing')
    start_time = time.time()
   # flat_inds = np.ravel_multi_index((its, rep_inds, state_ind), dims=nc_in.variables['energies'].shape)
   # energies = np.array(nc_in.variables['energies'][:, :, :].take(flat_inds))
   # np.save(energy_FN, energies)
    print(f'Completed in {time.time() - start_time} s')
    print('hello?')

unaligned_nc_FN = 'trajectory2.nc'
aligned_nc_FN = 'trajectory.nc'
if True:
    # Extract coordinates
    if True:
        print(f'Extracting coordinates from {args.thermo_state} state')
        start_time = time.time()
        #reshaped_coords = (nc_in.variables['positions'][:, :, :, :] * 10.).reshape(nits * nstates, natoms, 3)
        #flat_inds = np.ravel_multi_index((its, rep_inds), dims=(nits, nstates))
        #traj = np.array(reshaped_coords[flat_inds, :, :])
        #nc_out = init_AMBER_nc(unaligned_nc_FN, natoms)
        #nc_out.variables['coordinates'][:, :] = traj
        #nc_out.close()
        print(f'Completed in {time.time() - start_time} s')

    # Align coordinates
    if True:
        import sys
        print(f'Aligning trajectory to reference frame')
        nc_in = netcdf.Dataset('trajectory-ssh.nc')
        print(nc_in.variables)
        ref = mda.Universe(args.prmtop, args.inpcrd)
        nc_in.variables['positions'][:, :, :] = ref.coord.positions
        nc_in.close()
        #mobile = mda.Universe(args.prmtop, unaligned_nc_FN, in_memory=True)
        #alignment = align.AlignTraj(mobile, ref, select='protein and name CA')
        #alignment.run()  # Because the trajectory is in memory the coordinates are modified in place
        """
        w = mda.coordinates.TRJ.NCDFWriter(aligned_nc_FN, ref.atoms.n_atoms)
        w.write(ref)
        print(ref.coord.positions)
        print(w)
        w.trjfile.variables['positions'][:, :, :] = ref.coord.positions
        w.close()
        """