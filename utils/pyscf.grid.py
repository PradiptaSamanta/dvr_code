#!/usr/bin/env python

import numpy
from pyscf import gto, scf, dft
from pyscf import ao2mo, fci
from pyscf import mcscf
from pyscf.dft import numint
from pyscf.lo import orth
import pyscf.lib
import pyscf.ao2mo
import pyscf.fci



for grid_level in range(0,10):
 mol = gto.Mole()
 mol.build(
     atom =  'H 0 0 0',
     basis = 'cc-pVTZ',
     charge = 0,
     spin   = 1,
     verbose = 0,
 )
 mf = dft.RKS(mol)
 mf.grids.prune = None
# mf.grids.radi_method = dft.gauss_chebyshev
 mf.grids.level = grid_level
 e = mf.scf()
 coords = mf.grids.coords
 weights = mf.grids.weights
# print grid_level, weights.shape
 f = open('3dgrid.g'+str(mf.grids.level),"w+")
 ngrid = weights.shape[0]
# print weights.shape,coords.shape
 ao_value = numint.eval_ao(mol, coords, deriv=0)
 f.write("%i\n" %(ngrid))
 check  = 0 
 check1 = 0 
 check2 = 0 
 for i in range(ngrid):
  (x,y,z) = coords[i,:]
  r = numpy.sqrt(x**2+y**2+z**2)
  check  = check  + numpy.exp(-2*r**2)*weights[i]
  check2 = check2 + 2.71828182845904523536**(-2*abs(r))*weights[i]
  check1 = check1 + ao_value[i,2]**2*weights[i]
  f.write("%22.15e %22.15e %22.15e %22.15e \n" %(coords[i,0],coords[i,1],coords[i,2],weights[i]))
 print check,check1,check2,ngrid
 

