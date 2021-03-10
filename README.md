# ddPCM
A fast domain decomposition based implementation of the COSMO solvation model

[![DOI](https://zenodo.org/badge/63929685.svg)](https://zenodo.org/badge/latestdoi/63929685)

[![Bintray](https://img.shields.io/github/workflow/status/Nige91/ddPCM/Build/ddX)](https://img.shields.io/github/workflow/status/Nige91/ddPCM/Build/ddX)

COPYRIGHT (C) 2016 by Filippo Lipparini, Benjamin Stamm, Eric Canc�s,
Yvon Maday, Paolo Gatto, Jean-Philip Piquemal, Louis Lagard�re and 
Benedetta Mennucci.   
                          ALL RIGHT RESERVED.      

This code is governed by the LGPL license and abiding by the rules of 
distribution of free software.  
 
This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU Lesser General Public License for more details.

Users of this code are asked to include the following references in their
publications:

[1] E. Canc�s, Y. Maday, B. Stamm
    "Domain decomposition for implicit solvation models"
    J. Chem. Phys. 139, 054111 (2013)

[2] F. Lipparini, B. Stamm, E. Canc�s, Y. Maday, B. Mennucci
    "Fast Domain Decomposition Algorithm for Continuum Solvation Models: 
     Energy and First Derivatives"
    J. Chem. Theory Comput. 9, 3637�3648 (2013)

Also, include one of the three following reference depending on whether you
use this code in conjunction with a QM [3], Semiempirical [4] or Classical [5]
description of the solute:

[3] F. Lipparini, G. Scalmani, L. Lagard�re, B. Stamm, E. Canc�s, Y. Maday,
    J.-P. Piquemal, M. J. Frisch, B. Mennucci
    "Quantum, classical, and hybrid QM/MM calculations in solution: General 
     implementation of the ddCOSMO linear scaling strategy"
    J. Chem. Phys. 141, 184108 (2014)
    (for quantum mechanical models)

[4] F. Lipparini, L. Lagard�re, G. Scalmani, B. Stamm, E. Canc�s, Y. Maday,
    J.-P. Piquemal, M. J. Frisch, B. Mennucci
    "Quantum Calculations in Solution for Large to Very Large Molecules: 
     A New Linear Scaling QM/Continuum Approach"
    J. Phys. Chem. Lett. 5, 953-958 (2014)
    (for semiempirical models)

[5] F. Lipparini, L. Lagard�re, C. Raynaud, B. Stamm, E. Canc�s, B. Mennucci
    M. Schnieders, P. Ren, Y. Maday, J.-P. Piquemal
    "Polarizable Molecular Dynamics in a Polarizable Continuum Solvent"
    J. Chem. Theory Comput. 11, 623-634 (2015)
    (for classical models, including polarizable force fields

The users of this code should also include the appropriate reference to the
COSMO model. This distribution includes the routines to generate lebedev
grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the routines
are used, the following reference should also be included:

[6] V.I. Lebedev, and D.N. Laikov
    "A quadrature formula for the sphere of the 131st
     algebraic order of accuracy"
    Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.

Written by Filippo Lipparini, October 2015.

Please report problems and bug to filippo.lipparini@gmail.com

## build