/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of interMultiAdaptiveFoam.

    interMultiAdaptiveFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    interMultiAdaptiveFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with interMultiAdaptiveFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
// Update the paper-inspired AMR indicator before the dynamic mesh update.
#include "updateAmrIndicator.H"

// Update the dynamic mesh before assembling the fluid equations.  This is
// also where dynamicRefineFvMesh performs its refinement/unrefinement.
mesh.update();

if (mesh.changing())
  {
    // Do not apply a compression-flux correction from the old topology.
    if (mesh.topoChanging())
      {
        talphaPhi1Corr0.clear();
      }

    // Recalculate mesh-dependent gravity fields after the cells/faces move.
    gh = (g & mesh.C()) - ghRef;
    ghf = (g & mesh.Cf()) - ghRef;

    MRF.update();

    if (correctPhi)
      {
        // Reconstruct the absolute flux from the mapped face velocity.
        phi = mesh.Sf() & Uf();

        #include "correctPhi.H"

        // The solver equations use fluxes relative to mesh motion.
        fvc::makeRelative(phi, U);

        mixture.correct();
      }
  }

// Recalculate on the current topology so the written indicator and the next
// adaptation decision are consistent after refinement or unrefinement.
#include "updateAmrIndicator.H"

// if (frozenFlow)
//   {
// 		continue;
//   }

mixture.correct();

#include "alphaControls.H"
#include "alphaEqnSubCycle.H"

mixture.correct();

// Rebuild the capillary pressure flux after alpha and any mesh topology
// changes. It is intentionally excluded from automatic flux mapping below.
brackFS = fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1);

#include "calcF.H"

#include "TEqn.H"

#include "smoothQdot.H"

#include "UEqn.H"

// --- pressure loop
for (int corr=0; corr<nCorr; corr++)
	{
    #include "pEqn.H"
	}

turbulence.correct();
