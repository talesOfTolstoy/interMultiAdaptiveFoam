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
if (correctPhi)
	{
		rAU = 1.0/UEqn.A();
	}
 else
	 {
		 rAU = 1.0/UEqn.A();
	 }

surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
surfaceScalarField phiHbyA
(
 "phiHbyA",
 fvc::flux(HbyA)
 + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi,Uf))
 );
MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
	{
		fvc::makeRelative(phiHbyA, U);
		adjustPhi(phiHbyA, U, p_rgh);
		fvc::makeAbsolute(phiHbyA, U);
	}

surfaceScalarField phig
(
 (
	// mixture.surfaceTensionForce()
	brackFS
	- ghf*fvc::snGrad(rho)
	)*rAUf*mesh.magSf()
 );

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

while (piso.correctNonOrthogonal())
	{
		fvScalarMatrix p_rghEqn
			(
			 fvm::laplacian(rAUf, p_rgh)
			 ==
			 fvc::div(phiHbyA)
			 - mDot/rho
			 );

		p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

		p_rghEqn.solve
			(
			 mesh.solver
			 (
				p_rgh.select(bool(corr==(nCorr-1)))
				)
			 );

		if (piso.finalNonOrthogonalIter())
			{
				phi = phiHbyA - p_rghEqn.flux();

				p_rgh.relax();

				U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
				U.correctBoundaryConditions();
				fvOptions.correct(U);
			}
	}

#include "continuityError.H"

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

p == p_rgh + rho*gh;

if (p_rgh.needReference())
	{
		p += dimensionedScalar
			(
			 "p",
			 p.dimensions(),
			 pRefValue - getRefCellValue(p, pRefCell)
			 );
		p_rgh = p - rho*gh;
	}

if (!correctPhi)
	{
		rAU.clear();
	}
