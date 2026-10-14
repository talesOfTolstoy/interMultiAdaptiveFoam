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
volScalarField SArea = mag(fvc::grad(alpha1));

dimensionedScalar hi = (2.0*sigmaEvap)/(2.0-sigmaEvap)
  *(hlv * hlv)/((Tsat*(1.0/rho2-1.0/rho1)))
  *pow(1.0/(2.0*(constant::mathematical::pi)*R_g*Tsat),0.5);

// dimensionedScalar hi = (2.0*sigmaEvap)/(2.0-sigmaEvap)
// 	*(hlv * hlv)/((Tsat*(1.0/rho2)))
// 	*pow(1.0/(2.0*(constant::mathematical::pi)*R_g*Tsat),0.5);

kappa = alpha1*k1 + alpha2*k2;
surfaceScalarField kappaf(fvc::interpolate(kappa));

volScalarField rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2;
  
surfaceScalarField rhoPhiCp
("rhoPhiCp",
 fvc::interpolate(rhoCp)*phi);

 
qDot = min(alpha1*hi*(T - Tsat)*SArea, maxQdot);

volScalarField corrHlv = qDot / hlv * (cp1-cp2) * (T-Tsat);

for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
	{

		fvScalarMatrix TEqn
			(
			 fvm::ddt(rhoCp,T)
			 + fvm::div(rhoPhiCp,T)
			 - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
			 - fvm::laplacian(kappaf,T)
			 ==
			 - qDot
			 - corrHlv
			 );

		TEqn.relax();
		TEqn.solve();

		//update continuity source with new T
	}

qDot = min(alpha1*hi*(T - Tsat)*SArea, maxQdot);

Info << "Min(T) = " << min(T).value() << " Max(T) = " << max(T).value() << endl;
