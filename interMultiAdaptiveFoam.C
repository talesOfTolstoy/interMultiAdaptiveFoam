/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
		Copyright (C) 2026 Borja Diaz-Guardamino Sobrin
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

Application
    interMultiAdaptiveFoam

Description
    Transient solver for multi-region incompressible two-phase interface
		tracking simulations with phase change.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
// #include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fvcSmooth.H"

#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "incompressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
// #include "radiationModel.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "processorFvPatchFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for multi-region two-phase "
				"interface tracking phase change simulations."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "incompressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        // #include "readPISOControls.H"

        #include "incompressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // if (nOuterCorr != 1)
        // {
        //     forAll(fluidRegions, i)
        //     {
        //         #include "storeOldFluidFields.H"
        //     }
        // }

				forAll(fluidRegions, i)
					{
						Info<< "\nSolving for fluid region "
								<< fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionPISOControls.H"
            #include "solveFluid.H"
					}

				forAll(solidRegions, i)
					{
						Info<< "\nSolving for solid region "
								<< solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"
            #include "readSolidMultiRegionPISOControls.H"
            #include "solveSolid.H"
					}
				
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
