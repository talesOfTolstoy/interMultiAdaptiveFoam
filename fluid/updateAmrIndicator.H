/*---------------------------------------------------------------------------*\
  Paper-inspired AMR indicator.

  The Basilisk reference case adapts on wavelet errors in the VOF field and
  velocity, and explicitly maintains refinement at the solid boundary.  This
  OpenFOAM solver has a scalar dynamicRefineFvMesh selector rather than a
  multi-field wavelet estimator, so the corresponding flags are combined into
  one binary indicator field.
\*---------------------------------------------------------------------------*/

#include <cmath>

{
amrIndicator = scalar(0);

// Equivalent to the paper's uemax = 1e-2, using Delta*|grad(U)| as a local
// finite-volume proxy for the wavelet detail of the velocity field.
tmp<volTensorField> tGradU = fvc::grad(U);
const volTensorField& gradU = tGradU();
const scalarField& cellVolume = mesh.V();

forAll(amrIndicator, celli)
{
    const scalar alphaValue = alpha1[celli];
    const scalar delta = std::cbrt(cellVolume[celli]);
    const scalar velocityError = delta*mag(gradU[celli]);

    // Equivalent to the paper's femax = 1e-6 for a sharp VOF field:
    // mixed cells are the cells carrying the liquid-vapour interface.
    if
    (
        (alphaValue > 1.0e-6 && alphaValue < 1.0 - 1.0e-6)
     || velocityError > 1.0e-2
    )
    {
        amrIndicator[celli] = 1.0;
    }
}

// The paper explicitly refines the solid boundary.  In this two-region case
// the fluid-side owner cells of the mapped wall are the corresponding flag.
forAll(mesh.boundary(), patchi)
{
    const fvPatch& patch = mesh.boundary()[patchi];

    if (patch.name() == "fluid_to_solid")
    {
        const labelUList& faceCells = patch.faceCells();

        forAll(faceCells, facei)
        {
            amrIndicator[faceCells[facei]] = 1.0;
        }
    }
}

amrIndicator.correctBoundaryConditions();
}
