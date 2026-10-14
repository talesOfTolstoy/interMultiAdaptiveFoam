/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2011 OpenFOAM Foundation
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Description
    Custom axisymmetric refinement with a physical-size candidate filter.
    Assumes the rotation axis is z through x = y = 0.
    targetCellSize is required in dynamicRefineFvMeshAxiCoeffs, in metres.
    Existing consistency, flux mapping and unrefinement logic is retained.
    This source has not been compiled or runtime-tested by its editor.
\*---------------------------------------------------------------------------*/

#include "dynamicRefineFvMeshAxi.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "polyMesh.H"
#include "meshTools.H"

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineFvMeshAxi, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineFvMeshAxi, IOobject);
}

Foam::label Foam::dynamicRefineFvMeshAxi::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val) { n++; }
    }
    return n;
}

void Foam::dynamicRefineFvMeshAxi::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    if (protectedCell_.empty())
    {
        unrefineableCell.clear();
        return;
    }
    const labelList& cellLevel = meshCutter_.cellLevel();
    unrefineableCell = protectedCell_;
    labelList neiLevel(nFaces()-nInternalFaces());
    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        neiLevel[faceI-nInternalFaces()] = cellLevel[faceOwner()[faceI]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);
    while (true)
    {
        boolList seedFace(nFaces(), false);
        forAll(faceNeighbour(), faceI)
        {
            label own = faceOwner()[faceI];
            label nei = faceNeighbour()[faceI];
            bool ownProtected = unrefineableCell.get(own);
            bool neiProtected = unrefineableCell.get(nei);
            if (ownProtected && cellLevel[nei] > cellLevel[own])
            {
                seedFace[faceI] = true;
            }
            else if (neiProtected && cellLevel[own] > cellLevel[nei])
            {
                seedFace[faceI] = true;
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            label own = faceOwner()[faceI];
            if
            (
                unrefineableCell.get(own)
             && neiLevel[faceI-nInternalFaces()] > cellLevel[own]
            )
            {
                seedFace[faceI] = true;
            }
        }
        syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());
        bool hasExtended = false;
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            if (seedFace[faceI])
            {
                label own = faceOwner()[faceI];
                label nei = faceNeighbour()[faceI];
                if (!unrefineableCell.get(own))
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
                if (!unrefineableCell.get(nei))
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            if (seedFace[faceI])
            {
                label own = faceOwner()[faceI];
                if (!unrefineableCell.get(own))
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
            }
        }
        if (!returnReduce(hasExtended, orOp<bool>())) { break; }
    }
}

void Foam::dynamicRefineFvMeshAxi::readDict()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict", time().constant(), *this,
                IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE, false
            )
        ).subDict(typeName + "Coeffs")
    );
    List<Pair<word> > fluxVelocities(refineDict.lookup("correctFluxes"));
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }
    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}

Foam::autoPtr<Foam::mapPolyMesh> Foam::dynamicRefineFvMeshAxi::refine
(
    const labelList& cellsToRefine
)
{
    polyTopoChange meshMod(*this);
    meshCutter_.setRefinement(cellsToRefine, meshMod);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);
    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;
    if (debug)
    {
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            label oldFaceI = map().faceMap()[faceI];
            if (oldFaceI >= nInternalFaces())
            {
                FatalErrorIn("dynamicRefineFvMeshAxi::refine(const labelList&)")
                    << "New internal face:" << faceI
                    << " fc:" << faceCentres()[faceI]
                    << " originates from boundary oldFace:" << oldFaceI
                    << abort(FatalError);
            }
        }
    }
    updateMesh(map);
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();
        labelHashSet masterFaces(4*cellsToRefine.size());
        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];
            if (oldFaceI >= 0)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];
                if (masterFaceI < 0)
                {
                    FatalErrorIn("dynamicRefineFvMeshAxi::refine(const labelList&)")
                        << "Problem: should not have removed faces when refining."
                        << nl << "face:" << faceI << abort(FatalError);
                }
                else if (masterFaceI != faceI)
                {
                    masterFaces.insert(masterFaceI);
                }
            }
        }
        if (debug)
        {
            Info<< "Found " << returnReduce(masterFaces.size(), sumOp<label>())
                << " split faces " << endl;
        }
        HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
        forAllConstIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningIn("dynamicRefineFvMeshAxi::refine(const labelList&)")
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning." << endl;
                continue;
            }
            const word& UName = correctFluxes_[iter.key()];
            if (UName == "none") { continue; }
            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName << endl;
            }
            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf()
            );
            for (label faceI = 0; faceI < nInternalFaces(); faceI++)
            {
                label oldFaceI = faceMap[faceI];
                if (oldFaceI == -1 || reverseFaceMap[oldFaceI] != faceI)
                {
                    phi[faceI] = phiU[faceI];
                }
            }
            surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();
            forAll(bphi, patchI)
            {
                fvsPatchScalarField& patchPhi = bphi[patchI];
                const fvsPatchScalarField& patchPhiU = phiU.boundaryField()[patchI];
                label faceI = patchPhi.patch().start();
                forAll(patchPhi, i)
                {
                    label oldFaceI = faceMap[faceI];
                    if (oldFaceI == -1 || reverseFaceMap[oldFaceI] != faceI)
                    {
                        patchPhi[i] = patchPhiU[i];
                    }
                    faceI++;
                }
            }
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label faceI = iter.key();
                if (isInternalFace(faceI)) { phi[faceI] = phiU[faceI]; }
                else
                {
                    label patchI = boundaryMesh().whichPatch(faceI);
                    label i = faceI - boundaryMesh()[patchI].start();
                    bphi[patchI][i] = phiU.boundaryField()[patchI][i];
                }
            }
        }
    }
    meshCutter_.updateMesh(map);
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());
        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
        }
        protectedCell_.transfer(newProtectedCell);
    }
    meshCutter_.checkRefinementLevels(-1, labelList(0));
    return map;
}

Foam::autoPtr<Foam::mapPolyMesh> Foam::dynamicRefineFvMeshAxi::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(*this);
    meshCutter_.setUnrefinement(splitPoints, meshMod);
    Map<label> faceToSplitPoint(3*splitPoints.size());
    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];
        const labelList& pEdges = pointEdges()[pointI];
        forAll(pEdges, j)
        {
            label otherPointI = edges()[pEdges[j]].otherVertex(pointI);
            vector coord2 = points()[otherPointI];
            const labelList& pFaces1 = pointFaces()[pointI];
            forAll(pFaces1, pFaceI)
            {
                faceToSplitPoint.insert(pFaces1[pFaceI], pointI);
            }
            // Original coordinate-dependent logic retained, not validated
            // for the current wedge orientation. Keep unrefinement disabled.
            if (coord2[2] > 0)
            {
                const labelList& pFaces2 = pointFaces()[otherPointI];
                forAll(pFaces2, pFaceJ)
                {
                    faceToSplitPoint.insert(pFaces2[pFaceJ], otherPointI);
                }
            }
        }
    }
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);
    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;
    updateMesh(map);
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();
        HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
        forAllConstIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningIn("dynamicRefineFvMeshAxi::unrefine(const labelList&)")
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning." << endl;
                continue;
            }
            const word& UName = correctFluxes_[iter.key()];
            if (UName == "none") { continue; }
            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName << endl;
            }
            surfaceScalarField& phi = *iter();
            surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();
            const surfaceScalarField phiU
            (
                fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf()
            );
            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFaceI = iter.key();
                label oldPointI = iter();
                if (reversePointMap[oldPointI] < 0)
                {
                    label faceI = reverseFaceMap[oldFaceI];
                    if (faceI >= 0)
                    {
                        if (isInternalFace(faceI)) { phi[faceI] = phiU[faceI]; }
                        else
                        {
                            label patchI = boundaryMesh().whichPatch(faceI);
                            if (!isA<emptyPolyPatch>(boundaryMesh()[patchI]))
                            {
                                label i = faceI - boundaryMesh()[patchI].start();
                                bphi[patchI][i] = phiU.boundaryField()[patchI][i];
                            }
                        }
                    }
                }
            }
        }
    }
    meshCutter_.updateMesh(map);
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());
        forAll(newProtectedCell, cellI)
        {
            label oldCellI = map().cellMap()[cellI];
            if (oldCellI >= 0)
            {
                newProtectedCell.set(cellI, protectedCell_.get(oldCellI));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }
    meshCutter_.checkRefinementLevels(-1, labelList(0));
    return map;
}

Foam::scalarField Foam::dynamicRefineFvMeshAxi::maxPointField
(
    const scalarField& pFld
) const
{
    scalarField vFld(nCells(), -GREAT);
    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];
        forAll(pCells, i) { vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointI]); }
    }
    return vFld;
}

Foam::scalarField Foam::dynamicRefineFvMeshAxi::minCellField
(
    const volScalarField& vFld
) const
{
    scalarField pFld(nPoints(), GREAT);
    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];
        forAll(pCells, i) { pFld[pointI] = min(pFld[pointI], vFld[pCells[i]]); }
    }
    return pFld;
}

Foam::scalarField Foam::dynamicRefineFvMeshAxi::cellToPoint
(
    const scalarField& vFld
) const
{
    scalarField pFld(nPoints());
    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];
        scalar sum = 0.0;
        forAll(pCells, i) { sum += vFld[pCells[i]]; }
        pFld[pointI] = sum/pCells.size();
    }
    return pFld;
}

Foam::scalarField Foam::dynamicRefineFvMeshAxi::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);
    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);
        if (err >= 0) { c[i] = err; }
    }
    return c;
}

void Foam::dynamicRefineFvMeshAxi::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    scalarField cellError
    (
        maxPointField(error(cellToPoint(vFld), lowerRefineLevel, upperRefineLevel))
    );
    forAll(cellError, cellI)
    {
        if (cellError[cellI] > 0) { candidateCell.set(cellI, 1); }
    }
}

Foam::labelList Foam::dynamicRefineFvMeshAxi::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Original conservative budget estimate retained. The cutter actually
    // creates four children; consistency expansion can exceed this budget.
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;
    const labelList& cellLevel = meshCutter_.cellLevel();
    PackedBoolList unrefineableCell;
    calculateProtectedCells(unrefineableCell);
    label nCandidates = returnReduce(count(candidateCell, 1), sumOp<label>());
    DynamicList<label> candidates(nCells());
    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, cellI)
        {
            if
            (
                cellLevel[cellI] < maxRefinement
             && candidateCell.get(cellI)
             && (unrefineableCell.empty() || !unrefineableCell.get(cellI))
            )
            {
                candidates.append(cellI);
            }
        }
    }
    else
    {
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCell, cellI)
            {
                if
                (
                    cellLevel[cellI] == level
                 && candidateCell.get(cellI)
                 && (unrefineableCell.empty() || !unrefineableCell.get(cellI))
                )
                {
                    candidates.append(cellI);
                }
            }
            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }
    labelList consistentSet
    (
        meshCutter_.consistentRefinement(candidates.shrink(), true)
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;
    return consistentSet;
}

Foam::labelList Foam::dynamicRefineFvMeshAxi::selectUnrefinePoints
(
    const label axis,
    const scalar axisVal,
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    const labelList splitPoints(meshCutter_.getSplitPoints(axis, axisVal));
    DynamicList<label> newSplitPoints(splitPoints.size());
    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];
        if (pFld[pointI] < unrefineLevel)
        {
            const labelList& pCells = pointCells()[pointI];
            bool hasMarked = false;
            forAll(pCells, pCellI)
            {
                if (markedCell.get(pCells[pCellI]))
                {
                    hasMarked = true;
                    break;
                }
            }
            if (!hasMarked) { newSplitPoints.append(pointI); }
        }
    }
    newSplitPoints.shrink();
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement(newSplitPoints, false)
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>()) << "." << endl;
    return consistentSet;
}

void Foam::dynamicRefineFvMeshAxi::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    boolList markedFace(nFaces(), false);
    forAll(markedCell, cellI)
    {
        if (markedCell.get(cellI))
        {
            const cell& cFaces = cells()[cellI];
            forAll(cFaces, i) { markedFace[cFaces[i]] = true; }
        }
    }
    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());
    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
            markedCell.set(faceNeighbour()[faceI], 1);
        }
    }
    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        if (markedFace[faceI]) { markedCell.set(faceOwner()[faceI], 1); }
    }
}

Foam::dynamicRefineFvMeshAxi::dynamicRefineFvMeshAxi(const IOobject& io)
:
    dynamicFvMesh(io),
    meshCutter_(*this),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCell_(nCells(), 0)
{
    readDict();
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    labelList nAnchors(nCells(), 0);
    label nProtected = 0;
    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];
        forAll(pCells, i)
        {
            label cellI = pCells[i];
            if (!protectedCell_.get(cellI) && pointLevel[pointI] <= cellLevel[cellI])
            {
                nAnchors[cellI]++;
                if (nAnchors[cellI] > 8)
                {
                    protectedCell_.set(cellI, 1);
                    nProtected++;
                }
            }
        }
    }
    {
        labelList neiLevel(nFaces());
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceNeighbour()[faceI]];
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[faceOwner()[faceI]];
        }
        syncTools::swapFaceList(*this, neiLevel);
        boolList protectedFace(nFaces(), false);
        forAll(faceOwner(), faceI)
        {
            label faceLevel = max(cellLevel[faceOwner()[faceI]], neiLevel[faceI]);
            const face& f = faces()[faceI];
            label nAnchors = 0;
            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;
                    if (nAnchors > 4)
                    {
                        protectedFace[faceI] = true;
                        break;
                    }
                }
            }
        }
        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
                protectedCell_.set(faceNeighbour()[faceI], 1);
                nProtected++;
            }
        }
        for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(faceOwner()[faceI], 1);
                nProtected++;
            }
        }
    }
    if (returnReduce(nProtected, sumOp<label>()) == 0) { protectedCell_.clear(); }
}

Foam::dynamicRefineFvMeshAxi::~dynamicRefineFvMeshAxi()
{}

bool Foam::dynamicRefineFvMeshAxi::update()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict", time().constant(), *this,
                IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE, false
            )
        ).subDict(typeName + "Coeffs")
    );
    label refineInterval = readLabel(refineDict.lookup("refineInterval"));
    bool hasChanged = false;
    if (refineInterval == 0)
    {
        topoChanging(hasChanged);
        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorIn("dynamicRefineFvMeshAxi::update()")
            << "Illegal refineInterval " << refineInterval << nl
            << "refineInterval must be zero or positive." << nl
            << exit(FatalError);
    }
    // Existing behaviour: no refinement at time index zero.
    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        label maxCells = readLabel(refineDict.lookup("maxCells"));
        if (maxCells <= 0)
        {
            FatalErrorIn("dynamicRefineFvMeshAxi::update()")
                << "maxCells must be positive." << exit(FatalError);
        }
        label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));
        if (maxRefinement <= 0)
        {
            FatalErrorIn("dynamicRefineFvMeshAxi::update()")
                << "maxRefinement must be positive." << exit(FatalError);
        }
        const word fieldName(refineDict.lookup("field"));
        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);
        const scalar lowerRefineLevel = readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel = readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar targetCellSize = readScalar(refineDict.lookup("targetCellSize"));
        if (targetCellSize <= 0)
        {
            FatalErrorIn("dynamicRefineFvMeshAxi::update()")
                << "targetCellSize must be positive. Received "
                << targetCellSize << nl << exit(FatalError);
        }
        const scalar unrefineLevel = readScalar(refineDict.lookup("unrefineLevel"));
        const label nBufferLayers = readLabel(refineDict.lookup("nBufferLayers"));
        const label nBufferLayersR = readLabel(refineDict.lookup("nBufferLayersR"));
        const label axis = readLabel(refineDict.lookup("axis"));
        const scalar axisVal = readScalar(refineDict.lookup("axisVal"));
        PackedBoolList refineCell(nCells());
        if (globalData().nTotalCells() < maxCells)
        {
            selectRefineCandidates(lowerRefineLevel, upperRefineLevel, vFld, refineCell);
            for (label i = 0; i < nBufferLayersR; i++)
            {
                extendMarkedCells(refineCell);
            }

            // Preserve interface/buffer markings for unrefinement protection.
            // Only the separate refinement candidate list is size-filtered.
            PackedBoolList sizeCandidates(refineCell);
            const pointField& meshPoints = points();
            const scalar sizeLimit = targetCellSize*(1.0 + 1e-8);
            forAll(sizeCandidates, cellI)
            {
                if (!sizeCandidates.get(cellI)) { continue; }
                const labelList& cPoints = cellPoints()[cellI];
                scalar rMin = GREAT;
                scalar rMax = -GREAT;
                scalar zMin = GREAT;
                scalar zMax = -GREAT;
                forAll(cPoints, pointIndex)
                {
                    const point& p = meshPoints[cPoints[pointIndex]];
                    // Rotation axis: z, through x = y = 0.
                    const scalar radius = Foam::sqrt(p.x()*p.x() + p.y()*p.y());
                    rMin = min(rMin, radius);
                    rMax = max(rMax, radius);
                    zMin = min(zMin, p.z());
                    zMax = max(zMax, p.z());
                }
                if (rMax-rMin <= sizeLimit && zMax-zMin <= sizeLimit)
                {
                    sizeCandidates.unset(cellI);
                }
            }
            labelList cellsToRefine
            (
                selectRefineCells(maxCells, maxRefinement, sizeCandidates)
            );
            label nCellsToRefine = returnReduce(cellsToRefine.size(), sumOp<label>());
            if (nCellsToRefine > 0)
            {
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();
                    PackedBoolList newRefineCell(cellMap.size());
                    forAll(cellMap, cellI)
                    {
                        label oldCellI = cellMap[cellI];
                        if (oldCellI < 0 || reverseCellMap[oldCellI] != cellI)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else
                        {
                            newRefineCell.set(cellI, refineCell.get(oldCellI));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMarkedCells(refineCell);
                }
                hasChanged = true;
            }
        }
        {
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    axis, axisVal, unrefineLevel, refineCell, minCellField(vFld)
                )
            );
            label nSplitPoints = returnReduce(pointsToUnrefine.size(), sumOp<label>());
            if (nSplitPoints > 0)
            {
                unrefine(pointsToUnrefine);
                hasChanged = true;
            }
        }
        if ((nRefinementIterations_ % 10) == 0)
        {
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
    }
    topoChanging(hasChanged);
    return hasChanged;
}

bool Foam::dynamicRefineFvMeshAxi::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    const_cast<hexRef8axi&>(meshCutter_).setInstance(time().timeName());
    bool writeOk =
    (
        dynamicFvMesh::writeObject(streamOpt, writeOnProc)
     && meshCutter_.write()
    );
    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel", time().timeName(), *this,
                IOobject::NO_READ, IOobject::AUTO_WRITE, false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );
        const labelList& cellLevel = meshCutter_.cellLevel();
        forAll(cellLevel, cellI) { scalarCellLevel[cellI] = cellLevel[cellI]; }
        writeOk = writeOk && scalarCellLevel.write();
    }
    return writeOk;
}

// ************************************************************************* //

