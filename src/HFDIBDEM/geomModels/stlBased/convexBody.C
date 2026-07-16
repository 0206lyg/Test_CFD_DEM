/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "convexBody.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "treeBoundBox.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
// create immersed body for convex body
void convexBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{    
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    // find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    if(!isBBoxInMesh())
    {
        return;
    }

    label cellInIB = getCellInBody(octreeField);
    if(cellInIB == -1)
    {
        return;
    }

    // get the list of cell centroids
    const pointField& cp = mesh_.C();

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellInIB));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label tableSize = 128;
    if(cachedNeighbours_.valid()  && getRefineBuffers() == 0)
    {
        tableSize = cachedNeighbours_().toc().size()*1.5;
    }
    else
    {
        cachedNeighbours_ = new HashTable<const labelList&, label, Hash<label>>;
    }

    HashTable<bool, label, Hash<label>> cellInside(tableSize);

    label iterCount(0);label iterMax(mesh_.nCells());
    while (nextToCheck().size() > 0 and iterCount++ < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            label cCell = nextToCheck()[cellToCheck];
            if (!cellInside.found(cCell))
            {
                iterCount++;

                if(pointInside(cp[cCell]))
                {
                    cellInside.set(cCell, true);

                    if(cachedNeighbours_.valid() && cachedNeighbours_().found(cCell))
                    {
                        auxToCheck().append(cachedNeighbours_()[cCell]);
                    }
                    else
                    {
                        const labelList& neigh = mesh_.cellCells(cCell);
                        cachedNeighbours_().insert(cCell, neigh);
                        auxToCheck().append(neigh);
                    }
                }
                else
                {
                    cellInside.set(cCell, false);
                }
            }
        }
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }

    DynamicLabelList keyToErase;
    for(auto it = cachedNeighbours_().begin(); it != cachedNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    cachedNeighbours_().erase(keyToErase);

    DynamicLabelList potentSurfCells = 
        getPotentSurfCells(
            body,
            cellInside,
            cellPoints
        );

    correctSurfCells
    (
        body,
        potentSurfCells,
        cellInside,
        cellPoints
    );

    if(intCells_
    [Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = min(intCells_
        [Pstream::myProcNo()]);
    }
}
//---------------------------------------------------------------------------//
//
// Create immersed body on a moving/deformed mesh.
//
// Unlike the legacy reconstruction, this path does not require
// a cell-centre-inside seed and does not use a centroid-connected flood-fill.
// Candidate cells are obtained directly from the current mesh cell octree.
//
// The existing centre/vertex sampling in correctSurfCells() is retained.
//
void convexBody::createImmersedBodyDynamic
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    const label procI = Pstream::myProcNo();

    // octreeField is retained in the interface for consistency with
    // the legacy reconstruction, but it is not needed in this branch.
    (void)octreeField;

    // Clear reconstruction data from the previous mesh configuration.
    intCells_[procI].clear();
    surfCells_[procI].clear();

    // This quantity is used later to select the processor containing
    // the largest part of the body.
    ibPartialVolume_[procI] = 0;

    // Current transformed body bounding box.
    const boundBox bodyBounds(getBounds());

    // The OpenFOAM cell tree is invalidated automatically by movePoints()
    // and rebuilt on first access after mesh motion.
    const treeBoundBox searchBox
    (
        bodyBounds.min(),
        bodyBounds.max()
    );

    // Broad-phase search using the current deformed mesh geometry.
    const labelList candidateCells
    (
        mesh_.cellTree().findBox(searchBox)
    );

    label cachedCandidatesLocal = candidateCells.size();

    // correctSurfCells() requires the centre-inside status of every
    // cell supplied in potentSurfCells.
    const label tableSize =
        max(label(128), 2*candidateCells.size());

    HashTable<bool, label, Hash<label>> cellInside(tableSize);

    DynamicLabelList potentSurfCells(candidateCells.size());

    const pointField& cellCentres = mesh_.C();

    forAll(candidateCells, candidateI)
    {
        const label cCell = candidateCells[candidateI];

        cellInside.insert
        (
            cCell,
            pointInside(cellCentres[cCell])
        );

        // Every cell whose AABB overlaps the body AABB is examined
        // using the existing centre/vertex lambda approximation.
        potentSurfCells.append(cCell);
    }

    // Retain the original lambda calculation:
    //
    // cBody =
    //     0.5 * centreInside
    //   + 0.5/Nvertices * numberOfVerticesInside
    //
    // This function also fills intCells_, surfCells_,
    // ibPartialVolume_, and body[cell].
    correctSurfCells
    (
        body,
        potentSurfCells,
        cellInside,
        cellPoints
    );

// Number of reconstructed cells on this processor
label positiveLocal =
    intCells_[Pstream::myProcNo()].size()
  + surfCells_[Pstream::myProcNo()].size();

label cachedCandidatesGlobal = cachedCandidatesLocal;
label positiveGlobal = positiveLocal;

reduce(cachedCandidatesGlobal, sumOp<label>());
reduce(positiveGlobal, sumOp<label>());


// Print only when reconstruction has completely failed globally.
// Under normal operation this block produces no output.
if (positiveGlobal == 0)
{
    // Independent search object constructed from the current mesh state.
    // It is created only on a failure event, so its cost is negligible.
    meshSearch freshSearch(mesh_);

    const labelList freshCandidates
    (
        freshSearch.cellTree().findBox(searchBox)
    );

    label freshCandidatesGlobal = freshCandidates.size();
    reduce(freshCandidatesGlobal, sumOp<label>());

    // Is the particle centre actually inside a current mesh cell?
    const label centreCellLocal =
        freshSearch.findCell(getCoM());

    bool centreInsideGlobal = (centreCellLocal >= 0);
    reduce(centreInsideGlobal, orOp<bool>());

    const point bbMid =
        0.5*(bodyBounds.min() + bodyBounds.max());

    const scalar bbCentreError =
        mag(bbMid - getCoM());

//    if (Pstream::master())
//    {
//        Info
//            << "DYN_RECON_FAIL"
//            << " time=" << mesh_.time().value()
//            << " CoM=" << getCoM()
//            << " bbMin=" << bodyBounds.min()
//            << " bbMax=" << bodyBounds.max()
//            << " bbCentreError=" << bbCentreError
//            << " cachedCandidates=" << cachedCandidatesGlobal
//            << " freshCandidates=" << freshCandidatesGlobal
//            << " centreInsideMesh=" << centreInsideGlobal
//            << " positiveCells=" << positiveGlobal
//            << endl;
//    }
}

    // Keep a useful starting cell for any later legacy/static
    // reconstruction of this body.
    if (intCells_[procI].size() > 0)
    {
        cellToStartInCreateIB_ = intCells_[procI][0];
    }
    else if (surfCells_[procI].size() > 0)
    {
        cellToStartInCreateIB_ = surfCells_[procI][0];
    }
}

//---------------------------------------------------------------------------//
// Find first cell with center inside the body
label convexBody::getCellInBody
(
    Field<label>& octreeField
)
{
    // octreeField *= 0;
    labelHashSet checkedCells;
    // get the list of cell centroids
    const pointField& cp = mesh_.C();

    if(cellToStartInCreateIB_ >= octreeField.size())
        cellToStartInCreateIB_ = 0;

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellToStartInCreateIB_));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label iterCount(0);label iterMax(mesh_.nCells());

    while (nextToCheck().size() > 0 and iterCount < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            if (!checkedCells.found(nextToCheck()[cellToCheck]))
            {
                checkedCells.insert(nextToCheck()[cellToCheck]);
                iterCount++;

                if(pointInside(cp[nextToCheck()[cellToCheck]]))
                {
                    return nextToCheck()[cellToCheck];
                }
                else
                {
                    auxToCheck().append(mesh_.cellCells()[nextToCheck()[cellToCheck]]);
                }
            }
        }
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return -1;
}
//---------------------------------------------------------------------------//
