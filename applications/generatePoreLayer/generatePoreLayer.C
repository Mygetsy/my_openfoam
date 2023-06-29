/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    generatePoreLayer

Description
    Set initial conditions for alpha field, which represent porosity in the
     substrate. Works with dynamicRefineFvMesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using constant::mathematical::pi;

// Positive if outside ball, Negative if inside ball
scalarField generateBall
(
    const vectorField& points,
    const vector centre,
    const scalar radius
)
{
    return mag(points - centre) - radius;
}

tmp<volScalarField> VolumeOfFluid(const fvMesh& mesh, scalarField& f)
{
    cutCellIso cutCell(mesh, f);
    auto tres = volScalarField::New("result", mesh, dimensionedScalar());
    auto& res = tres.ref();

    forAll(res, cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI, Zero);

        if (cellStatus == -1)
        {
            res[cellI] = 1;
        }
        else if (cellStatus == 1)
        {
            res[cellI] = 0;
        }
        else if (cellStatus == 0)
        {
            if (mag(cutCell.faceArea()) != 0)
            {
                res[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0);
            }
        }
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addArgument("field");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    word alphaName = args.get<word>(1);

    Info<< "Reading field " << alphaName << endl;
    volScalarField alpha
    (
        IOobject
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const word dictName("poreLayerProperties");
    Info<< "Reading " << dictName << endl;
    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Dimensionless parameters of the pore distribution
    const label seed(dict.lookupOrDefault<label>("seed", 0));
    const label numberOfPores(dict.lookupOrDefault<label>("numberOfPores", 1));
    const scalar amplitudeRadius(dict.lookupOrDefault<scalar>("amplitudeRadius", 0));
    
    const vector poreLayerBoxMin(dict.lookup("poreLayerBoxMin"));
    const vector poreLayerBoxMax(dict.lookup("poreLayerBoxMax"));

    // Dimensioned parameters
    const dimensionedScalar poreRadius("poreRadius", dict);
    const dimensionedScalar substratePosition("substratePosition", dict);

    // Auxiliary constants
    const vector substrateNormal(0, 0, 1);
    const scalar domainVolume = gSum(mesh.V());
    const boundBox& bounds = mesh.bounds();

    label prevMeshSize;
    label timeIndex = 1;

    if (isA<dynamicRefineFvMesh>(mesh))
    {
        dictionary refineDict
        (
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ).optionalSubDict(mesh.typeName + "Coeffs")
        );
        timeIndex = refineDict.get<label>("refineInterval");
    }

    // NB: mesh.update() works only for timeIndex > 0 && timeIndex % refineInterval == 0
    runTime.setTime(0, timeIndex);

    do
    {
        scalar exactVolumeFraction = 0;  // for theoretical prediction

        // --- Generate substrate
        {
            scalarField f = substratePosition.value() - (mesh.points() & substrateNormal);
            alpha = VolumeOfFluid(mesh, f);
            exactVolumeFraction = alpha.weightedAverage(mesh.Vsc()).value();
        }
        // --- Generate pores
        Random random(seed);
        for (int i=0; i<numberOfPores; i++){
            const scalar R = poreRadius.value()*(1 + amplitudeRadius*random.position(-1., 1.));
            const scalar X = 0.5 * (poreLayerBoxMax.x() - poreLayerBoxMin.x())*random.position(-1., 1.) + 0.5 * (poreLayerBoxMax.x() + poreLayerBoxMin.x());
            const scalar Y = 0.5 * (poreLayerBoxMax.y() - poreLayerBoxMin.y())*random.position(-1., 1.) + 0.5 * (poreLayerBoxMax.y() + poreLayerBoxMin.y());
            const scalar Z = 0.5 * (poreLayerBoxMax.z() - poreLayerBoxMin.z())*random.position(-1., 1.) + 0.5 * (poreLayerBoxMax.z() + poreLayerBoxMin.z());
            
            scalarField f = -generateBall(mesh.points(), vector(X, Y, Z), R);
            alpha -= VolumeOfFluid(mesh, f);
            
            if
            (
                bounds.min().x() < X && X < bounds.max().x()
             && bounds.min().y() < Y && Y < bounds.max().y()
            )
            {
                exactVolumeFraction += 4./3*pi*pow(R, 3)/domainVolume;
            }
        }
        
        alpha.clip(0, 1);
        alpha.correctBoundaryConditions();

        // --- Analyze the result
        const scalar volumeFraction = alpha.weightedAverage(mesh.Vsc()).value();
        Info<< nl << "Volume fraction = " << volumeFraction
            << " theoretical = " << exactVolumeFraction
            << " error = " << (volumeFraction - exactVolumeFraction)/exactVolumeFraction
            << nl << endl;

        prevMeshSize = mesh.cells().size();
        mesh.update();
    }
    while (mesh.changing() && mesh.cells().size() > prevMeshSize);

    // --- Save the result
    Info<< "Writing field " << alphaName << endl;
    ISstream::defaultPrecision(18);
    runTime.setTime(0, 0);
    runTime.writeNow();

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
