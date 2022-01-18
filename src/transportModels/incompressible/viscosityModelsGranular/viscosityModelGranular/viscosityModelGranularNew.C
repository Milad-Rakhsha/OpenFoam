/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "viscosityModelGranular.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::viscosityModelGranular> Foam::viscosityModelGranular::New
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& p


)
{
    const word modelType(viscosityProperties.lookup("transportModel"));

    Info<< "Selecting incompressible transport model " << modelType << endl;

//    dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(modelType);
//
//    if (cstrIter == dictionaryConstructorTablePtr_->end())
//    {
//        FatalErrorInFunction
//            << "Unknown viscosityModel type "
//            << modelType << nl << nl
//            << "Valid viscosityModels are : " << endl
//            << dictionaryConstructorTablePtr_->sortedToc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<viscosityModelGranular>
//        (cstrIter()(name, viscosityProperties, U, phi,p));


    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
        	viscosityProperties,
            "viscosityModelGranular",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<viscosityModelGranular>(ctorPtr(name, viscosityProperties, U, phi, p));


}


// ************************************************************************* //
