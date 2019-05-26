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

#include "singlePhaseTransportModelGranular.H"
#include "../viscosityModelsGranular/viscosityModelGranular/viscosityModelGranular.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseTransportModelGranular, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModelGranular::singlePhaseTransportModelGranular
(
    const volVectorField& U,
    const surfaceScalarField& phi,
	const volScalarField& p

)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    viscosityModelPtr_(viscosityModelGranular::New("nu", *this, U, phi, p ))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModelGranular::~singlePhaseTransportModelGranular()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseTransportModelGranular::nu() const
{
    return viscosityModelPtr_->nu();
}


//Foam::tmp<Foam::volScalarField>
//Foam::singlePhaseTransportModelGranular::strainRate() const
//{
//    return viscosityModelPtr_->strainRate();
//}



Foam::tmp<Foam::scalarField>
Foam::singlePhaseTransportModelGranular::nu(const label patchi) const
{
    return viscosityModelPtr_->nu(patchi);
}


void Foam::singlePhaseTransportModelGranular::correct()
{
    viscosityModelPtr_->correct();
}


bool Foam::singlePhaseTransportModelGranular::read()
{
    if (regIOobject::read())
    {
        return viscosityModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
