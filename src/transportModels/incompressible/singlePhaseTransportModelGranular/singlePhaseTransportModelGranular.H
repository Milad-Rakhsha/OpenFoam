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

Class
    Foam::singlePhaseTransportModelGranular

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    singlePhaseTransportModelGranular.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseTransportModelGranular_H
#define singlePhaseTransportModelGranular_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class viscosityModelGranular;

/*---------------------------------------------------------------------------*\
                Class singlePhaseTransportModelGranular Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseTransportModelGranular
:
    public IOdictionary,
    public transportModel
{
    // Private Data

        autoPtr<viscosityModelGranular> viscosityModelPtr_;


    // Private Member Functions

        //- Disallow copy construct
        singlePhaseTransportModelGranular(const singlePhaseTransportModelGranular&);

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseTransportModelGranular&);


public:

    //- Runtime type information
    TypeName("singlePhaseTransportModelGranular");


    // Constructors

        //- Construct from components
        singlePhaseTransportModelGranular
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
			const volScalarField& p
        );


    //- Destructor
    virtual ~singlePhaseTransportModelGranular();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

//        virtual tmp<volScalarField> strainRate() const;


        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
