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

Namespace
    Foam::viscosityModelGranulars

Description
    A namespace for various incompressible viscosityModelGranular implementations.

Class
    Foam::viscosityModelGranular

Description
    An abstract base class for incompressible viscosityModelGranulars.

    The strain rate is defined by:

        mag(symm(grad(U)))


SourceFiles
    viscosityModelGranular.C
    viscosityModelGranularNew.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityModelGranular_H
#define viscosityModelGranular_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class viscosityModelGranular Declaration
\*---------------------------------------------------------------------------*/

class viscosityModelGranular
{

protected:

    // Protected data

        word name_;
        dictionary viscosityProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        const volScalarField& p_;


    // Private Member Functions

        //- Disallow copy construct
        viscosityModelGranular(const viscosityModelGranular&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityModelGranular&);


public:

    //- Runtime type information
    TypeName("viscosityModelGranular");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscosityModelGranular,
            dictionary,
            (
                const word& name,
                const dictionary& viscosityProperties,
                const volVectorField& U,
                const surfaceScalarField& phi,
			    const volScalarField& p


            ),
            (name, viscosityProperties, U, phi, p)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<viscosityModelGranular> New
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
		    const volScalarField& p


        );


    // Constructors

        //- Construct from components
        viscosityModelGranular
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
		    const volScalarField& p


        );


    //- Destructor
    virtual ~viscosityModelGranular()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }

        //- Return the strain rate
        tmp<volScalarField> strainRate() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;


        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
