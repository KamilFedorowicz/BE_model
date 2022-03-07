#include "BerisEdwards.H"
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(BerisEdwards, 0);
    addToRunTimeSelectionTable(constitutiveEq2, BerisEdwards, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::BerisEdwards::BerisEdwards
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq2(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Q_
    (
        IOobject
        (
            "Q" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    sigma2_
    (
        IOobject
        (
            "sigma2" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),


    tauTotal_
    (
        IOobject
        (
            "tauTotal" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),


    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    a_(dict.lookup("a")),
    b_(dict.lookup("b")),
    c_(dict.lookup("c")),
    L_(dict.lookup("L")),
    xi_(dict.lookup("xi")),
    gamma_(dict.lookup("gamma"))

{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::BerisEdwards::correct()
{
  // Velocity gradient tensor
    volTensorField gradU = fvc::grad(U());

    // Symmetric velocity gradient tensor
    volTensorField D = ( gradU+T(gradU) )/2;

    // Antisymmetric velocity gradient tensor
    volTensorField omega = skew(gradU);

tensor I(1,0,0,0,1,0,0,0,1);

	volTensorField S=( (xi_*D - omega) & (Q_ + I/3 ) ) 
			+ ( (Q_ + I/3 ) & (xi_*D + omega) )
			- ( 2*xi_ *  (Q_ + I/3 )*tr(Q_ & gradU) );


	volTensorField H=fvc::laplacian(L_, Q_) - a_*Q_ 
				+ b_*( (Q_ & Q_) - (tr(Q_ & Q_ ) )* I/3 ) - c_ * Q_ * ( tr(Q_ & Q_) );



//Q-tensor transport equation
    fvTensorMatrix QEqn
    (
	fvm::ddt(Q_) 
	+ fvm::div(phi(), Q_)
	==

	S + gamma_ * H

    );

	QEqn.relax();
	QEqn.solve();


	Q_=( Q_ + T(Q_) )/2;

forAll(Q_, cellI)
{
Q_[cellI][8]=-Q_[cellI][0]-Q_[cellI][4];
}


// calculate gradients of Q


	volScalarField Q11=Q_.component(0);
	volScalarField Q12=Q_.component(1);
	volScalarField Q13=Q_.component(2);
	volScalarField Q21=Q_.component(3);
	volScalarField Q22=Q_.component(4);
	volScalarField Q23=Q_.component(5);
	volScalarField Q31=Q_.component(6);
	volScalarField Q32=Q_.component(7);
	volScalarField Q33=Q_.component(8);



	volVectorField gradQ11=fvc::grad(Q11);
	volVectorField gradQ12=fvc::grad(Q12);
	volVectorField gradQ13=fvc::grad(Q13);
	volVectorField gradQ21=fvc::grad(Q21);
	volVectorField gradQ22=fvc::grad(Q22);
	volVectorField gradQ23=fvc::grad(Q23);
	volVectorField gradQ31=fvc::grad(Q31);
	volVectorField gradQ32=fvc::grad(Q32);
	volVectorField gradQ33=fvc::grad(Q33);



tensor Ixx(1,0,0,0,0,0,0,0,0);
tensor Ixy(0,1,0,0,0,0,0,0,0);
tensor Ixz(0,0,1,0,0,0,0,0,0);
tensor Iyx(0,0,0,1,0,0,0,0,0);
tensor Iyy(0,0,0,0,1,0,0,0,0);
tensor Iyz(0,0,0,0,0,1,0,0,0);
tensor Izx(0,0,0,0,0,0,1,0,0);
tensor Izy(0,0,0,0,0,0,0,1,0);
tensor Izz(0,0,0,0,0,0,0,0,1);

//elastic stress tensor

	volScalarField sigmaExx=
			gradQ11.component(0)*gradQ11.component(0)
			+gradQ12.component(0)*gradQ12.component(0)
			+gradQ13.component(0)*gradQ13.component(0)
			+gradQ21.component(0)*gradQ21.component(0)
			+gradQ22.component(0)*gradQ22.component(0)
			+gradQ23.component(0)*gradQ23.component(0)
			+gradQ31.component(0)*gradQ31.component(0)
			+gradQ32.component(0)*gradQ32.component(0)
			+gradQ33.component(0)*gradQ33.component(0)	;

volScalarField	sigmaExy=
			gradQ11.component(0)*gradQ11.component(1)
			+gradQ12.component(0)*gradQ12.component(1)
			+gradQ13.component(0)*gradQ13.component(1)
			+gradQ21.component(0)*gradQ21.component(1)
			+gradQ22.component(0)*gradQ22.component(1)
			+gradQ23.component(0)*gradQ23.component(1)
			+gradQ31.component(0)*gradQ31.component(1)
			+gradQ32.component(0)*gradQ32.component(1)
			+gradQ33.component(0)*gradQ33.component(1)	;

volScalarField	sigmaExz=
			gradQ11.component(0)*gradQ11.component(2)
			+gradQ12.component(0)*gradQ12.component(2)
			+gradQ13.component(0)*gradQ13.component(2)
			+gradQ21.component(0)*gradQ21.component(2)
			+gradQ22.component(0)*gradQ22.component(2)
			+gradQ23.component(0)*gradQ23.component(2)
			+gradQ31.component(0)*gradQ31.component(2)
			+gradQ32.component(0)*gradQ32.component(2)
			+gradQ33.component(0)*gradQ33.component(2)	;

volScalarField	sigmaEyx=
			gradQ11.component(1)*gradQ11.component(0)
			+gradQ12.component(1)*gradQ12.component(0)
			+gradQ13.component(1)*gradQ13.component(0)
			+gradQ21.component(1)*gradQ21.component(0)
			+gradQ22.component(1)*gradQ22.component(0)
			+gradQ23.component(1)*gradQ23.component(0)
			+gradQ31.component(1)*gradQ31.component(0)
			+gradQ32.component(1)*gradQ32.component(0)
			+gradQ33.component(1)*gradQ33.component(0)	;


volScalarField	sigmaEyy=
			gradQ11.component(1)*gradQ11.component(1)
			+gradQ12.component(1)*gradQ12.component(1)
			+gradQ13.component(1)*gradQ13.component(1)
			+gradQ21.component(1)*gradQ21.component(1)
			+gradQ22.component(1)*gradQ22.component(1)
			+gradQ23.component(1)*gradQ23.component(1)
			+gradQ31.component(1)*gradQ31.component(1)
			+gradQ32.component(1)*gradQ32.component(1)
			+gradQ33.component(1)*gradQ33.component(1)	;

volScalarField	sigmaEyz=
			gradQ11.component(1)*gradQ11.component(2)
			+gradQ12.component(1)*gradQ12.component(2)
			+gradQ13.component(1)*gradQ13.component(2)
			+gradQ21.component(1)*gradQ21.component(2)
			+gradQ22.component(1)*gradQ22.component(2)
			+gradQ23.component(1)*gradQ23.component(2)
			+gradQ31.component(1)*gradQ31.component(2)
			+gradQ32.component(1)*gradQ32.component(2)
			+gradQ33.component(1)*gradQ33.component(2)	;

volScalarField	sigmaEzx=
			gradQ11.component(2)*gradQ11.component(0)
			+gradQ12.component(2)*gradQ12.component(0)
			+gradQ13.component(2)*gradQ13.component(0)
			+gradQ21.component(2)*gradQ21.component(0)
			+gradQ22.component(2)*gradQ22.component(0)
			+gradQ23.component(2)*gradQ23.component(0)
			+gradQ31.component(2)*gradQ31.component(0)
			+gradQ32.component(2)*gradQ32.component(0)
			+gradQ33.component(2)*gradQ33.component(0)	;

volScalarField	sigmaEzy=
			gradQ11.component(2)*gradQ11.component(1)
			+gradQ12.component(2)*gradQ12.component(1)
			+gradQ13.component(2)*gradQ13.component(1)
			+gradQ21.component(2)*gradQ21.component(1)
			+gradQ22.component(2)*gradQ22.component(1)
			+gradQ23.component(2)*gradQ23.component(1)
			+gradQ31.component(2)*gradQ31.component(1)
			+gradQ32.component(2)*gradQ32.component(1)
			+gradQ33.component(2)*gradQ33.component(1)	;


volScalarField	sigmaEzz=
			gradQ11.component(2)*gradQ11.component(2)
			+gradQ12.component(2)*gradQ12.component(2)
			+gradQ13.component(2)*gradQ13.component(2)
			+gradQ21.component(2)*gradQ21.component(2)
			+gradQ22.component(2)*gradQ22.component(2)
			+gradQ23.component(2)*gradQ23.component(2)
			+gradQ31.component(2)*gradQ31.component(2)
			+gradQ32.component(2)*gradQ32.component(2)
			+gradQ33.component(2)*gradQ33.component(2)	;


sigma2_= Ixx*sigmaExx + Ixy*sigmaExy + Ixz*sigmaExz +
	Iyx*sigmaEyx + Iyy*sigmaEyy + Iyz*sigmaEyz + 
	Izx*sigmaEzx + Izy*sigmaEzy + Izz*sigmaEzz	;

//viscous stress tensor
	volTensorField sigma_= (H & Q_) - (Q_ & H)
				- xi_ * ( (Q_ + I/3 ) & H ) - xi_ * ( H & (Q_ + I/3 ) )
				+ 2 * xi_ * ( (Q_ + I/3 ) * (Q_ && H ) ) ;

// elastic+viscous non-Newtonian stresses

	tau_ = ( -L_* sigma2_ + sigma_ );


// total stress
tauTotal_=tau_+2*etaS_*D;

  tau_.correctBoundaryConditions();

}


// ************************************************************************* //
