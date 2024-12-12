/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "palabos3D.h"
#include "palabos3D.hh"

#include "ibm3D.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#if __cplusplus >= 202002L
#include <numbers> // std::numbers::pi_v<T>
#endif

#include "input.h"
#include "lammps.h"
#include "lammpsWrapper.h"
#include "library.h"
#include "mpi.h"

#include "latticeDecomposition.h"
// #include "nearestTwoNeighborLattices3D.h"

using namespace plb;
using namespace std;

namespace {

using T = double;

template<typename T>
using DESCRIPTOR = descriptors::ForcedD3Q19Descriptor<T>;
//               = descriptors::ForcedN2D3Q19Descriptor;  <- does not exist in plb

template<typename T>
using DYNAMICS = GuoExternalForceBGKdynamics<T, DESCRIPTOR>;
// Or:         = BGKdynamics<T, DESCRIPTOR>;

constexpr unsigned int NMAX = 150U;

#if __cplusplus >= 202002L
constexpr T pi = std::numbers::pi_v<T>;
#else
constexpr T pi = T(M_PI);
#endif

[[nodiscard]]
T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T Nx = parameters.getNx() - 1;
    const T Ny = parameters.getNy() - 1;

    const T Nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum += ((T) 1 / (std::pow(twoNplusOne, (T) 3) *
                         std::cosh(twoNplusOne * pi * Ny / ((T) 2 * Nx))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum -= ((T) 1 / (std::pow(twoNplusOne, (T) 3) *
                         std::cosh(twoNplusOne * pi * Ny / ((T) 2 * Nx))));
    }

    T alpha = -(T) 8 * uMax * pi * pi * pi /
              (Nx * Nx * (pi * pi * pi - (T) 32 * sum)); // alpha = -dp/dz / mu
    T deltaP = -(alpha * Nu);
    return deltaP;
}

[[maybe_unused, nodiscard]]
T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const &parameters,
                     plint maxN)
{
    const T Nx = parameters.getNx() - 1;
    const T Ny = parameters.getNy() - 1;

    const T x = (T) iX - Nx / (T) 2;
    const T y = (T) iY - Ny / (T) 2;

    const T alpha = -poiseuillePressure(parameters, maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum +=
            (std::cos(twoNplusOne * pi * x / Nx) * std::cosh(twoNplusOne * pi * y / Nx) /
             (std::pow(twoNplusOne, (T) 3) *
              std::cosh(twoNplusOne * pi * Ny / ((T) 2 * Nx))));
    }

    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum -=
            (std::cos(twoNplusOne * pi * x / Nx) * std::cosh(twoNplusOne * pi * y / Nx) /
             (std::pow(twoNplusOne, (T) 3) *
              std::cosh(twoNplusOne * pi * Ny / ((T) 2 * Nx))));
    }

    sum *= ((T) 4 * alpha * Nx * Nx / std::pow(pi, (T) 3));
    sum += (alpha / (T) 2 * (x * x - Nx * Nx / (T) 4));

    return sum;
}

template<typename T>
class SquarePoiseuilleDensityAndVelocity
{
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const &parameters_,
                                       plint maxN_)
        : parameters(parameters_), maxN(maxN_)
    {}

    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T, 3> &u) const
    {
        rho = (T) 1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }

private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template<typename T>
class SquarePoiseuilleVelocity
{
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const &parameters_, plint maxN_)
        : parameters(parameters_), maxN(maxN_)
    {}

    void operator()(plint iX, plint iY, plint iZ, Array<T, 3> &u) const
    {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }

private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template<typename T>
class ShearTopVelocity
{
public:
    ShearTopVelocity(IncomprFlowParam<T> const &parameters_, plint maxN_)
        : parameters(parameters_), maxN(maxN_)
    {}

    void operator()(plint iX, plint iY, plint iZ, Array<T, 3> &u) const
    {
        u[0] = T();
        u[1] = T();
        u[2] = parameters.getLatticeU();
    }

private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template<typename T>
class ShearBottomVelocity
{
public:
    ShearBottomVelocity(IncomprFlowParam<T> const &parameters_, plint maxN_)
        : parameters(parameters_), maxN(maxN_)
    {}

    void operator()(plint iX, plint iY, plint iZ, Array<T, 3> &u) const
    {
        u[0] = T();
        u[1] = T();
        u[2] = T();
    }

private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

void squarePoiseuilleSetup(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
                           IncomprFlowParam<T> const &parameters,
                           OnLatticeBoundaryCondition3D<T, DESCRIPTOR> &boundaryCondition)
{
    const plint Nx = parameters.getNx();
    const plint Ny = parameters.getNy();
    const plint Nz = parameters.getNz();
    Box3D top = Box3D(0, Nx - 1, Ny - 1, Ny - 1, 0, Nz - 1);
    Box3D bottom = Box3D(0, Nx - 1, 0, 0, 0, Nz - 1);

    Box3D inlet = Box3D(0, Nx - 1, 1, Ny - 2, 0, 0);
    Box3D outlet = Box3D(0, Nx - 1, 1, Ny - 2, Nz - 1, Nz - 1);

    Box3D left = Box3D(0, 0, 1, Ny - 2, 1, Nz - 2);
    Box3D right = Box3D(Nx - 1, Nx - 1, 1, Ny - 2, 1, Nz - 2);
    // shear flow top bottom surface
    /*
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet,
    boundary::outflow ); boundaryCondition.setVelocityConditionOnBlockBoundaries (
    lattice, outlet, boundary::outflow );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left,
    boundary::outflow ); boundaryCondition.setVelocityConditionOnBlockBoundaries (
    lattice, right, boundary::outflow );

    setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet,
    boundary::outflow ); boundaryCondition.setVelocityConditionOnBlockBoundaries (
    lattice, outlet, boundary::outflow );
    */
    // channel flow
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, outlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, top);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, bottom);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, left);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, right);

    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));

    setBoundaryVelocity(lattice, top, Array<T, 3>((T) 0.0, (T) 0.0, (T) 0.0));
    setBoundaryVelocity(lattice, bottom, Array<T, 3>((T) 0.0, (T) 0.0, (T) 0.0));
    setBoundaryVelocity(lattice, left, Array<T, 3>((T) 0.0, (T) 0.0, (T) 0.0));
    setBoundaryVelocity(lattice, right, Array<T, 3>((T) 0.0, (T) 0.0, (T) 0.0));

    // initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
    // SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T) 1.0,
                            Array<T, 3>(0.0, 0.0, 0.0));

    lattice.initialize();
}

[[maybe_unused, nodiscard]]
T computeRMSerror(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
                  IncomprFlowParam<T> const &parameters)
{
    MultiTensorField3D<T, 3> analyticalVelocity(lattice);
    setToFunction(analyticalVelocity, analyticalVelocity.getBoundingBox(),
                  SquarePoiseuilleVelocity<T>(parameters, NMAX));
    MultiTensorField3D<T, 3> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

    // Divide by lattice velocity to normalize the error
    return 1. / parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
           std::sqrt(computeAverage(
               *computeNormSqr(*subtract(analyticalVelocity, numericalVelocity))));
}

void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
              IncomprFlowParam<T> const &parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx / dt);
    vtkOut.writeData<3, float>(*computeVelocity(lattice), "velocity", dx / dt);
    vtkOut.writeData<3, float>(*computeVorticity(*computeVelocity(lattice)), "vorticity",
                               1. / dt);
}

} // namespace

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    /*
        if (argc != 2) {
            pcout << "Error the parameters are wrong. The structure must be :\n";
            pcout << "1 : N\n";
            exit(1);
        }
    */

    // Lattic resolution
    const plint N = 1;
    // or         = atoi(argv[1]);

    // Reynolds number
    const T Re = 5e-3;

    [[maybe_unused]] const plint Nref = 50;
    [[maybe_unused]] const T uMaxRef = 0.01;
    // Characteristic velocity in lattice units (proportional to Mach number)
    const T uMax = 0.00075;
    // or        = uMaxRef / (T)N * (T) Nref; <- Needed to avoid compressibility errors.

    IncomprFlowParam<T> parameters(uMax, Re, N,
                                   20., // x-length (in dimensionless units)
                                   20., // y-length ditto
                                   80.  // z-length ditto
    );

    const T maxT = 100; // 6.6e4; // (T)0.01;

    plint iSave = 10; // 2000; // 10;
    plint iCheck = 10 * iSave;

    writeLogFile(parameters, "3D square Poiseuille");

    LammpsWrapper wrapper(argv, global::mpi().getGlobalCommunicator());
    char *inlmp = argv[1];
    wrapper.execFile(inlmp);

    // MultiTensorField3D<T,3>
    // vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    pcout << "Nx,Ny,Nz " << parameters.getNx() << " " << parameters.getNy() << " "
          << parameters.getNz() << endl;
    LatticeDecomposition lDec(parameters.getNx(), parameters.getNy(), parameters.getNz(),
                              wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution *threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 2;

    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        MultiBlockManagement3D(blockStructure, threadAttribution, envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
        new DYNAMICS<T>(parameters.getOmega()));

    // Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout << "dx " << parameters.getDeltaX() << " dt  " << parameters.getDeltaT()
          << " tau " << parameters.getTau() << endl;
    // pcout<<"51 works"<<endl;

    /*
        MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new DYNAMICS<T>(parameters.getOmega()) );*/

    // Use periodic boundary conditions.
    // lattice.periodicity().toggle(2,true);

    auto boundaryCondition = std::unique_ptr<OnLatticeBoundaryCondition3D<T, DESCRIPTOR>>(
        createLocalBoundaryCondition3D<T, DESCRIPTOR>());

    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(), parameters.getResolution(),
                                  1.0e-3);
    // coupling between lammps and palabos

    for (plint iT = 0; iT < 4e3; iT++)
    {
        lattice.collideAndStream();
    }

    T timeduration = T();
    global::timer("mainloop").start();

    for (plint iT = 0; iT < maxT; ++iT)
    {
        // for (plint iT=0; iT<2; ++iT) {
        if (iT % iSave == 0 && iT > 0)
        {
            pcout << "Saving VTK file..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT % iCheck == 0 && iT > 0)
        {
            pcout << "Timestep " << iT << " Saving checkPoint file..." << endl;
            saveBinaryBlock(lattice, "checkpoint.dat");
        }

        // lammps to calculate force
        wrapper.execCommand("run 1 pre no post no");
        // Clear and spread fluid force
        Array<T, 3> force(0, 0., 0);
        setExternalVector(lattice, lattice.getBoundingBox(),
                          DESCRIPTOR<T>::ExternalField::forceBeginsAt, force);

        ////-----classical ibm coupling-------------//
        // spreadForce3D(lattice,wrapper);
        ////// Lattice Boltzmann iteration step.
        // lattice.collideAndStream();
        ////// Interpolate and update solid position
        // interpolateVelocity3D(lattice,wrapper);
        //-----force FSI ibm coupling-------------//

        forceCoupling3D(lattice, wrapper);
        lattice.collideAndStream();
    }

    timeduration = global::timer("mainloop").stop();
    pcout << "total execution time " << timeduration << endl;
}
