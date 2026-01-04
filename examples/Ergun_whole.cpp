/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2025 Wuhan University China
 *
 * Author: Siyuan Wang (wsyseuswjtu@163.com)
 */

/*
This case was used as a benchmark. It consists of rectangular duct,
where the middle part is filled with particles. The aim of this example is
to determain the correctness of dense particle model, and which resolution has
the best accuracy compared with Ergun Equation.

Ergun Equation:
  Delta_P/L = A*niu_f/d^2*(1-n)^2/n^3*u + B/d*rho_f*(1-n)/n^3*u^2
  A = 150, B = 1.72 (Re = 1~2400)
  u: the velocity of duct inlet; niu_f: the fluid dynamic viscosity; rho_f: the fluid density
  n: average porosity; d: particle diameter
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include "plb_ibDense.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "library_cfd_coupling.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::D3Q19Descriptor
#define BASEDYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS IBcompositeDynamics<T, DESCRIPTOR>( new BASEDYNAMICS )

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters) {
    T x = (T)iX / parameters.getResolution() - 7.5;
    T y = (T)iY / parameters.getResolution();
    T f2 = 0.0;
    T f1 = 0.0;
    T h = 15;
    T w = 15;
    for (plint n_se=1; n_se<=100; n_se=n_se+2) {
        f2 = f2 + 1/std::pow(n_se,3)*(1-std::cosh(n_se*3.1415926*x/h)/std::cosh(n_se*3.1415926*w/2/h))*std::sin(n_se*3.1415926*y/h);
    }
    for (plint n_se=1; n_se<=100; n_se=n_se+2) {
        f1 = f1 + 1/std::pow(n_se,5)*192/std::pow(3.1415926,5)*h/w*std::tanh(n_se*3.1415926*w/2/h);
    }
    f1 = 1-f1;
    return 48 / std::pow(3.1415926,3) * parameters.getLatticeU() * f2 / f1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {
        u[0] = T();
        u[1] = T();
        u[2] = -poiseuilleVelocity(iX, iY, parameters);
    }
private:
    IncomprFlowParam<T> parameters;
};


void ErgunSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D inlet = Box3D(0, nx-1, 0, ny-1, nz-1, nz-1);
    Box3D outlet = Box3D(0, nx-1, 0, ny-1, 0, 0);

    Box3D side1(0, 0, 0, ny-1, 0, nz-1);
    Box3D side2(nx-1, nx-1, 0, ny-1, 0, nz-1);
    Box3D side3(0, nx-1, 0, 0, 0, nz-1);
    Box3D side4(0, nx-1, ny-1, ny-1, 0, nz-1);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
    setBoundaryVelocity(lattice, inlet, PoiseuilleVelocity<T>(parameters) );

    boundaryCondition.setPressureConditionOnBlockBoundaries(lattice, outlet);
    setBoundaryDensity(lattice, outlet, (T) 0.99);


}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{

  MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());

  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;

  std::string fname(createFileName("vtk", iter, 8));

  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));

  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact );

  IBscalarQuantity sf = SolidFraction;
  applyProcessingFunctional(new GetScalarQuantityFromDynamicsFunctional<T,DESCRIPTOR,T>(sf),
                           lattice.getBoundingBox(),lattice,p);

  vtkOut.writeData<float>(p,"solidfraction",1. );


  pcout << "wrote " << fname << std::endl;
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    plint N(0), demSubsteps(10);
    T vInlet(0.), maxT(0.), u_lattice(0.), interval(1000), l_p(0.1), insertT(0.1);
    std::string outDir;

    try {
        global::argv(1).read(N);
        global::argv(2).read(vInlet); // physical unit
        global::argv(3).read(interval);
        global::argv(4).read(maxT);
        global::argv(5).read(insertT);
        global::argv(6).read(u_lattice); // lattice unit
        global::argv(7).read(l_p);
        global::argv(8).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : grid points along particle diameter\n";
        pcout << "2 : inlet velocity\n";
        pcout << "3 : interval\n";
        pcout << "4 : maximal run time\n";
        pcout << "5 : particle insert time\n";
        pcout << "6 : expected maximum lattice velocity\n";
        pcout << "7 : characteristic length\n";
        pcout << "8 : outDir\n";
        exit(1);
    }

// Parameters from fluid
    const T g = 9.81;
    const T mu_f = 0.01; // dyne/cm^2
    const T rho_f = 1.00;
    const T nu_f = mu_f/rho_f; // cm^2/s
    const T lx = 1.5, ly = 1.5, lz = 4.;

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);
    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

    wrapper.setVariable("rho_fluid",rho_f);
    wrapper.setVariable("dmp_dir",demOutDir);
    wrapper.setVariable("lx",lx);
    wrapper.setVariable("ly",ly);
    wrapper.setVariable("lz",lz);
    wrapper.execFile("in.lbdem");

    PhysUnits3D<T> units(l_p,vInlet,nu_f,lx,ly,lz,N,u_lattice,rho_f);
    units.setLbOffset(0.,0.,0.);

    IncomprFlowParam<T> parameters(units.getLbParam());

    const T vtkT = 1.0/interval;
    const T logT = 0.01/interval;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint insertSteps = units.getLbSteps(insertT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D Ergun test");

    pcout << "-----------------------------------\n";
    pcout << "grid size: "
          << parameters.getNx() << " "
          << parameters.getNy() << " "
          << parameters.getNz() << std::endl;
    pcout << "-----------------------------------" << std::endl;

    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);

    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();

    plint envelopeWidth = 1;

    MultiBlockLattice3D<T, DESCRIPTOR>
      lattice (MultiBlockManagement3D (blockStructure,
                                       threadAttribution,
                                       envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    defineDynamics(lattice,lattice.getBoundingBox(),new DYNAMICS);
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createInterpBoundaryCondition3D<T,DESCRIPTOR>();

    lattice.periodicity().toggleAll(false);
    lattice.initialize();

    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n"
          << "dt_phys: " << dt_phys << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "vtkT: " << vtkT << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "logT: " << logT << "\n"
          << "logSteps: " << logSteps << "\n"
          << "grid size: "
          << parameters.getNx() << " "
          << parameters.getNy() << " "
          << parameters.getNz() << std::endl;

    clock_t start = clock();

    T dt_dem = dt_phys/(T)demSubsteps;
    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
    wrapper.execFile("in2.lbdem");
    wrapper.runUpto(demSubsteps-1);

    util::ValueTracer<T> converge(u_lattice,lx,5.0e-3);
    plint iT=0;
// Particle Insertion
    for (; iT<insertSteps; ++iT) {

      wrapper.run(demSubsteps);

      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = ((T)difftime(end,start))/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*((T)logSteps)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();

        pcout << "Physical time is " << units.getPhysTime(iT) << endl;

      }
    }
    wrapper.execCommand("unfix ins");

// Particle Deposition & Stabiliztion
    for (; iT<maxSteps; ++iT) {

      setSpheresOnLattice(lattice,wrapper,units,false);

      if(iT%vtkSteps == 0){
        writeVTK(lattice,parameters,units,iT);
      }

      lattice.collideAndStream();
      getForcesFromLattice(lattice,wrapper,units);
      wrapper.run(demSubsteps);

      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = ((T)difftime(end,start))/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*((T)logSteps)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();

        pcout << "Physical time is " << units.getPhysTime(iT) << endl;

      }
      converge.takeValue(getStoredAverageEnergy(lattice),true);
      if (converge.hasConverged()){
          break;
      }
    }
    writeVTK(lattice,parameters,units,2);

// Ergun Test
    ErgunSetup(lattice, parameters, *boundaryCondition);
    for (; iT<maxSteps; ++iT) {

      setSpheresOnLattice(lattice,wrapper,units,false);

      if(iT%vtkSteps == 0){
        writeVTK(lattice,parameters,units,iT);
      }

      lattice.collideAndStream();
      getForcesFromLattice(lattice,wrapper,units);
      wrapper.run(demSubsteps);

      if(iT%logSteps == 0){
        clock_t end = clock();
        T time = ((T)difftime(end,start))/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*((T)logSteps)))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s" << std::endl;
        start = clock();

        pcout << "Physical time is " << units.getPhysTime(iT) << endl;

      }

    }
    writeVTK(lattice,parameters,units,3);

    delete boundaryCondition;
}
