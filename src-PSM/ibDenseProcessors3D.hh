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
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

#include "ibDef.h"
#include "ibDenseCompositeDynamics3D.h"
#include "ibDenseDynamicsParticleData.h"
#include "utilsDense.h"

#include "lammps.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

namespace plb{

  template<typename T, template<typename U> class Descriptor>
  class ClearParticleData : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    ClearParticleData() {}
  
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
          for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          
            IBdynamicsParticleData<T,Descriptor>* particleData =
              getParticleDataFromCell<T,Descriptor>(cell);
          
            if(particleData) {
            // 设置颗粒数据为0
              particleData->clear();
            }
          }
        }
      }
      global::mpi().barrier();
    }
  
    virtual ClearParticleData<T,Descriptor>* clone() const {
      return new ClearParticleData<T,Descriptor>(*this);
    }
  
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
      modified[0] = modif::staticVariables;  // 确保数据在处理器边界上更新
    }
  };

  /*
   * implementation of SetSingleSphere3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    // modifies stuff, but we don't want updates running
    // after each particle set... should be modif::dataStructure
    // modified[0] = modif::dataStructure;
    modified[0] = modif::nothing; 
  }

  template<typename T, template<typename U> class Descriptor>
  SetSingleSphere3D<T,Descriptor>* SetSingleSphere3D<T,Descriptor>::clone() const 
  {
    return new SetSingleSphere3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          
          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;

          // this one actually helps, believe it or not
          __builtin_prefetch(particleData,1);
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

                    
          T const dx = xGlobal - x[0];
          T const dy = yGlobal - x[1];
          T const dz = zGlobal - x[2];

          T const dx_com = xGlobal - com[0];
          T const dy_com = yGlobal - com[1];
          T const dz_com = zGlobal - com[2];

          T const sf = calcSolidFraction(dx,dy,dz,r);
          if (sf <= SOLFRAC_MIN) continue;

          // Check the existance of this particle
          bool particleExists = false;
          for(plint i=0;i<particleData->particles.size();++i){
            if(particleData->particles[i].partId==id){
              // Update particle info
              setValues(particleData->particles[i],sf,dx_com,dy_com,dz_com);
              particleExists = true;
              break;
            }
          }

          // If the particle ID does not exist, add new particle data
          if(!particleExists){
            typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo newParticle;
            newParticle.partId = id;
            setValues(newParticle,sf,dx_com,dy_com,dz_com);
            particleData->particles.push_back(newParticle);
          }

          // if desired, initialize interior of sphere with sphere velocity
          if(initVelFlag && sf > SOLFRAC_MAX){
                        //For the case of multiple particles, we use the velocity of the particle with the highest solid fraction
            Array<T,Descriptor<T>::d> maxSfVelocity;
            T maxSf = 0;
            for(plint i=0;i<particleData->particles.size();++i){
              if (particleData->particles[i].solidFraction>maxSf){
                maxSf = particleData->particles[i].solidFraction;
                maxSfVelocity = particleData->particles[i].uPart;
              }
            }
            if(maxSf>0){
              cell.defineVelocity(maxSfVelocity);
            }
          }

        }
      }
    }
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_)
  {
    // return calcSolidFractionRec(dx_,dy_,dz_,r_);

    static plint const slicesPerDim = 5;
    static T const sliceWidth = 1./((T)slicesPerDim);
    static T const fraction = 1./((T)(slicesPerDim*slicesPerDim*slicesPerDim));
    
    // should be sqrt(3.)/2.
    // add a little to avoid roundoff errors
    static const T sqrt3half = (T) sqrt(3.1)/2.; 

    T const dist = dx_*dx_ + dy_*dy_ + dz_*dz_;

    T const r_p = r_ + sqrt3half;
    if (dist > r_p*r_p) return 0;

    T const r_m = r_ - sqrt3half;
    if (dist < r_m*r_m) return 1;

    T const r_sq = r_*r_;
    T dx_sq[slicesPerDim],dy_sq[slicesPerDim],dz_sq[slicesPerDim];

    // pre-calculate d[xyz]_sq for efficiency
    for(plint i=0;i<slicesPerDim;i++){
      T const delta = -0.5 + ((T)i+0.5)*sliceWidth;
      T const dx = dx_+delta; dx_sq[i] = dx*dx;
      T const dy = dy_+delta; dy_sq[i] = dy*dy;
      T const dz = dz_+delta; dz_sq[i] = dz*dz;
    }

    pluint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      for(plint j=0;j<slicesPerDim;j++){
        for(plint k=0;k<slicesPerDim;k++){
          n += (dx_sq[i] + dy_sq[j] + dz_sq[k] < r_sq);
        }
      }
    }

    return fraction*((T)n);
  }



  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setValues(typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo &p,
                                                  T const sf, T const dx, T const dy, T const dz)
  {    
    p.uPart.from_cArray(v);
    if(omega != 0){
      p.uPart[0] += omega[1]*dz - omega[2]*dy;
      p.uPart[1] += -omega[0]*dz + omega[2]*dx; 
      p.uPart[2] += omega[0]*dy - omega[1]*dx;
    }
    p.solidFraction = sf;
    p.partId = id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setToZero(typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo &p)
  {
    p.uPart.resetToZero();
    p.solidFraction = 0;
    p.partId = 0;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x_,
                                                   T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_)
    : x(x_), force(force_), torque(torque_), wrapper(wrapper_)
  {}

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig)
    : BoxProcessingFunctional3D_L<T,Descriptor>(orig), x(orig.x),
      force(orig.force), torque(orig.torque), wrapper(orig.wrapper) {}
  
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();
    
    // "real" domain size is nx-2 etc
    plint nx = lattice.getNx()-2, ny = lattice.getNy()-2, nz = lattice.getNz()-2;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;
          if(particleData->particles.empty()) continue;
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

          // Process each particle in the cell
          for(plint i=0;i<particleData->particles.size();++i){
            // LIGGGHTS indices start at 1
            plint const id = particleData->particles[i].partId;
            if(id < 1) continue; // no particle here

            plint const ind = wrapper.lmp->atom->map(id);
            if(ind < 0 || ind >= static_cast<plint>(x.size())) continue;

            T dx = xGlobal - x[ind][0];
            T dy = yGlobal - x[ind][1];
            T dz = zGlobal - x[ind][2];

            // minimum image convention, needed if
            // (1) PBC are used and
            // (2) both ends of PBC lie on the same processor
            if(dx > nx/2) dx -= nx;
            else if(dx < -nx/2) dx += nx;
            if(dy > ny/2) dy -= ny;
            else if(dy < -ny/2) dy += ny;
            if(dz > nz/2) dz -= nz;
            else if(dz < -nz/2) dz += nz;

            T const forceX = particleData->particles[i].hydrodynamicForce[0];
            T const forceY = particleData->particles[i].hydrodynamicForce[1];
            T const forceZ = particleData->particles[i].hydrodynamicForce[2];

            T const torqueX = dy*forceZ - dz*forceY;
            T const torqueY = -dx*forceZ + dz*forceX;
            T const torqueZ = dx*forceY - dy*forceX;

            addForce(ind,0,forceX);
            addForce(ind,1,forceY);
            addForce(ind,2,forceZ);

            addTorque(ind,0,torqueX);
            addTorque(ind,1,torqueY);
            addTorque(ind,2,torqueZ);
          }
        }
      }
    }
  }
        

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint const partId, plint const coord, T const value)
  {
    force[3*partId+coord] += value;
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint const partId, plint const coord, T const value)
  {
    torque[3*partId+coord] += value;
  }

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>* SumForceTorque3D<T,Descriptor>::clone() const
  { 
    return new SumForceTorque3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
  }

  
};
