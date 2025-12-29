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

#ifndef IB_COMPOSITE_DYNAMICS_HH_LBDEM
#define IB_COMPOSITE_DYNAMICS_HH_LBDEM

#include "ibDef.h"
#include "ibDenseProcessors3D.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/mrtTemplates.h"
#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/smagorinskyDynamics.h"


namespace plb {

/* *************** Class IBcompositeDynamics *********************************************** */
  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBcompositeDynamics<T,Descriptor> >("IBcomposite");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fPre =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
							 bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(const IBcompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::~IBcompositeDynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>* IBcompositeDynamics<T,Descriptor>::clone() const {
    return new IBcompositeDynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {
    serializer.addValue<plint>(particleData.particles.size());
    for(plint iPart=0;iPart<static_cast<plint>(particleData.particles.size());++iPart) {
      const typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.hydrodynamicForce[i]);

      serializer.addValue(particle.solidFraction);
      serializer.addValue<plint>(particle.partId);
    }


    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.particles.clear();

    plint numParticles;
    unserializer.readValue(numParticles);

    for(plint iPart=0;iPart<numParticles;++iPart) {
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo particle;

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.hydrodynamicForce[i]);

      unserializer.readValue(particle.solidFraction);
      unserializer.readValue<plint>(particle.partId);

     particleData.particles.push_back(particle);
    }
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.getTotalSolidFraction() > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u); 
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }
  
  /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    prepareCollision(cell);

    for(plint iPart=0;iPart<particleData.particles.size();++iPart) {
      particleData.particles[iPart].hydrodynamicForce.resetToZero();
    }

    T totalSolidFraction = particleData.getTotalSolidFraction();

    if(totalSolidFraction <= SOLFRAC_MAX)
      CompositeDynamics<T,Descriptor>::collide(cell,statistics);
      // this->getBaseDynamics().collide(cell,statistics);
    
    if(totalSolidFraction < SOLFRAC_MIN)
      return;
    
    if(particleData.particles.empty())
      return;

    T rhoBar(0.);
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);

    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j); 
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         j, jSqr, fEq );

    T omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getOmega();

    if(particleData.particles.size()==1){ /* Single particle */
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[0];

      Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
      T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
      dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                          uPart, uPartSqr, fEqSolid );

      /* Weighting coefficient and collision operators.
      B1 = particle.solidFraction
      B2 = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo)
      C1 = (fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      C2 = (fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      C3 = (fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      The original version is B2C3.
      The convergence of C1 is terrible.
      */

      if(particle.solidFraction > SOLFRAC_MAX){

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2
        //T const bb0 = 0;//C1

        cell[0] = fPre[0] + bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = ((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = ((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = ((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = ((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = ((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = ((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] = fPre[iPop] + bb;
          cell[iOpp] = fPre[iOpp] + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      } else { /* particle.solidFraction < SOLFRAC_MAX */
        T const ooo = 1./omega - 0.5;

        T const B = particle.solidFraction*ooo/((1.- particle.solidFraction) + ooo);//B2
        //T const B = particle.solidFraction;//B1

        T const oneMinB = 1. - B;

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2
        //T const bb0 = 0;//C1

        cell[0] = fPre[0] + oneMinB*(cell[0] - fPre[0]) + B*bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = B*((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = B*((fPre[iPop] - fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = B*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = B*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = B*((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = B*((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] = fPre[iPop] + oneMinB*(cell[iPop] - fPre[iPop]) + bb;
          cell[iOpp] = fPre[iOpp] + oneMinB*(cell[iOpp] - fPre[iOpp]) + bbOpp;

         for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
             -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      }
    } else { /* Multiple particle */
      T const ooo = 1./omega - 0.5;

      for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
        cell[iPop] = fPre[iPop] + omega*(fEq[iPop] - fPre[iPop]);
      }
      for(plint iPart=0;iPart<particleData.particles.size();++iPart){
        typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];
        Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
        T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
        dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                              uPart, uPartSqr, fEqSolid );

        T const Bns = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
        //T const Bns = particle.solidFraction;//B1

        //cell[0] += -Bns*omega*(fEq[0] - fPre[0]) + Bns*(-fEq[0] + fEqSolid[0]);//C3
        cell[0] += -Bns*omega*(fEq[0] 
                  - fPre[0]) + Bns*((fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]));//C2
        //cell[0] += -Bns*omega*(fEq[0] - fPre[0]);//C1

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = Bns*((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = Bns*((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = Bns*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = Bns*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = Bns*((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = Bns*((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] += -Bns*omega*(fEq[iPop] - fPre[iPop]) + bb;
          cell[iOpp] += -Bns*omega*(fEq[iOpp] - fPre[iOpp]) + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
          }
        }
      } /* Multiple particles */

  } /* IBcompositeDynamics::collide() */

/* *************** Class IBSmagorinskycompositeDynamics *********************************************** */
  template<typename T, template<typename U> class Descriptor>
  int IBSmagorinskycompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBSmagorinskycompositeDynamics<T,Descriptor> >("IBSmagorinskycomposite");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskycompositeDynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskycompositeDynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskycompositeDynamics<T,Descriptor>::fPre =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskycompositeDynamics<T,Descriptor>::IBSmagorinskycompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
               bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_)
  { }

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskycompositeDynamics<T,Descriptor>::IBSmagorinskycompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskycompositeDynamics<T,Descriptor>::IBSmagorinskycompositeDynamics(const IBSmagorinskycompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig)
  { }

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskycompositeDynamics<T,Descriptor>::~IBSmagorinskycompositeDynamics() {}

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskycompositeDynamics<T,Descriptor>* IBSmagorinskycompositeDynamics<T,Descriptor>::clone() const {
    return new IBSmagorinskycompositeDynamics<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  int IBSmagorinskycompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskycompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {
    serializer.addValue<plint>(particleData.particles.size());
    for(plint iPart=0;iPart<static_cast<plint>(particleData.particles.size());++iPart) {
      const typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.hydrodynamicForce[i]);

      serializer.addValue(particle.solidFraction);
      serializer.addValue<plint>(particle.partId);
    }


    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskycompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.particles.clear();

    plint numParticles;
    unserializer.readValue(numParticles);

    for(plint iPart=0;iPart<numParticles;++iPart) {
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo particle;

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.hydrodynamicForce[i]);

      unserializer.readValue(particle.solidFraction);
      unserializer.readValue<plint>(particle.partId);

     particleData.particles.push_back(particle);
    }
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskycompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.getTotalSolidFraction() > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskycompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell,
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u);
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }

  /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskycompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    prepareCollision(cell);

    for(plint iPart=0;iPart<particleData.particles.size();++iPart) {
      particleData.particles[iPart].hydrodynamicForce.resetToZero();
    }

    T totalSolidFraction = particleData.getTotalSolidFraction();

    //if(totalSolidFraction <= SOLFRAC_MAX)
      //CompositeDynamics<T,Descriptor>::collide(cell,statistics);
      // this->getBaseDynamics().collide(cell,statistics);

    if(totalSolidFraction < SOLFRAC_MIN)
      return;

    if(particleData.particles.empty())
      return;

    T rhoBar(0.);
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);

    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         j, jSqr, fEq );

    T omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getDynamicParameter(1011,cell);

    if(particleData.particles.size()==1){ /* Single particle */
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[0];

      Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
      T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
      dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                          uPart, uPartSqr, fEqSolid );

      /* Weighting coefficient and collision operators.
      B1 = particle.solidFraction
      B2 = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo)
      C1 = (fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      C2 = (fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      C3 = (fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop])
      The original version is B2C3.
      The convergence of C1 is terrible.
      */

      if(particle.solidFraction > SOLFRAC_MAX){

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2
        //T const bb0 = 0;//C1

        cell[0] = fPre[0] + bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = ((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = ((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = ((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = ((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = ((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = ((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] = fPre[iPop] + bb;
          cell[iOpp] = fPre[iOpp] + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      } else { /* particle.solidFraction < SOLFRAC_MAX */
        T const ooo = 1./omega - 0.5;

        T const B = particle.solidFraction*ooo/((1.- particle.solidFraction) + ooo);//B2
        //T const B = particle.solidFraction;//B1

        T const oneMinB = 1. - B;

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2
        //T const bb0 = 0;//C1

        cell[0] = fPre[0] + oneMinB*omega*(fEq[0] - fPre[0]) + B*bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = B*((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = B*((fPre[iPop] - fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = B*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = B*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = B*((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = B*((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] = fPre[iPop] + oneMinB*omega*(fEq[iPop] - fPre[iPop]) + bb;
          cell[iOpp] = fPre[iOpp] + oneMinB*omega*(fEq[iOpp] - fPre[iOpp]) + bbOpp;

         for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
             -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      }
    } else { /* Multiple particle */
      T const ooo = 1./omega - 0.5;

      for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
        cell[iPop] = fPre[iPop] + omega*(fEq[iPop] - fPre[iPop]);
      }
      for(plint iPart=0;iPart<particleData.particles.size();++iPart){
        typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];
        Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
        T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
        dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                              uPart, uPartSqr, fEqSolid );

        T const Bns = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
        //T const Bns = particle.solidFraction;//B1

        //cell[0] += -Bns*omega*(fEq[0] - fPre[0]) + Bns*(-fEq[0] + fEqSolid[0]);//C3
        cell[0] += -Bns*omega*(fEq[0]
                  - fPre[0]) + Bns*((fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]));//C2
        //cell[0] += -Bns*omega*(fEq[0] - fPre[0]);//C1

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = Bns*((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = Bns*((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = Bns*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = Bns*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2
          //T const bb = Bns*((fPre[iOpp]-fEqSolid[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C1
          //T const bbOpp = Bns*((fPre[iPop]-fEqSolid[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C1

          cell[iPop] += -Bns*omega*(fEq[iPop] - fPre[iPop]) + bb;
          cell[iOpp] += -Bns*omega*(fEq[iOpp] - fPre[iOpp]) + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
          }
        }
      } /* Multiple particles */

  } /* IBSmagorinskycompositeDynamics::collide() */

/* *************** Class IBMRTcompositeDynamics *********************************************** */
  template<typename T, template<typename U> class Descriptor>
  int IBMRTcompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBMRTcompositeDynamics<T,Descriptor> >("IBMRTcomposite");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBMRTcompositeDynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBMRTcompositeDynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBMRTcompositeDynamics<T,Descriptor>::fPre =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  IBMRTcompositeDynamics<T,Descriptor>::IBMRTcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
               bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBMRTcompositeDynamics<T,Descriptor>::IBMRTcompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBMRTcompositeDynamics<T,Descriptor>::IBMRTcompositeDynamics(const IBMRTcompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBMRTcompositeDynamics<T,Descriptor>::~IBMRTcompositeDynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBMRTcompositeDynamics<T,Descriptor>* IBMRTcompositeDynamics<T,Descriptor>::clone() const {
    return new IBMRTcompositeDynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBMRTcompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBMRTcompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {
    serializer.addValue<plint>(particleData.particles.size());
    for(plint iPart=0;iPart<static_cast<plint>(particleData.particles.size());++iPart) {
      const typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.hydrodynamicForce[i]);

      serializer.addValue(particle.solidFraction);
      serializer.addValue<plint>(particle.partId);
    }


    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBMRTcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.particles.clear();

    plint numParticles;
    unserializer.readValue(numParticles);

    for(plint iPart=0;iPart<numParticles;++iPart) {
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo particle;

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.hydrodynamicForce[i]);

      unserializer.readValue(particle.solidFraction);
      unserializer.readValue<plint>(particle.partId);

     particleData.particles.push_back(particle);
    }
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBMRTcompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.getTotalSolidFraction() > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  template<typename T, template<typename U> class Descriptor>
  void IBMRTcompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u); 
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }
  
  /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBMRTcompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    prepareCollision(cell);

    for(plint iPart=0;iPart<particleData.particles.size();++iPart) {
      particleData.particles[iPart].hydrodynamicForce.resetToZero();
    }

    T totalSolidFraction = particleData.getTotalSolidFraction();

    if(totalSolidFraction <= SOLFRAC_MAX)
      CompositeDynamics<T,Descriptor>::collide(cell,statistics);
      // this->getBaseDynamics().collide(cell,statistics);
    
    if(totalSolidFraction < SOLFRAC_MIN)
      return;
    
    if(particleData.particles.empty())
      return;

    T rhoBar(0.);
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);

    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j); 
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         j, jSqr, fEq );

    T omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getOmega();

    if(particleData.particles.size()==1){ /* Single particle */
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[0];

      Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
      T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
      dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                          uPart, uPartSqr, fEqSolid );

      if(particle.solidFraction > SOLFRAC_MAX){

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2

        cell[0] = fPre[0] + bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = ((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = ((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = ((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = ((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2


          cell[iPop] = fPre[iPop] + bb;
          cell[iOpp] = fPre[iOpp] + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      } else { /* particle.solidFraction < SOLFRAC_MAX */

        T const ooo = 1./omega - 0.5;

        T const B = particle.solidFraction*ooo/((1.- particle.solidFraction) + ooo);//B2
        //T const B = particle.solidFraction;//B1

        T const oneMinB = 1. - B;

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2
        
        cell[0] = fPre[0]+oneMinB*(cell[0]-fPre[0])+B*bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = B*((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = B*((fPre[iPop] - fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = B*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = B*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2

          cell[iPop] = fPre[iPop]+oneMinB*(cell[iPop]-fPre[iPop])+bb;
          cell[iOpp] = fPre[iOpp]+oneMinB*(cell[iOpp]-fPre[iOpp])+bbOpp;

         for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
             -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      }
    } else { /* Multiple particle */

      T const ooo = 1./omega - 0.5;
      T const Bn = totalSolidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
      //T const Bn = totalSolidFraction;//B1
      T const oneMinBn = 1. - Bn;

      for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
        cell[iPop] = fPre[iPop]+oneMinBn*(cell[iPop]-fPre[iPop]);
      }

      for(plint iPart=0;iPart<particleData.particles.size();++iPart){
        typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];
        Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
        T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
        dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                              uPart, uPartSqr, fEqSolid );

        T const Bns = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
        //T const Bns = particle.solidFraction;//B1

        //cell[0] += Bns*(-fEq[0] + fEqSolid[0]);//C3
        cell[0] += Bns*((fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]));//C2

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = Bns*((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = Bns*((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = Bns*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = Bns*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2

          cell[iPop] += bb;
          cell[iOpp] += bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
          }
        }
      } /* Multiple particles */

  } /* IBMRTcompositeDynamics::collide() */

/* *************** Class IBSmagorinskyMRTcompositeDynamics *********************************************** */
  template<typename T, template<typename U> class Descriptor>
  int IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBSmagorinskyMRTcompositeDynamics<T,Descriptor> >("IBSmagorinskyMRTcomposite");

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::fPre =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::IBSmagorinskyMRTcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
               bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::IBSmagorinskyMRTcompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::IBSmagorinskyMRTcompositeDynamics(const IBSmagorinskyMRTcompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig)
  { }
  
  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::~IBSmagorinskyMRTcompositeDynamics() {}
  
  template<typename T, template<typename U> class Descriptor>
  IBSmagorinskyMRTcompositeDynamics<T,Descriptor>* IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::clone() const {
    return new IBSmagorinskyMRTcompositeDynamics<T,Descriptor>(*this);
  }
  
  template<typename T, template<typename U> class Descriptor>
  int IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {
    serializer.addValue<plint>(particleData.particles.size());
    for(plint iPart=0;iPart<static_cast<plint>(particleData.particles.size());++iPart) {
      const typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        serializer.addValue(particle.hydrodynamicForce[i]);

      serializer.addValue(particle.solidFraction);
      serializer.addValue<plint>(particle.partId);
    }


    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    particleData.particles.clear();

    plint numParticles;
    unserializer.readValue(numParticles);

    for(plint iPart=0;iPart<numParticles;++iPart) {
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo particle;

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.uPart[i]);

      for(plint i=0;i<Descriptor<T>::d;i++)
        unserializer.readValue(particle.hydrodynamicForce[i]);

      unserializer.readValue(particle.solidFraction);
      unserializer.readValue<plint>(particle.partId);

     particleData.particles.push_back(particle);
    }
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }
  
  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.getTotalSolidFraction() > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, 
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u); 
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }
  
  /// Implementation of the collision step
  template<typename T, template<typename U> class Descriptor>
  void IBSmagorinskyMRTcompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    prepareCollision(cell);

    for(plint iPart=0;iPart<particleData.particles.size();++iPart) {
      particleData.particles[iPart].hydrodynamicForce.resetToZero();
    }

    T totalSolidFraction = particleData.getTotalSolidFraction();

    if(totalSolidFraction <= SOLFRAC_MAX)
      CompositeDynamics<T,Descriptor>::collide(cell,statistics);
      // this->getBaseDynamics().collide(cell,statistics);
    
    if(totalSolidFraction < SOLFRAC_MIN)
      return;
    
    if(particleData.particles.empty())
      return;

    T rhoBar(0.);
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);

    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j); 
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         j, jSqr, fEq );

    T omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getDynamicParameter(1011,cell);
    //T omega = CompositeDynamics<T,Descriptor>::getBaseDynamics().getOmega();

    if(particleData.particles.size()==1){ /* Single particle */
      typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[0];

      Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
      T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
      dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                          uPart, uPartSqr, fEqSolid );

      if(particle.solidFraction > SOLFRAC_MAX){

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2

        cell[0] = fPre[0] + bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = ((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = ((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = ((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = ((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2


          cell[iPop] = fPre[iPop] + bb;
          cell[iOpp] = fPre[iOpp] + bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      } else { /* particle.solidFraction < SOLFRAC_MAX */

        T const ooo = 1./omega - 0.5;

        T const B = particle.solidFraction*ooo/((1.- particle.solidFraction) + ooo);//B2
        //T const B = particle.solidFraction;//B1

        T const oneMinB = 1. - B;

        //T const bb0 = -fEq[0] + fEqSolid[0];//C3
        T const bb0 = (fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]);//C2

        cell[0] = fPre[0] + oneMinB*(cell[0]-fPre[0]) + B*bb0;

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = B*((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = B*((fPre[iPop] - fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = B*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = B*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2

          cell[iPop] = fPre[iPop] + oneMinB*(cell[iPop]-fPre[iPop]) + bb;
          cell[iOpp] = fPre[iOpp] + oneMinB*(cell[iOpp]-fPre[iOpp]) + bbOpp;

         for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
             -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
        }
      }
    } else { /* Multiple particle */
      T const ooo = 1./omega - 0.5;
      T const Bn = totalSolidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
      //T const Bn = totalSolidFraction;//B1
      T const oneMinBn = 1. - Bn;

      for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
        cell[iPop] = fPre[iPop] + oneMinBn*(cell[iPop]-fPre[iPop]);
      }
      for(plint iPart=0;iPart<particleData.particles.size();++iPart){
        typename IBdynamicsParticleData<T,Descriptor>::ParticleInfo& particle = particleData.particles[iPart];
        Array<T,Descriptor<T>::d> uPart = particle.uPart*(1.+rhoBar);
        T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
        dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                              uPart, uPartSqr, fEqSolid );

        T const Bns = particle.solidFraction*ooo/((1.- totalSolidFraction) + ooo);//B2
        //T const Bns = particle.solidFraction;//B1

        //cell[0] += Bns*(-fEq[0] + fEqSolid[0]);//C3
        cell[0] += Bns*((fEqSolid[0] - fPre[0]) + (1-omega)*(fPre[0]-fEq[0]));//C2

        for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
          plint const iOpp = iPop + Descriptor<T>::q/2;

          //T const bb = Bns*((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));//C3
          //T const bbOpp = Bns*((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));//C3
          T const bb = Bns*((fEqSolid[iPop]-fPre[iPop]) + (1-omega)*(fPre[iPop]-fEq[iPop]));//C2
          T const bbOpp = Bns*((fEqSolid[iOpp]-fPre[iOpp]) + (1-omega)*(fPre[iOpp]-fEq[iOpp]));//C2

          cell[iPop] += bb;
          cell[iOpp] += bbOpp;

          for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
            particle.hydrodynamicForce[iDim]
              -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
          }
        }
      } /* Multiple particles */

  } /* IBSmagorinskyMRTcompositeDynamics::collide() */

}; /* namespace plb */

#endif /* IB_COMPOSITE_DYNAMICS_HH_LBDEM */
