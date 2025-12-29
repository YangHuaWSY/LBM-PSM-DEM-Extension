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


#ifndef UTILS_H_LBDEM
#define UTILS_H_LBDEM

#include "ibDenseDynamicsParticleData.h"

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell);

  
}; /* namespace plb */

// --------------------------
// implementaton starts here

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell)
  {
    // 支持所有IB dynamics类型
    static IBcompositeDynamics<T,Descriptor> const cdyn(new NoDynamics<T,Descriptor>(),false);
    static IBSmagorinskycompositeDynamics<T,Descriptor> const smagdyn(new NoDynamics<T,Descriptor>(),false);
    static IBMRTcompositeDynamics<T,Descriptor> const mrtdyn(new NoDynamics<T,Descriptor>(),false);
    static IBSmagorinskyMRTcompositeDynamics<T,Descriptor> const smagmrtdyn(new NoDynamics<T,Descriptor>(),false);
    
    static plint const cdynId = cdyn.getId();
    static plint const smagdynId = smagdyn.getId();
    static plint const mrtdynId = mrtdyn.getId();
    static plint const smagmrtdynId = smagmrtdyn.getId();
    
    Dynamics<T,Descriptor> *dyn = &( cell.getDynamics() );
    
    // 检查所有IB类型
    if(dyn->getId() == cdynId){
      return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }
    if(dyn->getId() == smagdynId){
      return &( (static_cast< IBSmagorinskycompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }
    if(dyn->getId() == mrtdynId){
      return &( (static_cast< IBMRTcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }
    if(dyn->getId() == smagmrtdynId){
      return &( (static_cast< IBSmagorinskyMRTcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }
    
    // 递归检查composite dynamics
    while(dyn->isComposite()){
      dyn = &(static_cast<CompositeDynamics<T,Descriptor>* >(dyn))->getBaseDynamics();
      if(dyn->getId() == cdynId)
          return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
      if(dyn->getId() == smagdynId)
          return &( (static_cast< IBSmagorinskycompositeDynamics<T,Descriptor>* >(dyn))->particleData );
      if(dyn->getId() == mrtdynId)
          return &( (static_cast< IBMRTcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
      if(dyn->getId() == smagmrtdynId)
          return &( (static_cast< IBSmagorinskyMRTcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }
    return 0;
}

}; /* namespace plb */

#endif /* UTILS_H_LBDEM */
