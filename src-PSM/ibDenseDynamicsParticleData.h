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

#ifndef IB_DYNAMICS_PARTICLE_DATA_H_LBDEM
#define IB_DYNAMICS_PARTICLE_DATA_H_LBDEM

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct IBdynamicsParticleData {
public:
    IBdynamicsParticleData() {}

    // Info for each particle
    struct ParticleInfo {
        plint partId;
        T solidFraction;
        Array<T,Descriptor<T>::d> uPart;
        Array<T,Descriptor<T>::d> hydrodynamicForce;

        ParticleInfo()
            : partId(0), solidFraction(0.){
                uPart.resetToZero();
                hydrodynamicForce.resetToZero();
            }
    };

    std::vector<ParticleInfo> particles;

    T getTotalSolidFraction(){
        T totalFraction = 0;
        for(plint i=0; i<particles.size(); i++) {
            totalFraction += particles[i].solidFraction;
        }
        // Cap at 1.0 to avoid invalid solid fractions
        if (totalFraction > 1.) {
            for(plint i=0; i<particles.size(); i++) {
                particles[i].solidFraction=particles[i].solidFraction-(totalFraction-1.0)/particles.size();
            }
        }
        return (totalFraction > 1.) ? 1. : totalFraction;
    }

    // Add or update particle data
    void addOrUpdateParticle(plint id, T fraction, Array<T,Descriptor<T>::d> const& velocity) {
        for(plint i=0; i<particles.size(); ++i) {
            if(particles[i].partId == id) {
                // Update existing particle data
                particles[i].solidFraction = fraction;
                particles[i].uPart = velocity;
                return;
            }
        }
        // Add new particle if not found
        ParticleInfo newParticle;
        newParticle.partId = id;
        newParticle.solidFraction = fraction;
        newParticle.uPart = velocity;
        particles.push_back(newParticle);
    }

    void clear() {
        particles.clear();
    }
  };

}; /* namespace plb */ 

#endif /* IB_DYNAMICS_PARTICLE_DATA_H_LBDEM */
