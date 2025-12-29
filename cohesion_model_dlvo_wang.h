/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Siyuan Wang, Peng Hou (School of Civil Engineering, Wuhan University)
   Attention: This script is only suitable for cgs unit.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   1. Van der Waals force:
   W_vdw = -H(D)/6*{2R1R2/(2R1+2R2+D)/D+2R1R2/(2R1+D)/(2R2+D)+ln[(2R1+2R2+D)D/(2R1+D)/(2R2+D)]}; //particle-particle
   W_vdw = -H(D)/6*{R/D+R/(2R+D)+ln[D/(2R+D)]}; //particle-wall
   f(D) = W_vdw/H(D);
   F_vdw = -dW_vdw/dD;
   F_vdw = -H'(D)*f(D)-H(D)*f'(D)
   2. Electrostatic force:
   F_edl = kappa*R1R2/(R1+R2)*Z*exp(-kappa*D); //particle-particle
   F_edl = kappa*R*Z*exp(-kappa*D); //particle-wall
   Z =64*pi*epsilon0*epsilon*(kT/e)^2*tanh(z*e*psi_i/4/k/T)*tanh(z*e*psi_j/4/k/T)

   D is s, seperated distance
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_DLVO_WANG,dlvo/wang,11)
#else
#ifndef COHESION_MODEL_DLVO_WANG_H_
#define COHESION_MODEL_DLVO_WANG_H_
#include "contact_models.h"
#include "math.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_DLVO_WANG> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), cohEnergyDens(NULL),
    sMin(NULL), AL(NULL)
    {
      if(comm->me==0) printf("COHESION_DLVO_WANG loaded\n");
      const double chiStiff = lmp->force->chiStiffnessScaling();
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
      registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model dlvo/wang");
      registry.registerProperty("sMin", &MODEL_PARAMS::createSMin);
      registry.connect("sMin", sMin, "cohesion_model dlvo/wang");
      registry.registerProperty("AL", &MODEL_PARAMS::createAL);
      registry.connect("AL", AL, "cohesion_model dlvo/wang");
      registry.registerProperty("EDLz", &MODEL_PARAMS::createEDLz);
      registry.connect("EDLz", EDLz, "cohesion_model dlvo/wang");
      registry.registerProperty("InvDebye", &MODEL_PARAMS::createInverseDebye);
      registry.connect("InvDebye", InvDebye, "cohesion_model dlvo/wang");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model dlvo/wang");

    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      //r is the distance between the sphere's centers
      const double r = cdata.r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double smin = sMin[cdata.itype][cdata.jtype]; //a minimum cutoff is used to avoit singularity.

      // Stifness scaling
      const double chiStiff = lmp->force->chiStiffnessScaling();
      const double realH_KD = cohEnergyDens[cdata.itype][cdata.jtype];
      const double realH_L = AL[cdata.itype][cdata.jtype];
      const double Inv_D = InvDebye[cdata.itype][cdata.jtype];
      const double Inter_Z = EDLz[cdata.itype][cdata.jtype];

      double Fn_coh;
      double f_s;
      double H_s;
      double f_sd; // "d" represents derivative formulation
      double H_sd;
      if(cdata.is_wall) {
        f_s = -(ri/smin+ri/(2*ri+smin)+log(smin/(2*ri+smin)))/6.0;
        H_s = (realH_KD*exp(-Inv_D*smin)+realH_L/(1.+14*smin/100e-7))/chiStiff;
        f_sd = ri/6.0/smin/smin;
        H_sd = -Inv_D*exp(-Inv_D*smin)*realH_KD/chiStiff-14/(100e-7)/pow(1.+14*smin/100e-7,2)*realH_L/chiStiff;

        Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*Inter_Z*exp(-Inv_D*smin)/chiStiff;

      } else {
        f_s = -(2*ri*rj/(2*ri+2*rj+smin)/smin+2*ri*rj/(2*ri+smin)/(2*rj+smin)+log((2*ri+2*rj+smin)*smin/(2*ri+smin)/(2*rj+smin)))/6.0;
        H_s = (realH_KD*exp(-Inv_D*smin)+realH_L/(1.+11*smin/100e-7))/chiStiff;
        f_sd = (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin))
           * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)
           * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)/3.0;
        H_sd = -Inv_D*exp(-Inv_D*smin)*realH_KD/chiStiff-11/(100e-7)/pow(1.+11*smin/100e-7,2)*realH_L/chiStiff;

        Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*rj/(ri+rj)*Inter_Z*exp(-Inv_D*smin)/chiStiff;
      }

      // Question for normal force: cohesion force should not impact for tangential force, therefore we remove the following statement,
      // cdata.Fn += Fn_coh;

      // apply normal force
      if(cdata.is_wall) {
        i_forces.delta_F[0] += Fn_coh * cdata.en[0];
        i_forces.delta_F[1] += Fn_coh * cdata.en[1];
        i_forces.delta_F[2] += Fn_coh * cdata.en[2];
      }
      else {
        const double fx = Fn_coh * cdata.en[0];
        const double fy = Fn_coh * cdata.en[1];
        const double fz = Fn_coh * cdata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
      }

      if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    void noCollision(ContactData& cdata, ForceData & i_forces, ForceData & j_forces)
    {
      double **v = atom->v;
      int *type = atom->type;
      //r is the distance between the sphere's centers
      const double r = sqrt(cdata.rsq);
      const double rinv = 1/r;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      double s, smax;
      if(cdata.is_wall) {
        s = r - ri; // separating distance between the surfaces of the two interacting particles
        smax = ri / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the DLVO force is not considered.
        //smax = ri; //XL
      } else{
        s = r - (ri + rj); // separating distance between the surfaces of the two interacting particles
        smax = (ri + rj) / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the DLVO force is not considered.
        //smax = ri*rj/(ri + rj); //XL
      }

      //only particle-particle interactions are considered here
      const int itype = type[cdata.i];
      const int jtype = type[cdata.j];
      const double smin = sMin[itype][jtype]; //1 nm Yu2012 Chemial Engineering Science

      // Stifness scaling
      const double chiStiff = lmp->force->chiStiffnessScaling();
      const double realH_KD = cohEnergyDens[cdata.itype][cdata.jtype];
      const double realH_L = AL[cdata.itype][cdata.jtype];
      const double Inv_D = InvDebye[cdata.itype][cdata.jtype];
      const double Inter_Z = EDLz[cdata.itype][cdata.jtype];

      double Fn_coh;
      double f_s;
      double H_s;
      double f_sd; // "d" represents derivative formulation
      double H_sd;

      if (s < smax) {
        if (s < smin) {
          if(cdata.is_wall) {
            f_s = -(ri/smin+ri/(2*ri+smin)+log(smin/(2*ri+smin)))/6.0;
            H_s = (realH_KD*exp(-Inv_D*smin)+realH_L/(1.+14*smin/100e-7))/chiStiff;
            f_sd = ri/6.0/smin/smin;
            H_sd = -Inv_D*exp(-Inv_D*smin)*realH_KD/chiStiff-14/(100e-7)/pow(1.+14*smin/100e-7,2)*realH_L/chiStiff;

            Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*Inter_Z*exp(-Inv_D*smin);
          } else {
            f_s = -(2*ri*rj/(2*ri+2*rj+smin)/smin+2*ri*rj/(2*ri+smin)/(2*rj+smin)+log((2*ri+2*rj+smin)*smin/(2*ri+smin)/(2*rj+smin)))/6.0;
            H_s = (realH_KD*exp(-Inv_D*smin)+realH_L/(1.+11*smin/100e-7))/chiStiff;
            f_sd = (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin))
              * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)
              * ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)/3.0;
            H_sd = -Inv_D*exp(-Inv_D*smin)*realH_KD/chiStiff-11/(100e-7)/pow(1.+11*smin/100e-7,2)*realH_L/chiStiff;

            Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*rj/(ri+rj)*Inter_Z*exp(-Inv_D*smin);
          }
        } else {
          if(cdata.is_wall) {
            f_s = -(ri/s+ri/(2*ri+s)+log(s/(2*ri+s)))/6.0;
            H_s = (realH_KD*exp(-Inv_D*s)+realH_L/(1.+14*s/100e-7));
            f_sd = ri/6.0/s/s;
            H_sd = -Inv_D*exp(-Inv_D*s)*realH_KD-14/(100e-7)/pow(1.+14*s/100e-7,2)*realH_L;

            Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*Inter_Z*exp(-Inv_D*s);
          } else {
            f_s = -(2*ri*rj/(2*ri+2*rj+s)/s+2*ri*rj/(2*ri+s)/(2*rj+s)+log((2*ri+2*rj+s)*s/(2*ri+s)/(2*rj+s)))/6.0;
            H_s = (realH_KD*exp(-Inv_D*s)+realH_L/(1.+11*s/100e-7));
            f_sd = (2.0*ri*rj*(ri+rj+s))/(s*s*(2*ri + 2*rj + s)*(2*ri + 2*rj + s))
              * ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1)
              * ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1)/3.0;
            H_sd = -Inv_D*exp(-Inv_D*s)*realH_KD-11/(100e-7)/pow(1.+11*s/100e-7,2)*realH_L;

            Fn_coh = -(f_s*H_sd)-(f_sd*H_s)+Inv_D*ri*rj/(ri+rj)*Inter_Z*exp(-Inv_D*s);
          }
        }

        if(cdata.is_wall) {
          i_forces.delta_F[0] += Fn_coh * cdata.delta[0] * rinv;
          i_forces.delta_F[1] += Fn_coh * cdata.delta[1] * rinv;
          i_forces.delta_F[2] += Fn_coh * cdata.delta[2] * rinv;
        } else {
          const double fx = Fn_coh * cdata.delta[0] * rinv;
          const double fy = Fn_coh * cdata.delta[1] * rinv;
          const double fz = Fn_coh * cdata.delta[2] * rinv;

          i_forces.delta_F[0] += fx;
          i_forces.delta_F[1] += fy;
          i_forces.delta_F[2] += fz;

          j_forces.delta_F[0] -= fx;
          j_forces.delta_F[1] -= fy;
          j_forces.delta_F[2] -= fz;

          cdata.has_force_update = true;
          if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
        }
      }
    }


  private:
    double ** cohEnergyDens;
    double ** sMin;
    double ** AL;
    double ** EDLz;
    double ** InvDebye;
  };
}
}


#endif // COHESION_MODEL_DLVO_WANG_H_
#endif

