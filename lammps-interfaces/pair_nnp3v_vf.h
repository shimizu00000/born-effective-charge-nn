/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(nnp3v_vf,PairNNP3v_vf)

#else

#ifndef LMP_PAIR_NNP3V_VF_H
#define LMP_PAIR_NNP3V_VF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNNP3v_vf : public Pair {
 public:
  PairNNP3v_vf(class LAMMPS *);
  virtual ~PairNNP3v_vf();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);

 private:
  bigint num_atoms;
  bigint num_atoms_a,num_atoms_b,num_atoms_c;

 protected:
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;              // added for TersoffZBL
    double ZBLcut,ZBLexpscale;
    double c5,ca1,ca4;           // added for TersoffMOD
    double powern_del;
    double c0;                   // added for TersoffMODC
  };

  Param *params;                // parameter set for an I-J-K interaction
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  double rc, eforce_z;
  int nweight_a,nweight_b,nweight_c;
  double *weight_a,*weight_b,*weight_c;
  int *g_a,*g_b,*g_c;
  int nlayer_a,nlayer_b,nlayer_c;
  int *layer_a,*layer_b,*layer_c;

  //(natom_x,num_g2_x_x)
  double *g2_a_a, *g2_a_b, *g2_a_c;
  double *g5_a_aa,*g5_a_ab,*g5_a_ac,
         *g5_a_bb,*g5_a_bc,
         *g5_a_cc;
  double *g2_b_a, *g2_b_b, *g2_b_c;
  double *g5_b_aa,*g5_b_ab,*g5_b_ac,
         *g5_b_bb,*g5_b_bc,
         *g5_b_cc;
  double *g2_c_a, *g2_c_b, *g2_c_c;
  double *g5_c_aa,*g5_c_ab,*g5_c_ac,
         *g5_c_bb,*g5_c_bc,
         *g5_c_cc;
  //(natom_x,natom,3,num_g2_x_x)
  double ***g2_deriv_a_a, ***g2_deriv_a_b, ***g2_deriv_a_c;
  double ***g5_deriv_a_aa,***g5_deriv_a_ab,***g5_deriv_a_ac,
         ***g5_deriv_a_bb,***g5_deriv_a_bc,
         ***g5_deriv_a_cc;
  double ***g2_deriv_b_a, ***g2_deriv_b_b, ***g2_deriv_b_c;
  double ***g5_deriv_b_aa,***g5_deriv_b_ab,***g5_deriv_b_ac,
         ***g5_deriv_b_bb,***g5_deriv_b_bc,
         ***g5_deriv_b_cc;
  double ***g2_deriv_c_a, ***g2_deriv_c_b, ***g2_deriv_c_c;
  double ***g5_deriv_c_aa,***g5_deriv_c_ab,***g5_deriv_c_ac,
         ***g5_deriv_c_bb,***g5_deriv_c_bc,
         ***g5_deriv_c_cc;
  double *g2_a_a_max, *g2_a_b_max, *g2_a_c_max;
  double *g5_a_aa_max,*g5_a_ab_max,*g5_a_ac_max,
         *g5_a_bb_max,*g5_a_bc_max,
         *g5_a_cc_max;
  double *g2_b_a_max, *g2_b_b_max, *g2_b_c_max;
  double *g5_b_aa_max,*g5_b_ab_max,*g5_b_ac_max,
         *g5_b_bb_max,*g5_b_bc_max,
         *g5_b_cc_max;
  double *g2_c_a_max, *g2_c_b_max, *g2_c_c_max;
  double *g5_c_aa_max,*g5_c_ab_max,*g5_c_ac_max,
         *g5_c_bb_max,*g5_c_bc_max,
         *g5_c_cc_max;
  double *g2_a_a_min, *g2_a_b_min, *g2_a_c_min;
  double *g5_a_aa_min,*g5_a_ab_min,*g5_a_ac_min,
         *g5_a_bb_min,*g5_a_bc_min,
         *g5_a_cc_min;
  double *g2_b_a_min, *g2_b_b_min, *g2_b_c_min;
  double *g5_b_aa_min,*g5_b_ab_min,*g5_b_ac_min,
         *g5_b_bb_min,*g5_b_bc_min,
         *g5_b_cc_min;
  double *g2_c_a_min, *g2_c_b_min, *g2_c_c_min;
  double *g5_c_aa_min,*g5_c_ab_min,*g5_c_ac_min,
         *g5_c_bb_min,*g5_c_bc_min,
         *g5_c_cc_min;
  double *eta2_a_a,*eta2_a_b,*eta2_a_c;
  double *eta2_b_a,*eta2_b_b,*eta2_b_c;
  double *eta2_c_a,*eta2_c_b,*eta2_c_c;
  double *rs2_a_a,*rs2_a_b,*rs2_a_c;
  double *rs2_b_a,*rs2_b_b,*rs2_b_c;
  double *rs2_c_a,*rs2_c_b,*rs2_c_c;
  double *eta5_a_aa,*eta5_a_ab,*eta5_a_ac,
         *eta5_a_bb,*eta5_a_bc,
         *eta5_a_cc;
  double *eta5_b_aa,*eta5_b_ab,*eta5_b_ac,
         *eta5_b_bb,*eta5_b_bc,
         *eta5_b_cc;
  double *eta5_c_aa,*eta5_c_ab,*eta5_c_ac,
         *eta5_c_bb,*eta5_c_bc,
         *eta5_c_cc;
  double *theta5_a_aa,*theta5_a_ab,*theta5_a_ac,
         *theta5_a_bb,*theta5_a_bc,
         *theta5_a_cc;
  double *theta5_b_aa,*theta5_b_ab,*theta5_b_ac,
         *theta5_b_bb,*theta5_b_bc,
         *theta5_b_cc;
  double *theta5_c_aa,*theta5_c_ab,*theta5_c_ac,
         *theta5_c_bb,*theta5_c_bc,
         *theta5_c_cc;
  double *zeta5_a_aa,*zeta5_a_ab,*zeta5_a_ac,
         *zeta5_a_bb,*zeta5_a_bc,
         *zeta5_a_cc;
  double *zeta5_b_aa,*zeta5_b_ab,*zeta5_b_ac,
         *zeta5_b_bb,*zeta5_b_bc,
         *zeta5_b_cc;
  double *zeta5_c_aa,*zeta5_c_ab,*zeta5_c_ac,
         *zeta5_c_bb,*zeta5_c_bc,
         *zeta5_c_cc;
  double *R_s_a_aa,*R_s_a_ab,*R_s_a_ac,
         *R_s_a_bb,*R_s_a_bc,
         *R_s_a_cc;
  double *R_s_b_aa,*R_s_b_ab,*R_s_b_ac,
         *R_s_b_bb,*R_s_b_bc,
         *R_s_b_cc;
  double *R_s_c_aa,*R_s_c_ab,*R_s_c_ac,
         *R_s_c_bb,*R_s_c_bc,
         *R_s_c_cc;
  double *input_a,*input_b,*input_c;
  double *dGdx_a,*dGdy_a,*dGdz_a;
  double *dGdx_b,*dGdy_b,*dGdz_b;
  double *dGdx_c,*dGdy_c,*dGdz_c;
  int    nhidden_a,nhidden_b,nhidden_c;
  double *hidden_a,*hidden_b,*hidden_c;

  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  void arraynnp();
  void sf_g2(int, int, double,
             double, double, double,
             double, double, double,
             double *&, double ***&, int,
             double *&, double *& );
  void sf_g5(int, int, int, double, double,
             double, double, double,
             double, double, double,
             double, double, double,
             double *&, double ***&, int,
             double *&, double *&, double *& );
  void hdnnp_e_a(double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, double & );
  void hdnnp_f_a(double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 int, double *&, double &, double &, double & );
  void hdnnp_e_b(double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, double & );
  void hdnnp_f_b(double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 int, double *&, double &, double &, double & );
  void hdnnp_e_c(double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, int,
                 double *&, double & );
  void hdnnp_f_c(double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 double ***&, int,
                 int, double *&, double &, double &, double & );

  virtual double fc(double);
  virtual double fc_deriv(double,double,double);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Tersoff requires atom IDs

This is a requirement to use the Tersoff potential.

E: Pair style Tersoff requires newton pair on

See the newton command.  This is a restriction to use the Tersoff
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/
