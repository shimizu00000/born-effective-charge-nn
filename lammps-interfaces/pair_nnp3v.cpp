/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_nnp3v.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 2048 //1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairNNP3v::PairNNP3v(LAMMPS *lmp) : Pair(lmp)
{

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  maxshort = 10;
  neighshort = NULL;

  g_a = NULL;
  g_b = NULL;
  g_c = NULL;

  layer_a = NULL;
  layer_b = NULL;
  layer_c = NULL;

  weight_a = NULL;
  weight_b = NULL;
  weight_c = NULL;

  eta2_a_a = NULL;
  eta2_a_b = NULL;
  eta2_a_c = NULL;
  rs2_a_a = NULL;
  rs2_a_b = NULL;
  rs2_a_c = NULL;
  eta5_a_aa = NULL;
  eta5_a_ab = NULL;
  eta5_a_ac = NULL;
  eta5_a_bb = NULL;
  eta5_a_bc = NULL;
  eta5_a_cc = NULL;
  theta5_a_aa = NULL;
  theta5_a_ab = NULL;
  theta5_a_ac = NULL;
  theta5_a_bb = NULL;
  theta5_a_bc = NULL;
  theta5_a_cc = NULL;
  zeta5_a_aa = NULL;
  zeta5_a_ab = NULL;
  zeta5_a_ac = NULL;
  zeta5_a_bb = NULL;
  zeta5_a_bc = NULL;
  zeta5_a_cc = NULL;
  lambda5_a_aa = NULL;
  lambda5_a_ab = NULL;
  lambda5_a_ac = NULL;
  lambda5_a_bb = NULL;
  lambda5_a_bc = NULL;
  lambda5_a_cc = NULL;

  eta2_b_a = NULL;
  eta2_b_b = NULL;
  eta2_b_c = NULL;
  rs2_b_a = NULL;
  rs2_b_b = NULL;
  rs2_b_c = NULL;
  eta5_b_aa = NULL;
  eta5_b_ab = NULL;
  eta5_b_ac = NULL;
  eta5_b_bb = NULL;
  eta5_b_bc = NULL;
  eta5_b_cc = NULL;
  theta5_b_aa = NULL;
  theta5_b_ab = NULL;
  theta5_b_ac = NULL;
  theta5_b_bb = NULL;
  theta5_b_bc = NULL;
  theta5_b_cc = NULL;
  zeta5_b_aa = NULL;
  zeta5_b_ab = NULL;
  zeta5_b_ac = NULL;
  zeta5_b_bb = NULL;
  zeta5_b_bc = NULL;
  zeta5_b_cc = NULL;
  lambda5_b_aa = NULL;
  lambda5_b_ab = NULL;
  lambda5_b_ac = NULL;
  lambda5_b_bb = NULL;
  lambda5_b_bc = NULL;
  lambda5_b_cc = NULL;

  eta2_c_a = NULL;
  eta2_c_b = NULL;
  eta2_c_c = NULL;
  rs2_c_a = NULL;
  rs2_c_b = NULL;
  rs2_c_c = NULL;
  eta5_c_aa = NULL;
  eta5_c_ab = NULL;
  eta5_c_ac = NULL;
  eta5_c_bb = NULL;
  eta5_c_bc = NULL;
  eta5_c_cc = NULL;
  theta5_c_aa = NULL;
  theta5_c_ab = NULL;
  theta5_c_ac = NULL;
  theta5_c_bb = NULL;
  theta5_c_bc = NULL;
  theta5_c_cc = NULL;
  zeta5_c_aa = NULL;
  zeta5_c_ab = NULL;
  zeta5_c_ac = NULL;
  zeta5_c_bb = NULL;
  zeta5_c_bc = NULL;
  zeta5_c_cc = NULL;
  lambda5_c_aa = NULL;
  lambda5_c_ab = NULL;
  lambda5_c_ac = NULL;
  lambda5_c_bb = NULL;
  lambda5_c_bc = NULL;
  lambda5_c_cc = NULL;

  g2_a_a_max = NULL;
  g2_a_b_max = NULL;
  g2_a_c_max = NULL;
  g5_a_aa_max = NULL;
  g5_a_ab_max = NULL;
  g5_a_ac_max = NULL;
  g5_a_bb_max = NULL;
  g5_a_bc_max = NULL;
  g5_a_cc_max = NULL;

  g2_b_a_max = NULL;
  g2_b_b_max = NULL;
  g2_b_c_max = NULL;
  g5_b_aa_max = NULL;
  g5_b_ab_max = NULL;
  g5_b_ac_max = NULL;
  g5_b_bb_max = NULL;
  g5_b_bc_max = NULL;
  g5_b_cc_max = NULL;

  g2_c_a_max = NULL;
  g2_c_b_max = NULL;
  g2_c_c_max = NULL;
  g5_c_aa_max = NULL;
  g5_c_ab_max = NULL;
  g5_c_ac_max = NULL;
  g5_c_bb_max = NULL;
  g5_c_bc_max = NULL;
  g5_c_cc_max = NULL;

  g2_a_a_min = NULL;
  g2_a_b_min = NULL;
  g2_a_c_min = NULL;
  g5_a_aa_min = NULL;
  g5_a_ab_min = NULL;
  g5_a_ac_min = NULL;
  g5_a_bb_min = NULL;
  g5_a_bc_min = NULL;
  g5_a_cc_min = NULL;

  g2_b_a_min = NULL;
  g2_b_b_min = NULL;
  g2_b_c_min = NULL;
  g5_b_aa_min = NULL;
  g5_b_ab_min = NULL;
  g5_b_ac_min = NULL;
  g5_b_bb_min = NULL;
  g5_b_bc_min = NULL;
  g5_b_cc_min = NULL;

  g2_c_a_min = NULL;
  g2_c_b_min = NULL;
  g2_c_c_min = NULL;
  g5_c_aa_min = NULL;
  g5_c_ab_min = NULL;
  g5_c_ac_min = NULL;
  g5_c_bb_min = NULL;
  g5_c_bc_min = NULL;
  g5_c_cc_min = NULL;

  g2_a_a = NULL;
  g2_a_b = NULL;
  g2_a_c = NULL;
  g5_a_aa = NULL;
  g5_a_ab = NULL;
  g5_a_ac = NULL;
  g5_a_bb = NULL;
  g5_a_bc = NULL;
  g5_a_cc = NULL;

  g2_b_a = NULL;
  g2_b_b = NULL;
  g2_b_c = NULL;
  g5_b_aa = NULL;
  g5_b_ab = NULL;
  g5_b_ac = NULL;
  g5_b_bb = NULL;
  g5_b_bc = NULL;
  g5_b_cc = NULL;

  g2_c_a = NULL;
  g2_c_b = NULL;
  g2_c_c = NULL;
  g5_c_aa = NULL;
  g5_c_ab = NULL;
  g5_c_ac = NULL;
  g5_c_bb = NULL;
  g5_c_bc = NULL;
  g5_c_cc = NULL;

  g2_deriv_a_a = NULL;
  g2_deriv_a_b = NULL;
  g2_deriv_a_c = NULL;
  g5_deriv_a_aa = NULL;
  g5_deriv_a_ab = NULL;
  g5_deriv_a_ac = NULL;
  g5_deriv_a_bb = NULL;
  g5_deriv_a_bc = NULL;
  g5_deriv_a_cc = NULL;

  g2_deriv_b_a = NULL;
  g2_deriv_b_b = NULL;
  g2_deriv_b_c = NULL;
  g5_deriv_b_aa = NULL;
  g5_deriv_b_ab = NULL;
  g5_deriv_b_ac = NULL;
  g5_deriv_b_bb = NULL;
  g5_deriv_b_bc = NULL;
  g5_deriv_b_cc = NULL;

  g2_deriv_c_a = NULL;
  g2_deriv_c_b = NULL;
  g2_deriv_c_c = NULL;
  g5_deriv_c_aa = NULL;
  g5_deriv_c_ab = NULL;
  g5_deriv_c_ac = NULL;
  g5_deriv_c_bb = NULL;
  g5_deriv_c_bc = NULL;
  g5_deriv_c_cc = NULL;

  hidden_a = NULL;
  hidden_b = NULL;
  hidden_c = NULL;

  input_a = NULL;
  input_b = NULL;
  input_c = NULL;

  dGdx_a = NULL;
  dGdy_a = NULL;
  dGdz_a = NULL;

  dGdx_b = NULL;
  dGdy_b = NULL;
  dGdz_b = NULL;

  dGdx_c = NULL;
  dGdy_c = NULL;
  dGdz_c = NULL;


}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairNNP3v::~PairNNP3v()
{

  if (copymode) return;

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];

  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  memory->destroy(g_a);
  memory->destroy(g_b);
  memory->destroy(g_c);

  memory->destroy(layer_a);
  memory->destroy(layer_b);
  memory->destroy(layer_c);

  memory->destroy(weight_a);
  memory->destroy(weight_b);
  memory->destroy(weight_c);

  memory->destroy(eta2_a_a);
  memory->destroy(eta2_a_b);
  memory->destroy(eta2_a_c);

  memory->destroy(eta2_b_a);
  memory->destroy(eta2_b_b);
  memory->destroy(eta2_b_c);

  memory->destroy(eta2_c_a);
  memory->destroy(eta2_c_b);
  memory->destroy(eta2_c_c);

  memory->destroy(rs2_a_a);
  memory->destroy(rs2_a_b);
  memory->destroy(rs2_a_c);

  memory->destroy(rs2_b_a);
  memory->destroy(rs2_b_b);
  memory->destroy(rs2_b_c);

  memory->destroy(rs2_c_a);
  memory->destroy(rs2_c_b);
  memory->destroy(rs2_c_c);

  memory->destroy(eta5_a_aa);
  memory->destroy(eta5_a_ab);
  memory->destroy(eta5_a_ac);
  memory->destroy(eta5_a_bb);
  memory->destroy(eta5_a_bc);
  memory->destroy(eta5_a_cc);

  memory->destroy(eta5_b_aa);
  memory->destroy(eta5_b_ab);
  memory->destroy(eta5_b_ac);
  memory->destroy(eta5_b_bb);
  memory->destroy(eta5_b_bc);
  memory->destroy(eta5_b_cc);

  memory->destroy(eta5_c_aa);
  memory->destroy(eta5_c_ab);
  memory->destroy(eta5_c_ac);
  memory->destroy(eta5_c_bb);
  memory->destroy(eta5_c_bc);
  memory->destroy(eta5_c_cc);

  memory->destroy(theta5_a_aa);
  memory->destroy(theta5_a_ab);
  memory->destroy(theta5_a_ac);
  memory->destroy(theta5_a_bb);
  memory->destroy(theta5_a_bc);
  memory->destroy(theta5_a_cc);

  memory->destroy(theta5_b_aa);
  memory->destroy(theta5_b_ab);
  memory->destroy(theta5_b_ac);
  memory->destroy(theta5_b_bb);
  memory->destroy(theta5_b_bc);
  memory->destroy(theta5_b_cc);

  memory->destroy(theta5_c_aa);
  memory->destroy(theta5_c_ab);
  memory->destroy(theta5_c_ac);
  memory->destroy(theta5_c_bb);
  memory->destroy(theta5_c_bc);
  memory->destroy(theta5_c_cc);

  memory->destroy(zeta5_a_aa);
  memory->destroy(zeta5_a_ab);
  memory->destroy(zeta5_a_ac);
  memory->destroy(zeta5_a_bb);
  memory->destroy(zeta5_a_bc);
  memory->destroy(zeta5_a_cc);

  memory->destroy(zeta5_b_aa);
  memory->destroy(zeta5_b_ab);
  memory->destroy(zeta5_b_ac);
  memory->destroy(zeta5_b_bb);
  memory->destroy(zeta5_b_bc);
  memory->destroy(zeta5_b_cc);

  memory->destroy(zeta5_c_aa);
  memory->destroy(zeta5_c_ab);
  memory->destroy(zeta5_c_ac);
  memory->destroy(zeta5_c_bb);
  memory->destroy(zeta5_c_bc);
  memory->destroy(zeta5_c_cc);

  memory->destroy(lambda5_a_aa);
  memory->destroy(lambda5_a_ab);
  memory->destroy(lambda5_a_ac);
  memory->destroy(lambda5_a_bb);
  memory->destroy(lambda5_a_bc);
  memory->destroy(lambda5_a_cc);

  memory->destroy(lambda5_b_aa);
  memory->destroy(lambda5_b_ab);
  memory->destroy(lambda5_b_ac);
  memory->destroy(lambda5_b_bb);
  memory->destroy(lambda5_b_bc);
  memory->destroy(lambda5_b_cc);

  memory->destroy(lambda5_c_aa);
  memory->destroy(lambda5_c_ab);
  memory->destroy(lambda5_c_ac);
  memory->destroy(lambda5_c_bb);
  memory->destroy(lambda5_c_bc);
  memory->destroy(lambda5_c_cc);

  memory->destroy(g2_a_a_max);
  memory->destroy(g2_a_b_max);
  memory->destroy(g2_a_c_max);
  memory->destroy(g5_a_aa_max);
  memory->destroy(g5_a_ab_max);
  memory->destroy(g5_a_ac_max);
  memory->destroy(g5_a_bb_max);
  memory->destroy(g5_a_bc_max);
  memory->destroy(g5_a_cc_max);

  memory->destroy(g2_b_a_max);
  memory->destroy(g2_b_b_max);
  memory->destroy(g2_b_c_max);
  memory->destroy(g5_b_aa_max);
  memory->destroy(g5_b_ab_max);
  memory->destroy(g5_b_ac_max);
  memory->destroy(g5_b_bb_max);
  memory->destroy(g5_b_bc_max);
  memory->destroy(g5_b_cc_max);

  memory->destroy(g2_c_a_max);
  memory->destroy(g2_c_b_max);
  memory->destroy(g2_c_c_max);
  memory->destroy(g5_c_aa_max);
  memory->destroy(g5_c_ab_max);
  memory->destroy(g5_c_ac_max);
  memory->destroy(g5_c_bb_max);
  memory->destroy(g5_c_bc_max);
  memory->destroy(g5_c_cc_max);

  memory->destroy(g2_a_a_min);
  memory->destroy(g2_a_b_min);
  memory->destroy(g2_a_c_min);
  memory->destroy(g5_a_aa_min);
  memory->destroy(g5_a_ab_min);
  memory->destroy(g5_a_ac_min);
  memory->destroy(g5_a_bb_min);
  memory->destroy(g5_a_bc_min);
  memory->destroy(g5_a_cc_min);

  memory->destroy(g2_b_a_min);
  memory->destroy(g2_b_b_min);
  memory->destroy(g2_b_c_min);
  memory->destroy(g5_b_aa_min);
  memory->destroy(g5_b_ab_min);
  memory->destroy(g5_b_ac_min);
  memory->destroy(g5_b_bb_min);
  memory->destroy(g5_b_bc_min);
  memory->destroy(g5_b_cc_min);

  memory->destroy(g2_c_a_min);
  memory->destroy(g2_c_b_min);
  memory->destroy(g2_c_c_min);
  memory->destroy(g5_c_aa_min);
  memory->destroy(g5_c_ab_min);
  memory->destroy(g5_c_ac_min);
  memory->destroy(g5_c_bb_min);
  memory->destroy(g5_c_bc_min);
  memory->destroy(g5_c_cc_min);

  memory->destroy(g2_a_a);
  memory->destroy(g2_a_b);
  memory->destroy(g2_a_c);
  memory->destroy(g5_a_aa);
  memory->destroy(g5_a_ab);
  memory->destroy(g5_a_ac);
  memory->destroy(g5_a_bb);
  memory->destroy(g5_a_bc);
  memory->destroy(g5_a_cc);

  memory->destroy(g2_b_a);
  memory->destroy(g2_b_b);
  memory->destroy(g2_b_c);
  memory->destroy(g5_b_aa);
  memory->destroy(g5_b_ab);
  memory->destroy(g5_b_ac);
  memory->destroy(g5_b_bb);
  memory->destroy(g5_b_bc);
  memory->destroy(g5_b_cc);

  memory->destroy(g2_c_a);
  memory->destroy(g2_c_b);
  memory->destroy(g2_c_c);
  memory->destroy(g5_c_aa);
  memory->destroy(g5_c_ab);
  memory->destroy(g5_c_ac);
  memory->destroy(g5_c_bb);
  memory->destroy(g5_c_bc);
  memory->destroy(g5_c_cc);

  memory->destroy(g2_deriv_a_a);
  memory->destroy(g2_deriv_a_b);
  memory->destroy(g2_deriv_a_c);
  memory->destroy(g5_deriv_a_aa);
  memory->destroy(g5_deriv_a_ab);
  memory->destroy(g5_deriv_a_ac);
  memory->destroy(g5_deriv_a_bb);
  memory->destroy(g5_deriv_a_bc);
  memory->destroy(g5_deriv_a_cc);

  memory->destroy(g2_deriv_b_a);
  memory->destroy(g2_deriv_b_b);
  memory->destroy(g2_deriv_b_c);
  memory->destroy(g5_deriv_b_aa);
  memory->destroy(g5_deriv_b_ab);
  memory->destroy(g5_deriv_b_ac);
  memory->destroy(g5_deriv_b_bb);
  memory->destroy(g5_deriv_b_bc);
  memory->destroy(g5_deriv_b_cc);

  memory->destroy(g2_deriv_c_a);
  memory->destroy(g2_deriv_c_b);
  memory->destroy(g2_deriv_c_c);
  memory->destroy(g5_deriv_c_aa);
  memory->destroy(g5_deriv_c_ab);
  memory->destroy(g5_deriv_c_ac);
  memory->destroy(g5_deriv_c_bb);
  memory->destroy(g5_deriv_c_bc);
  memory->destroy(g5_deriv_c_cc);

  memory->destroy(hidden_a);
  memory->destroy(hidden_b);
  memory->destroy(hidden_c);

  memory->destroy(input_a);
  memory->destroy(input_b);
  memory->destroy(input_c);

  memory->destroy(dGdx_a);
  memory->destroy(dGdy_a);
  memory->destroy(dGdz_a);

  memory->destroy(dGdx_b);
  memory->destroy(dGdy_b);
  memory->destroy(dGdz_b);

  memory->destroy(dGdx_c);
  memory->destroy(dGdy_c);
  memory->destroy(dGdz_c);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
    delete [] map;
  }


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::compute(int eflag, int vflag)
{
  int i,j,k,n,ii,jj,kk,nn,inum,jnum;
  int itype,jtype,ktype;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,srsq1,srsq2;
  double delr1[3],delr2[3],fx,fy,fz;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal; // # of atoms / core
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

  inum = list->inum; // # of atoms
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int nall = atom->nlocal + atom->nghost;
  memory->create(g2_deriv_a_a,nall,3,g_a[0],"pair:nnp1:g2_deriv_a_a");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[0]; j++) {
        g2_deriv_a_a[i][0][j] = 0.0;
        g2_deriv_a_a[i][1][j] = 0.0;
        g2_deriv_a_a[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_a_b,nall,3,g_a[1],"pair:nnp1:g2_deriv_a_b");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[1]; j++) {
        g2_deriv_a_b[i][0][j] = 0.0;
        g2_deriv_a_b[i][1][j] = 0.0;
        g2_deriv_a_b[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_a_c,nall,3,g_a[2],"pair:nnp1:g2_deriv_a_c");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[2]; j++) {
        g2_deriv_a_c[i][0][j] = 0.0;
        g2_deriv_a_c[i][1][j] = 0.0;
        g2_deriv_a_c[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_aa,nall,3,g_a[3],"pair:nnp1:g5_deriv_a_aa");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[3]; j++) {
        g5_deriv_a_aa[i][0][j] = 0.0;
        g5_deriv_a_aa[i][1][j] = 0.0;
        g5_deriv_a_aa[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_ab,nall,3,g_a[4],"pair:nnp1:g5_deriv_a_ab");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[4]; j++) {
        g5_deriv_a_ab[i][0][j] = 0.0;
        g5_deriv_a_ab[i][1][j] = 0.0;
        g5_deriv_a_ab[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_ac,nall,3,g_a[5],"pair:nnp1:g5_deriv_a_ac");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[5]; j++) {
        g5_deriv_a_ac[i][0][j] = 0.0;
        g5_deriv_a_ac[i][1][j] = 0.0;
        g5_deriv_a_ac[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_bb,nall,3,g_a[6],"pair:nnp1:g5_deriv_a_bb");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[6]; j++) {
        g5_deriv_a_bb[i][0][j] = 0.0;
        g5_deriv_a_bb[i][1][j] = 0.0;
        g5_deriv_a_bb[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_bc,nall,3,g_a[7],"pair:nnp1:g5_deriv_a_bc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[7]; j++) {
        g5_deriv_a_bc[i][0][j] = 0.0;
        g5_deriv_a_bc[i][1][j] = 0.0;
        g5_deriv_a_bc[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_a_cc,nall,3,g_a[8],"pair:nnp1:g5_deriv_a_cc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_a[8]; j++) {
        g5_deriv_a_cc[i][0][j] = 0.0;
        g5_deriv_a_cc[i][1][j] = 0.0;
        g5_deriv_a_cc[i][2][j] = 0.0;
  }}

  memory->create(g2_deriv_b_a,nall,3,g_b[0],"pair:nnp1:g2_deriv_b_a");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[0]; j++) {
        g2_deriv_b_a[i][0][j] = 0.0;
        g2_deriv_b_a[i][1][j] = 0.0;
        g2_deriv_b_a[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_b_b,nall,3,g_b[1],"pair:nnp1:g2_deriv_b_b");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[1]; j++) {
        g2_deriv_b_b[i][0][j] = 0.0;
        g2_deriv_b_b[i][1][j] = 0.0;
        g2_deriv_b_b[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_b_c,nall,3,g_b[2],"pair:nnp1:g2_deriv_b_c");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[2]; j++) {
        g2_deriv_b_c[i][0][j] = 0.0;
        g2_deriv_b_c[i][1][j] = 0.0;
        g2_deriv_b_c[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_aa,nall,3,g_b[3],"pair:nnp1:g5_deriv_b_aa");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[3]; j++) {
        g5_deriv_b_aa[i][0][j] = 0.0;
        g5_deriv_b_aa[i][1][j] = 0.0;
        g5_deriv_b_aa[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_ab,nall,3,g_b[4],"pair:nnp1:g5_deriv_b_ab");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[4]; j++) {
        g5_deriv_b_ab[i][0][j] = 0.0;
        g5_deriv_b_ab[i][1][j] = 0.0;
        g5_deriv_b_ab[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_ac,nall,3,g_b[5],"pair:nnp1:g5_deriv_b_ac");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[5]; j++) {
        g5_deriv_b_ac[i][0][j] = 0.0;
        g5_deriv_b_ac[i][1][j] = 0.0;
        g5_deriv_b_ac[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_bb,nall,3,g_b[6],"pair:nnp1:g5_deriv_b_bb");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[6]; j++) {
        g5_deriv_b_bb[i][0][j] = 0.0;
        g5_deriv_b_bb[i][1][j] = 0.0;
        g5_deriv_b_bb[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_bc,nall,3,g_b[7],"pair:nnp1:g5_deriv_b_bc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[7]; j++) {
        g5_deriv_b_bc[i][0][j] = 0.0;
        g5_deriv_b_bc[i][1][j] = 0.0;
        g5_deriv_b_bc[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_b_cc,nall,3,g_b[8],"pair:nnp1:g5_deriv_b_cc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_b[8]; j++) {
        g5_deriv_b_cc[i][0][j] = 0.0;
        g5_deriv_b_cc[i][1][j] = 0.0;
        g5_deriv_b_cc[i][2][j] = 0.0;
  }}

  memory->create(g2_deriv_c_a,nall,3,g_c[0],"pair:nnp1:g2_deriv_c_a");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[0]; j++) {
        g2_deriv_c_a[i][0][j] = 0.0;
        g2_deriv_c_a[i][1][j] = 0.0;
        g2_deriv_c_a[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_c_b,nall,3,g_c[1],"pair:nnp1:g2_deriv_c_b");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[1]; j++) {
        g2_deriv_c_b[i][0][j] = 0.0;
        g2_deriv_c_b[i][1][j] = 0.0;
        g2_deriv_c_b[i][2][j] = 0.0;
  }}
  memory->create(g2_deriv_c_c,nall,3,g_c[2],"pair:nnp1:g2_deriv_c_c");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[2]; j++) {
        g2_deriv_c_c[i][0][j] = 0.0;
        g2_deriv_c_c[i][1][j] = 0.0;
        g2_deriv_c_c[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_aa,nall,3,g_c[3],"pair:nnp1:g5_deriv_c_aa");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[3]; j++) {
        g5_deriv_c_aa[i][0][j] = 0.0;
        g5_deriv_c_aa[i][1][j] = 0.0;
        g5_deriv_c_aa[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_ab,nall,3,g_c[4],"pair:nnp1:g5_deriv_c_ab");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[4]; j++) {
        g5_deriv_c_ab[i][0][j] = 0.0;
        g5_deriv_c_ab[i][1][j] = 0.0;
        g5_deriv_c_ab[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_ac,nall,3,g_c[5],"pair:nnp1:g5_deriv_c_ac");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[5]; j++) {
        g5_deriv_c_ac[i][0][j] = 0.0;
        g5_deriv_c_ac[i][1][j] = 0.0;
        g5_deriv_c_ac[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_bb,nall,3,g_c[6],"pair:nnp1:g5_deriv_c_bb");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[6]; j++) {
        g5_deriv_c_bb[i][0][j] = 0.0;
        g5_deriv_c_bb[i][1][j] = 0.0;
        g5_deriv_c_bb[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_bc,nall,3,g_c[7],"pair:nnp1:g5_deriv_c_bc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[7]; j++) {
        g5_deriv_c_bc[i][0][j] = 0.0;
        g5_deriv_c_bc[i][1][j] = 0.0;
        g5_deriv_c_bc[i][2][j] = 0.0;
  }}
  memory->create(g5_deriv_c_cc,nall,3,g_c[8],"pair:nnp1:g5_deriv_c_cc");
  for (i = 0; i < nall; i++) {
      for (j = 0; j < g_c[8]; j++) {
        g5_deriv_c_cc[i][0][j] = 0.0;
        g5_deriv_c_cc[i][1][j] = 0.0;
        g5_deriv_c_cc[i][2][j] = 0.0;
  }}


  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]]; // element
    xtmp = x[i][0]; ytmp = x[i][1]; ztmp = x[i][2];

    // SF : array size(num_g2), SF_deriv : array size(natom,3,num_g2)
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK; //???
      delx = xtmp - x[j][0]; dely = ytmp - x[j][1]; delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz; // rc**2
      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }
    } // for jj

    // three-body interactions, skip immediately if I-J is not within cutoff
    for (jj = 0; jj < numshort; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      delr1[0] = x[j][0] - xtmp; delr1[1] = x[j][1] - ytmp; delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= cutshortsq) continue;
      srsq1 = sqrt(rsq1);

      // A-
      if (itype == 0) {
        // A-A
        if(jtype == 0) {
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_a_a, g2_deriv_a_a, g_a[0],
                eta2_a_a, rs2_a_a); 
        }
        // A-B
        else if(jtype == 1) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_a_b, g2_deriv_a_b, g_a[1],
                eta2_a_b, rs2_a_b); 
        }
        // A-C
        else if(jtype == 2) {
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_a_c, g2_deriv_a_c, g_a[2],
                eta2_a_c, rs2_a_c); 
        }
      } // A-
      else if(itype == 1) {
        // B-A
        if(jtype == 0) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_b_a, g2_deriv_b_a, g_b[0],
                eta2_b_a, rs2_b_a); 
        }
        // B-B
        else if(jtype == 1) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_b_b, g2_deriv_b_b, g_b[1],
                eta2_b_b, rs2_b_b); 
        }
        // B-C
        else if(jtype == 2) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_b_c, g2_deriv_b_c, g_b[2],
                eta2_b_c, rs2_b_c); 
        }
      } // B-
      else if(itype == 2) {
        // C-A
        if(jtype == 0) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_c_a, g2_deriv_c_a, g_c[0],
                eta2_c_a, rs2_c_a); 
        }
        // C-B 
        else if(jtype == 1) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_c_b, g2_deriv_c_b, g_c[1],
                eta2_c_b, rs2_c_b);  
        }
        // C-C
        else if(jtype == 2) { 
          sf_g2(i,j,srsq1,
                x[i][0],x[i][1],x[i][2],
                x[j][0],x[j][1],x[j][2],
                g2_c_c, g2_deriv_c_c, g_c[2],
                eta2_c_c, rs2_c_c); 
        }
      } // C-

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= cutshortsq) continue;
        srsq2 = sqrt(rsq2);
 
        // A-
        if (itype == 0) {
          // A-AA
          if (jtype == 0 && ktype == 0) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_aa,g5_deriv_a_aa,g_a[3],
                  eta5_a_aa,lambda5_a_aa,zeta5_a_aa); 
          }
          // A-AB
          else if (jtype == 0 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_ab,g5_deriv_a_ab,g_a[4],
                  eta5_a_ab,lambda5_a_ab,zeta5_a_ab); 
          }
          // A-AC
          else if (jtype == 0 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_ac,g5_deriv_a_ac,g_a[5],
                  eta5_a_ac,lambda5_a_ac,zeta5_a_ac); 
          }
          // A-BB
          else if (jtype == 1 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_bb,g5_deriv_a_bb,g_a[6],
                  eta5_a_bb,lambda5_a_bb,zeta5_a_bb); 
          }
          // A-BC
          else if (jtype == 1 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_bc,g5_deriv_a_bc,g_a[7],
                  eta5_a_bc,lambda5_a_bc,zeta5_a_bc); 
          }
          // A-CC
          else if (jtype == 2 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_a_cc,g5_deriv_a_cc,g_a[8],
                  eta5_a_cc,lambda5_a_cc,zeta5_a_cc); 
          }
          // A-BA,-CA,-CB
          else continue;
        } // A-
        // B-
        else if(itype == 1) {
          // B-AA
          if (jtype == 0 && ktype == 0) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_aa,g5_deriv_b_aa,g_b[3],
                  eta5_b_aa,lambda5_b_aa,zeta5_b_aa); 
          }
          // B-AB
          else if (jtype == 0 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_ab,g5_deriv_b_ab,g_b[4],
                  eta5_b_ab,lambda5_b_ab,zeta5_b_ab); 
          }
          // B-AC
          else if (jtype == 0 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_ac,g5_deriv_b_ac,g_b[5],
                  eta5_b_ac,lambda5_b_ac,zeta5_b_ac); 
          }
          // B-BB
          else if (jtype == 1 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_bb,g5_deriv_b_bb,g_b[6],
                  eta5_b_bb,lambda5_b_bb,zeta5_b_bb); 
          }
          // B-BC
          else if (jtype == 1 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_bc,g5_deriv_b_bc,g_b[7],
                  eta5_b_bc,lambda5_b_bc,zeta5_b_bc); 
          }
          // B-CC
          else if (jtype == 2 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_b_cc,g5_deriv_b_cc,g_b[8],
                  eta5_b_cc,lambda5_b_cc,zeta5_b_cc); 
          }
          else continue;
        }
        // C-
        else if(itype == 2) {
          // C-AA
          if (jtype == 0 && ktype == 0) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_aa,g5_deriv_c_aa,g_c[3],
                  eta5_c_aa,lambda5_c_aa,zeta5_c_aa); 
          }
          // C-AB
          else if (jtype == 0 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_ab,g5_deriv_c_ab,g_c[4],
                  eta5_c_ab,lambda5_c_ab,zeta5_c_ab); 
          }
          // C-AC
          else if (jtype == 0 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_ac,g5_deriv_c_ac,g_c[5],
                  eta5_c_ac,lambda5_c_ac,zeta5_c_ac); 
          }
          // C-BB
          else if (jtype == 1 && ktype == 1) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_bb,g5_deriv_c_bb,g_c[6],
                  eta5_c_bb,lambda5_c_bb,zeta5_c_bb); 
          }
          // C-BC
          else if (jtype == 1 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_bc,g5_deriv_c_bc,g_c[7],
                  eta5_c_bc,lambda5_c_bc,zeta5_c_bc); 
          }
          // C-CC
          else if (jtype == 2 && ktype == 2) {
            sf_g5(i,j,k,srsq1,srsq2,
                  x[i][0],x[i][1],x[i][2],
                  x[j][0],x[j][1],x[j][2],
                  x[k][0],x[k][1],x[k][2],
                  g5_c_cc,g5_deriv_c_cc,g_c[8],
                  eta5_c_cc,lambda5_c_cc,zeta5_c_cc); 
          }
          else continue;
        }

      } // for kk
    } // for jj

    // Energy i
    if (itype == 0) {
      hdnnp_e_a( g2_a_a,  g_a[0],
                 g2_a_b,  g_a[1],
                 g2_a_c,  g_a[2],
                 g5_a_aa, g_a[3],
                 g5_a_ab, g_a[4],
                 g5_a_ac, g_a[5],
                 g5_a_bb, g_a[6],
                 g5_a_bc, g_a[7],
                 g5_a_cc, g_a[8],
                 hidden_a, evdwl );
    }
    else if (itype == 1) {
      hdnnp_e_b( g2_b_a,  g_b[0],
                 g2_b_b,  g_b[1],
                 g2_b_c,  g_b[2],
                 g5_b_aa, g_b[3],
                 g5_b_ab, g_b[4],
                 g5_b_ac, g_b[5],
                 g5_b_bb, g_b[6],
                 g5_b_bc, g_b[7],
                 g5_b_cc, g_b[8],
                 hidden_b, evdwl );
    }
    else if (itype == 2) {
      hdnnp_e_c( g2_c_a,  g_c[0],
                 g2_c_b,  g_c[1],
                 g2_c_c,  g_c[2],
                 g5_c_aa, g_c[3],
                 g5_c_ab, g_c[4],
                 g5_c_ac, g_c[5],
                 g5_c_bb, g_c[6],
                 g5_c_bc, g_c[7],
                 g5_c_cc, g_c[8],
                 hidden_c, evdwl );
    }
    if (evflag) ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);

    // Force i
    if (itype == 0) {
      hdnnp_f_a( g2_deriv_a_a,  g_a[0],
                 g2_deriv_a_b,  g_a[1],
                 g2_deriv_a_c,  g_a[2],
                 g5_deriv_a_aa, g_a[3],
                 g5_deriv_a_ab, g_a[4],
                 g5_deriv_a_ac, g_a[5],
                 g5_deriv_a_bb, g_a[6],
                 g5_deriv_a_bc, g_a[7],
                 g5_deriv_a_cc, g_a[8],
                 i, hidden_a, fx, fy, fz );
    }
    else if (itype == 1) {
      hdnnp_f_b( g2_deriv_b_a,  g_b[0],
                 g2_deriv_b_b,  g_b[1],
                 g2_deriv_b_c,  g_b[2],
                 g5_deriv_b_aa, g_b[3],
                 g5_deriv_b_ab, g_b[4],
                 g5_deriv_b_ac, g_b[5],
                 g5_deriv_b_bb, g_b[6],
                 g5_deriv_b_bc, g_b[7],
                 g5_deriv_b_cc, g_b[8],
                 i, hidden_b, fx, fy, fz );
    }
    else if (itype == 2) {
      hdnnp_f_c( g2_deriv_c_a,  g_c[0],
                 g2_deriv_c_b,  g_c[1],
                 g2_deriv_c_c,  g_c[2],
                 g5_deriv_c_aa, g_c[3],
                 g5_deriv_c_ab, g_c[4],
                 g5_deriv_c_ac, g_c[5],
                 g5_deriv_c_bb, g_c[6],
                 g5_deriv_c_bc, g_c[7],
                 g5_deriv_c_cc, g_c[8],
                 i, hidden_c, fx, fy, fz );
    }
    f[i][0] += fx; f[i][1] += fy; f[i][2] += fz;

    // Force j
    for (jj = 0; jj < numshort; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      if (itype == 0) {
        hdnnp_f_a( g2_deriv_a_a,  g_a[0],
                   g2_deriv_a_b,  g_a[1],
                   g2_deriv_a_c,  g_a[2],
                   g5_deriv_a_aa, g_a[3],
                   g5_deriv_a_ab, g_a[4],
                   g5_deriv_a_ac, g_a[5],
                   g5_deriv_a_bb, g_a[6],
                   g5_deriv_a_bc, g_a[7],
                   g5_deriv_a_cc, g_a[8],
                   j, hidden_a, fx, fy, fz );
      }
      else if (itype == 1) {
        hdnnp_f_b( g2_deriv_b_a,  g_b[0],
                   g2_deriv_b_b,  g_b[1],
                   g2_deriv_b_c,  g_b[2],
                   g5_deriv_b_aa, g_b[3],
                   g5_deriv_b_ab, g_b[4],
                   g5_deriv_b_ac, g_b[5],
                   g5_deriv_b_bb, g_b[6],
                   g5_deriv_b_bc, g_b[7],
                   g5_deriv_b_cc, g_b[8],
                   j, hidden_b, fx, fy, fz );
      }
      else if (itype == 2) {
        hdnnp_f_c( g2_deriv_c_a,  g_c[0],
                   g2_deriv_c_b,  g_c[1],
                   g2_deriv_c_c,  g_c[2],
                   g5_deriv_c_aa, g_c[3],
                   g5_deriv_c_ab, g_c[4],
                   g5_deriv_c_ac, g_c[5],
                   g5_deriv_c_bb, g_c[6],
                   g5_deriv_c_bc, g_c[7],
                   g5_deriv_c_cc, g_c[8],
                   j, hidden_c, fx, fy, fz );
      }
      f[j][0] += fx; f[j][1] += fy; f[j][2] += fz;

    } // jj2

    // initialization
    for (i = 0; i < nall; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < g_a[0]; k++) g2_deriv_a_a[i][j][k] = 0.0;
      for (k = 0; k < g_a[1]; k++) g2_deriv_a_b[i][j][k] = 0.0;
      for (k = 0; k < g_a[2]; k++) g2_deriv_a_c[i][j][k] = 0.0;
      for (k = 0; k < g_a[3]; k++) g5_deriv_a_aa[i][j][k] = 0.0;
      for (k = 0; k < g_a[4]; k++) g5_deriv_a_ab[i][j][k] = 0.0;
      for (k = 0; k < g_a[5]; k++) g5_deriv_a_ac[i][j][k] = 0.0;
      for (k = 0; k < g_a[6]; k++) g5_deriv_a_bb[i][j][k] = 0.0;
      for (k = 0; k < g_a[7]; k++) g5_deriv_a_bc[i][j][k] = 0.0;
      for (k = 0; k < g_a[8]; k++) g5_deriv_a_cc[i][j][k] = 0.0;

      for (k = 0; k < g_b[0]; k++) g2_deriv_b_a[i][j][k] = 0.0;
      for (k = 0; k < g_b[1]; k++) g2_deriv_b_b[i][j][k] = 0.0;
      for (k = 0; k < g_b[2]; k++) g2_deriv_b_c[i][j][k] = 0.0;
      for (k = 0; k < g_b[3]; k++) g5_deriv_b_aa[i][j][k] = 0.0;
      for (k = 0; k < g_b[4]; k++) g5_deriv_b_ab[i][j][k] = 0.0;
      for (k = 0; k < g_b[5]; k++) g5_deriv_b_ac[i][j][k] = 0.0;
      for (k = 0; k < g_b[6]; k++) g5_deriv_b_bb[i][j][k] = 0.0;
      for (k = 0; k < g_b[7]; k++) g5_deriv_b_bc[i][j][k] = 0.0;
      for (k = 0; k < g_b[8]; k++) g5_deriv_b_cc[i][j][k] = 0.0;

      for (k = 0; k < g_c[0]; k++) g2_deriv_c_a[i][j][k] = 0.0;
      for (k = 0; k < g_c[1]; k++) g2_deriv_c_b[i][j][k] = 0.0;
      for (k = 0; k < g_c[2]; k++) g2_deriv_c_c[i][j][k] = 0.0;
      for (k = 0; k < g_c[3]; k++) g5_deriv_c_aa[i][j][k] = 0.0;
      for (k = 0; k < g_c[4]; k++) g5_deriv_c_ab[i][j][k] = 0.0;
      for (k = 0; k < g_c[5]; k++) g5_deriv_c_ac[i][j][k] = 0.0;
      for (k = 0; k < g_c[6]; k++) g5_deriv_c_bb[i][j][k] = 0.0;
      for (k = 0; k < g_c[7]; k++) g5_deriv_c_bc[i][j][k] = 0.0;
      for (k = 0; k < g_c[8]; k++) g5_deriv_c_cc[i][j][k] = 0.0;
    }}


  } // for ii


  if (vflag_fdotr) virial_fdotr_compute();


  memory->destroy(g2_deriv_a_a);
  memory->destroy(g2_deriv_a_b);
  memory->destroy(g2_deriv_a_c);
  memory->destroy(g5_deriv_a_aa);
  memory->destroy(g5_deriv_a_ab);
  memory->destroy(g5_deriv_a_ac);
  memory->destroy(g5_deriv_a_bb);
  memory->destroy(g5_deriv_a_bc);
  memory->destroy(g5_deriv_a_cc);

  memory->destroy(g2_deriv_b_a);
  memory->destroy(g2_deriv_b_b);
  memory->destroy(g2_deriv_b_c);
  memory->destroy(g5_deriv_b_aa);
  memory->destroy(g5_deriv_b_ab);
  memory->destroy(g5_deriv_b_ac);
  memory->destroy(g5_deriv_b_bb);
  memory->destroy(g5_deriv_b_bc);
  memory->destroy(g5_deriv_b_cc);

  memory->destroy(g2_deriv_c_a);
  memory->destroy(g2_deriv_c_b);
  memory->destroy(g2_deriv_c_c);
  memory->destroy(g5_deriv_c_aa);
  memory->destroy(g5_deriv_c_ab);
  memory->destroy(g5_deriv_c_ac);
  memory->destroy(g5_deriv_c_bb);
  memory->destroy(g5_deriv_c_bc);
  memory->destroy(g5_deriv_c_cc);


}

/* ---------------------------------------------------------------------- */

void PairNNP3v::allocate()
{

  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort,maxshort,"pair:neighshort");
  map = new int[n+1];

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNNP3v::settings(int narg, char **/*arg*/)
{

  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNNP3v::coeff(int narg, char **arg)
{

  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();


  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNNP3v::init_style()
{

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style NNP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style NNP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  arraynnp();

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNNP3v::init_one(int i, int j)
{

  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");


  return cutmax;

}

/* ---------------------------------------------------------------------- */

void PairNNP3v::read_file(char *file)
{
  //Read file decleared in pair_coeff tag

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open NNP potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // Find weight parameters from file
  //rewind(fp);
  char sen[MAXLINE],*mojiretsu;
  int ns,nsuji,dummyi;
  int end_of_file = 0;
  int end_of_weight_a = 0;
  int end_of_weight_b = 0;
  int end_of_weight_c = 0;

  while (1) {
    if (comm->me == 0) {
      mojiretsu = fgets(sen,MAXLINE,fp);
      if(mojiretsu == NULL) {
        end_of_file = 1;
        fclose(fp);
      } else ns = strlen(sen) + 1;
    }
    MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
    if (end_of_file) break;
    MPI_Bcast(&ns,1,MPI_INT,0,world);
    MPI_Bcast(sen,ns,MPI_CHAR,0,world);

    if ((mojiretsu = strchr(sen,'#'))) *mojiretsu = '\0';
    nsuji = atom->count_words(sen);
    // nwords = # of words per line
    if (nsuji == 0) continue;

    //char sen2[MAXLINE],*moji;
    //int nlayer_a,dummyi;
    //int end_of_file = 0;
    //int end_of_weight = 0;

    mojiretsu = strtok(sen," \t\n\r\f");
    while (mojiretsu) {
      if (strcmp(mojiretsu,"cutoff_distance") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if(mojiretsu == NULL) return;
        rc = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g_a") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        dummyi = atof(mojiretsu);
        g_a = new int[ dummyi ];
      }
      else if (strcmp(mojiretsu,"num_g2_a_a") == 0) { 
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[0] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_a_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_a_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[2] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_aa") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[3] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_ab") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[4] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_ac") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[5] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_bb") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[6] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_bc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[7] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_a_cc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_a[8] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        dummyi = atof(mojiretsu);
        g_b = new int[ dummyi ];
      }
      else if (strcmp(mojiretsu,"num_g2_b_a") == 0) { 
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[0] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_b_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_b_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[2] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_aa") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[3] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_ab") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[4] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_ac") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[5] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_bb") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[6] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_bc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[7] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_b_cc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_b[8] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        dummyi = atof(mojiretsu);
        g_c = new int[ dummyi ];
      }
      else if (strcmp(mojiretsu,"num_g2_c_a") == 0) { 
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[0] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_c_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g2_c_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[2] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_aa") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[3] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_ab") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[4] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_ac") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[5] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_bb") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[6] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_bc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[7] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"num_g5_c_cc") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        g_c[8] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"nlayer_a") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nlayer_a = atof(mojiretsu);
        layer_a = new int[nlayer_a];
        layer_a[0] =   g_a[0]  + g_a[1]  + g_a[2]  
                     + g_a[3]  + g_a[4]  + g_a[5]  
                     + g_a[6]  + g_a[7]
                     + g_a[8];
      }
      else if (strcmp(mojiretsu,"node_h1_a") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_a[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"node_h2_a") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_a[2] = atof(mojiretsu);
      } // node_h3_a, node_h4_a, need to revise!!
      else if (strcmp(mojiretsu,"node_out_a") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_a[nlayer_a-1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"nlayer_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nlayer_b = atof(mojiretsu);
        layer_b = new int[nlayer_b];
        layer_b[0] =   g_b[0] + g_b[1] + g_b[2]
                     + g_b[3] + g_b[4] + g_b[5]
                     + g_b[6] + g_b[7]
                     + g_b[8];
      }
      else if (strcmp(mojiretsu,"node_h1_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_b[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"node_h2_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_b[2] = atof(mojiretsu);
      } // node_h3_a, node_h4_a, need to revise!!
      else if (strcmp(mojiretsu,"node_out_b") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_b[nlayer_b-1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"nlayer_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nlayer_c = atof(mojiretsu);
        layer_c = new int[nlayer_c];
        layer_c[0] =   g_c[0] + g_c[1] + g_c[2]
                     + g_c[3] + g_c[4] + g_c[5]
                     + g_c[6] + g_c[7]
                     + g_c[8];
      }
      else if (strcmp(mojiretsu,"node_h1_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_c[1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"node_h2_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_c[2] = atof(mojiretsu);
      } // node_h3_a, node_h4_a, need to revise!!
      else if (strcmp(mojiretsu,"node_out_c") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        layer_c[nlayer_c-1] = atof(mojiretsu);
      }
      else if (strcmp(mojiretsu,"param_g2_a_a") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_a[0]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_a[0]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_a_a = new double[ g_a[0] ];
        rs2_a_a  = new double[ g_a[0] ];
        for (int i=0; i<g_a[0]; i++) {
          eta2_a_a[i] = atof( suji[2*i] );
          rs2_a_a[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_a_b") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_a[1]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_a[1]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_a_b = new double[ g_a[1] ];
        rs2_a_b  = new double[ g_a[1] ];
        for (int i=0; i<g_a[1]; i++) {
          eta2_a_b[i] = atof( suji[2*i] );
          rs2_a_b[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_a_c") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_a[2]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_a[2]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_a_c = new double[ g_a[2] ];
        rs2_a_c  = new double[ g_a[2] ];
        for (int i=0; i<g_a[2]; i++) {
          eta2_a_c[i] = atof( suji[2*i] );
          rs2_a_c[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_aa") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[3]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break; 
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[3]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_aa    = new double[ g_a[3] ];
        theta5_a_aa  = new double[ g_a[3] ];
        zeta5_a_aa   = new double[ g_a[3] ];
        lambda5_a_aa = new double[ g_a[3] ];
        for (int i=0; i<g_a[3]; i++) {
          eta5_a_aa[i]    = atof( suji[4*i] );
          theta5_a_aa[i]  = atof( suji[4*i+1] );
          zeta5_a_aa[i]   = atof( suji[4*i+2] );
          lambda5_a_aa[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_ab") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[4]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[4]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_ab    = new double[ g_a[4] ];
        theta5_a_ab  = new double[ g_a[4] ];
        zeta5_a_ab   = new double[ g_a[4] ];
        lambda5_a_ab = new double[ g_a[4] ];
        for (int i=0; i<g_a[4]; i++) {
          eta5_a_ab[i]    = atof( suji[4*i] );
          theta5_a_ab[i]  = atof( suji[4*i+1] );
          zeta5_a_ab[i]   = atof( suji[4*i+2] );
          lambda5_a_ab[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_ac") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[5]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[5]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_ac    = new double[ g_a[5] ];
        theta5_a_ac  = new double[ g_a[5] ];
        zeta5_a_ac   = new double[ g_a[5] ];
        lambda5_a_ac = new double[ g_a[5] ];
        for (int i=0; i<g_a[5]; i++) {
          eta5_a_ac[i]    = atof( suji[4*i] );
          theta5_a_ac[i]  = atof( suji[4*i+1] );
          zeta5_a_ac[i]   = atof( suji[4*i+2] );
          lambda5_a_ac[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_bb") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[6]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[6]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_bb    = new double[ g_a[6] ];
        theta5_a_bb  = new double[ g_a[6] ];
        zeta5_a_bb   = new double[ g_a[6] ];
        lambda5_a_bb = new double[ g_a[6] ];
        for (int i=0; i<g_a[6]; i++) {
          eta5_a_bb[i]    = atof( suji[4*i] );
          theta5_a_bb[i]  = atof( suji[4*i+1] );
          zeta5_a_bb[i]   = atof( suji[4*i+2] );
          lambda5_a_bb[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_bc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[7]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[7]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_bc    = new double[ g_a[7] ];
        theta5_a_bc  = new double[ g_a[7] ];
        zeta5_a_bc   = new double[ g_a[7] ];
        lambda5_a_bc = new double[ g_a[7] ];
        for (int i=0; i<g_a[7]; i++) {
          eta5_a_bc[i]    = atof( suji[4*i] );
          theta5_a_bc[i]  = atof( suji[4*i+1] );
          zeta5_a_bc[i]   = atof( suji[4*i+2] );
          lambda5_a_bc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_a_cc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_a[8]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_a[8]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_a_cc    = new double[ g_a[8] ];
        theta5_a_cc  = new double[ g_a[8] ];
        zeta5_a_cc   = new double[ g_a[8] ];
        lambda5_a_cc = new double[ g_a[8] ];
        for (int i=0; i<g_a[8]; i++) {
          eta5_a_cc[i]    = atof( suji[4*i] );
          theta5_a_cc[i]  = atof( suji[4*i+1] );
          zeta5_a_cc[i]   = atof( suji[4*i+2] );
          lambda5_a_cc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_b_a") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_b[0]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_b[0]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_b_a = new double[ g_b[0] ];
        rs2_b_a  = new double[ g_b[0] ];
        for (int i=0; i<g_b[0]; i++) {
          eta2_b_a[i] = atof( suji[2*i] );
          rs2_b_a[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_b_b") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_b[1]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_b[1]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_b_b = new double[ g_b[1] ];
        rs2_b_b  = new double[ g_b[1] ];
        for (int i=0; i<g_b[1]; i++) {
          eta2_b_b[i] = atof( suji[2*i] );
          rs2_b_b[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_b_c") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_b[2]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_b[2]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_b_c = new double[ g_b[2] ];
        rs2_b_c  = new double[ g_b[2] ];
        for (int i=0; i<g_b[2]; i++) {
          eta2_b_c[i] = atof( suji[2*i] );
          rs2_b_c[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_aa") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[3]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break; 
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[3]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_aa    = new double[ g_b[3] ];
        theta5_b_aa  = new double[ g_b[3] ];
        zeta5_b_aa   = new double[ g_b[3] ];
        lambda5_b_aa = new double[ g_b[3] ];
        for (int i=0; i<g_b[3]; i++) {
          eta5_b_aa[i]    = atof( suji[4*i] );
          theta5_b_aa[i]  = atof( suji[4*i+1] );
          zeta5_b_aa[i]   = atof( suji[4*i+2] );
          lambda5_b_aa[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_ab") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[4]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[4]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_ab    = new double[ g_b[4] ];
        theta5_b_ab  = new double[ g_b[4] ];
        zeta5_b_ab   = new double[ g_b[4] ];
        lambda5_b_ab = new double[ g_b[4] ];
        for (int i=0; i<g_b[4]; i++) {
          eta5_b_ab[i]    = atof( suji[4*i] );
          theta5_b_ab[i]  = atof( suji[4*i+1] );
          zeta5_b_ab[i]   = atof( suji[4*i+2] );
          lambda5_b_ab[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_ac") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[5]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[5]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_ac    = new double[ g_b[5] ];
        theta5_b_ac  = new double[ g_b[5] ];
        zeta5_b_ac   = new double[ g_b[5] ];
        lambda5_b_ac = new double[ g_b[5] ];
        for (int i=0; i<g_b[5]; i++) {
          eta5_b_ac[i]    = atof( suji[4*i] );
          theta5_b_ac[i]  = atof( suji[4*i+1] );
          zeta5_b_ac[i]   = atof( suji[4*i+2] );
          lambda5_b_ac[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_bb") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[6]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[6]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_bb    = new double[ g_b[6] ];
        theta5_b_bb  = new double[ g_b[6] ];
        zeta5_b_bb   = new double[ g_b[6] ];
        lambda5_b_bb = new double[ g_b[6] ];
        for (int i=0; i<g_b[6]; i++) {
          eta5_b_bb[i]    = atof( suji[4*i] );
          theta5_b_bb[i]  = atof( suji[4*i+1] );
          zeta5_b_bb[i]   = atof( suji[4*i+2] );
          lambda5_b_bb[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_bc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[7]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[7]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_bc    = new double[ g_b[7] ];
        theta5_b_bc  = new double[ g_b[7] ];
        zeta5_b_bc   = new double[ g_b[7] ];
        lambda5_b_bc = new double[ g_b[7] ];
        for (int i=0; i<g_b[7]; i++) {
          eta5_b_bc[i]    = atof( suji[4*i] );
          theta5_b_bc[i]  = atof( suji[4*i+1] );
          zeta5_b_bc[i]   = atof( suji[4*i+2] );
          lambda5_b_bc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_b_cc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_b[8]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_b[8]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_b_cc    = new double[ g_b[8] ];
        theta5_b_cc  = new double[ g_b[8] ];
        zeta5_b_cc   = new double[ g_b[8] ];
        lambda5_b_cc = new double[ g_b[8] ];
        for (int i=0; i<g_b[8]; i++) {
          eta5_b_cc[i]    = atof( suji[4*i] );
          theta5_b_cc[i]  = atof( suji[4*i+1] );
          zeta5_b_cc[i]   = atof( suji[4*i+2] );
          lambda5_b_cc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_c_a") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_c[0]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_c[0]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_c_a = new double[ g_c[0] ];
        rs2_c_a  = new double[ g_c[0] ];
        for (int i=0; i<g_c[0]; i++) {
          eta2_c_a[i] = atof( suji[2*i] );
          rs2_c_a[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_c_b") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_c[1]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_c[1]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_c_b = new double[ g_c[1] ];
        rs2_c_b  = new double[ g_c[1] ];
        for (int i=0; i<g_c[1]; i++) {
          eta2_c_b[i] = atof( suji[2*i] );
          rs2_c_b[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g2_c_c") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2*g_c[2]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2*g_c[2]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta2_c_c = new double[ g_c[2] ];
        rs2_c_c  = new double[ g_c[2] ];
        for (int i=0; i<g_c[2]; i++) {
          eta2_c_c[i] = atof( suji[2*i] );
          rs2_c_c[i] = atof( suji[2*i+1] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_aa") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[3]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break; 
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[3]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_aa    = new double[ g_c[3] ];
        theta5_c_aa  = new double[ g_c[3] ];
        zeta5_c_aa   = new double[ g_c[3] ];
        lambda5_c_aa = new double[ g_c[3] ];
        for (int i=0; i<g_c[3]; i++) {
          eta5_c_aa[i]    = atof( suji[4*i] );
          theta5_c_aa[i]  = atof( suji[4*i+1] );
          zeta5_c_aa[i]   = atof( suji[4*i+2] );
          lambda5_c_aa[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_ab") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[4]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[4]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_ab    = new double[ g_c[4] ];
        theta5_c_ab  = new double[ g_c[4] ];
        zeta5_c_ab   = new double[ g_c[4] ];
        lambda5_c_ab = new double[ g_c[4] ];
        for (int i=0; i<g_c[4]; i++) {
          eta5_c_ab[i]    = atof( suji[4*i] );
          theta5_c_ab[i]  = atof( suji[4*i+1] );
          zeta5_c_ab[i]   = atof( suji[4*i+2] );
          lambda5_c_ab[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_ac") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[5]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[5]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_ac    = new double[ g_c[5] ];
        theta5_c_ac  = new double[ g_c[5] ];
        zeta5_c_ac   = new double[ g_c[5] ];
        lambda5_c_ac = new double[ g_c[5] ];
        for (int i=0; i<g_c[5]; i++) {
          eta5_c_ac[i]    = atof( suji[4*i] );
          theta5_c_ac[i]  = atof( suji[4*i+1] );
          zeta5_c_ac[i]   = atof( suji[4*i+2] );
          lambda5_c_ac[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_bb") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[6]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[6]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_bb    = new double[ g_c[6] ];
        theta5_c_bb  = new double[ g_c[6] ];
        zeta5_c_bb   = new double[ g_c[6] ];
        lambda5_c_bb = new double[ g_c[6] ];
        for (int i=0; i<g_c[6]; i++) {
          eta5_c_bb[i]    = atof( suji[4*i] );
          theta5_c_bb[i]  = atof( suji[4*i+1] );
          zeta5_c_bb[i]   = atof( suji[4*i+2] );
          lambda5_c_bb[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_bc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[7]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[7]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_bc    = new double[ g_c[7] ];
        theta5_c_bc  = new double[ g_c[7] ];
        zeta5_c_bc   = new double[ g_c[7] ];
        lambda5_c_bc = new double[ g_c[7] ];
        for (int i=0; i<g_c[7]; i++) {
          eta5_c_bc[i]    = atof( suji[4*i] );
          theta5_c_bc[i]  = atof( suji[4*i+1] );
          zeta5_c_bc[i]   = atof( suji[4*i+2] );
          lambda5_c_bc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"param_g5_c_cc") == 0) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 4*g_c[8]) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 4*g_c[8]+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        eta5_c_cc    = new double[ g_c[8] ];
        theta5_c_cc  = new double[ g_c[8] ];
        zeta5_c_cc   = new double[ g_c[8] ];
        lambda5_c_cc = new double[ g_c[8] ];
        for (int i=0; i<g_c[8]; i++) {
          eta5_c_cc[i]    = atof( suji[4*i] );
          theta5_c_cc[i]  = atof( suji[4*i+1] );
          zeta5_c_cc[i]   = atof( suji[4*i+2] );
          lambda5_c_cc[i] = atof( suji[4*i+3] );
        }
        delete [] suji,sen2,moji;
      }
      else if (strcmp(mojiretsu,"g2_a_a_max") == 0) {
        g2_a_a_max = new double[ g_a[0] ];
        for (int i = 0; i < g_a[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_a_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_a_a_min") == 0) {
        g2_a_a_min = new double[ g_a[0] ];
        for (int i = 0; i < g_a[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_a_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_a_b_max") == 0) {
        g2_a_b_max = new double[ g_a[1] ];
        for (int i = 0; i < g_a[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_b_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_a_b_min") == 0) {
        g2_a_b_min = new double[ g_a[1] ];
        for (int i = 0; i < g_a[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_b_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_a_c_max") == 0) {
        g2_a_c_max = new double[ g_a[2] ];
        for (int i = 0; i < g_a[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_c_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_a_c_min") == 0) {
        g2_a_c_min = new double[ g_a[2] ];
        for (int i = 0; i < g_a[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_a_c_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_aa_max") == 0) {
        g5_a_aa_max = new double[ g_a[3] ];
        for (int i = 0; i < g_a[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_aa_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_aa_min") == 0) {
        g5_a_aa_min = new double[ g_a[3] ];
        for (int i = 0; i < g_a[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_aa_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_ab_max") == 0) {
        g5_a_ab_max = new double[ g_a[4] ];
        for (int i = 0; i < g_a[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_ab_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_ab_min") == 0) {
        g5_a_ab_min = new double[ g_a[4] ];
        for (int i = 0; i < g_a[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_ab_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_ac_max") == 0) {
        g5_a_ac_max = new double[ g_a[5] ];
        for (int i = 0; i < g_a[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_ac_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_ac_min") == 0) {
        g5_a_ac_min = new double[ g_a[5] ];
        for (int i = 0; i < g_a[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_ac_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_bb_max") == 0) {
        g5_a_bb_max = new double[ g_a[6] ];
        for (int i = 0; i < g_a[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_bb_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_bb_min") == 0) {
        g5_a_bb_min = new double[ g_a[6] ];
        for (int i = 0; i < g_a[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_bb_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_bc_max") == 0) {
        g5_a_bc_max = new double[ g_a[7] ];
        for (int i = 0; i < g_a[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_bc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_bc_min") == 0) {
        g5_a_bc_min = new double[ g_a[7] ];
        for (int i = 0; i < g_a[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_bc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_cc_max") == 0) {
        g5_a_cc_max = new double[ g_a[8] ];
        for (int i = 0; i < g_a[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_cc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_a_cc_min") == 0) {
        g5_a_cc_min = new double[ g_a[8] ];
        for (int i = 0; i < g_a[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_a_cc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_a_max") == 0) {
        g2_b_a_max = new double[ g_b[0] ];
        for (int i = 0; i < g_b[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_a_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_a_min") == 0) {
        g2_b_a_min = new double[ g_b[0] ];
        for (int i = 0; i < g_b[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_a_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_b_max") == 0) {
        g2_b_b_max = new double[ g_b[1] ];
        for (int i = 0; i < g_b[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_b_max[dummyi] = atof( suji[1] );
//printf("\n g2_b_b_max %15.9f",g2_b_b_max[dummyi]);
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_b_min") == 0) {
        g2_b_b_min = new double[ g_b[1] ];
        for (int i = 0; i < g_b[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_b_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_c_max") == 0) {
        g2_b_c_max = new double[ g_b[2] ];
        for (int i = 0; i < g_b[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_c_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_b_c_min") == 0) {
        g2_b_c_min = new double[ g_b[2] ];
        for (int i = 0; i < g_b[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_b_c_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_aa_max") == 0) {
        g5_b_aa_max = new double[ g_b[3] ];
        for (int i = 0; i < g_b[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_aa_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_aa_min") == 0) {
        g5_b_aa_min = new double[ g_b[3] ];
        for (int i = 0; i < g_b[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_aa_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_ab_max") == 0) {
        g5_b_ab_max = new double[ g_b[4] ];
        for (int i = 0; i < g_b[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_ab_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_ab_min") == 0) {
        g5_b_ab_min = new double[ g_b[4] ];
        for (int i = 0; i < g_b[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_ab_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_ac_max") == 0) {
        g5_b_ac_max = new double[ g_b[5] ];
        for (int i = 0; i < g_b[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_ac_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_ac_min") == 0) {
        g5_b_ac_min = new double[ g_b[5] ];
        for (int i = 0; i < g_b[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_ac_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_bb_max") == 0) {
        g5_b_bb_max = new double[ g_b[6] ];
        for (int i = 0; i < g_b[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_bb_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_bb_min") == 0) {
        g5_b_bb_min = new double[ g_b[6] ];
        for (int i = 0; i < g_b[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_bb_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_bc_max") == 0) {
        g5_b_bc_max = new double[ g_b[7] ];
        for (int i = 0; i < g_b[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_bc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_bc_min") == 0) {
        g5_b_bc_min = new double[ g_b[7] ];
        for (int i = 0; i < g_b[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_bc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_cc_max") == 0) {
        g5_b_cc_max = new double[ g_b[8] ];
        for (int i = 0; i < g_b[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_cc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_b_cc_min") == 0) {
        g5_b_cc_min = new double[ g_b[8] ];
        for (int i = 0; i < g_b[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_b_cc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_a_max") == 0) {
        g2_c_a_max = new double[ g_c[0] ];
        for (int i = 0; i < g_c[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_a_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_a_min") == 0) {
        g2_c_a_min = new double[ g_c[0] ];
        for (int i = 0; i < g_c[0]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_a_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_b_max") == 0) {
        g2_c_b_max = new double[ g_c[1] ];
        for (int i = 0; i < g_c[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_b_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_b_min") == 0) {
        g2_c_b_min = new double[ g_c[1] ];
        for (int i = 0; i < g_c[1]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_b_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_c_max") == 0) {
        g2_c_c_max = new double[ g_c[2] ];
        for (int i = 0; i < g_c[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_c_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g2_c_c_min") == 0) {
        g2_c_c_min = new double[ g_c[2] ];
        for (int i = 0; i < g_c[2]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
        dummyi = atof( suji[0] );
        g2_c_c_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_aa_max") == 0) {
        g5_c_aa_max = new double[ g_c[3] ];
        for (int i = 0; i < g_c[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_aa_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_aa_min") == 0) {
        g5_c_aa_min = new double[ g_c[3] ];
        for (int i = 0; i < g_c[3]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_aa_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_ab_max") == 0) {
        g5_c_ab_max = new double[ g_c[4] ];
        for (int i = 0; i < g_c[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_ab_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_ab_min") == 0) {
        g5_c_ab_min = new double[ g_c[4] ];
        for (int i = 0; i < g_c[4]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_ab_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_ac_max") == 0) {
        g5_c_ac_max = new double[ g_c[5] ];
        for (int i = 0; i < g_c[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_ac_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_ac_min") == 0) {
        g5_c_ac_min = new double[ g_c[5] ];
        for (int i = 0; i < g_c[5]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_ac_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_bb_max") == 0) {
        g5_c_bb_max = new double[ g_c[6] ];
        for (int i = 0; i < g_c[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_bb_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_bb_min") == 0) {
        g5_c_bb_min = new double[ g_c[6] ];
        for (int i = 0; i < g_c[6]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_bb_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_bc_max") == 0) {
        g5_c_bc_max = new double[ g_c[7] ];
        for (int i = 0; i < g_c[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_bc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_bc_min") == 0) {
        g5_c_bc_min = new double[ g_c[7] ];
        for (int i = 0; i < g_c[7]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_bc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_cc_max") == 0) {
        g5_c_cc_max = new double[ g_c[8] ];
        for (int i = 0; i < g_c[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1; 
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_cc_max[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"g5_c_cc_min") == 0) {
        g5_c_cc_min = new double[ g_c[8] ];
        for (int i = 0; i < g_c[8]; i++) {
        char sen2[MAXLINE],*moji;
        if (comm->me == 0) {
          moji = fgets(sen2,MAXLINE,fp);
          if (moji == NULL) {
            end_of_file = 1;
            fclose(fp);
          } else ns = strlen(sen2) + 1;
        }
        MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
        if (end_of_file) break;
        MPI_Bcast(&ns,1,MPI_INT,0,world);
        MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
        if ((moji = strchr(sen2,'#'))) *moji = '\0';
        nsuji = atom->count_words(sen2);
        if (nsuji == 0) continue;
        while (nsuji < 2) {
          ns = strlen(sen2);
          if (comm->me == 0) {
            moji = fgets(&sen2[ns],MAXLINE-ns,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
        }
        char **suji = new char*[ 2+1 ];
        nsuji = 0;
        suji[nsuji++] = strtok(sen2," \t\n\r\f");
        while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          g5_c_cc_min[dummyi] = atof( suji[1] );
        delete [] suji,sen2,moji;
        }
      }
      else if (strcmp(mojiretsu,"weight_a_params") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nweight_a = atof(mojiretsu);
        weight_a = new double[ nweight_a ];
      }
      else if (strcmp(mojiretsu,"weight_a") == 0) {
        for (int i = 0; i < nweight_a; i++) {
          char sen2[MAXLINE],*moji;
          if (comm->me == 0) {
            moji = fgets(sen2,MAXLINE,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
          if (nsuji == 0) continue;
          while (nsuji < 2) {
            ns = strlen(sen2);
            if (comm->me == 0) {
              moji = fgets(&sen2[ns],MAXLINE-ns,fp);
              if (moji == NULL) {
                end_of_file = 1;
                fclose(fp);
              } else ns = strlen(sen2) + 1;
            }
            MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
            if (end_of_file) break;
            MPI_Bcast(&ns,1,MPI_INT,0,world);
            MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
            if ((moji = strchr(sen2,'#'))) *moji = '\0';
            nsuji = atom->count_words(sen2);
          }
          char **suji = new char*[ 2+1 ];
          nsuji = 0;
          suji[nsuji++] = strtok(sen2," \t\n\r\f");
          while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          weight_a[dummyi] = atof( suji[1] );
          delete [] suji,sen2,moji;
          if (dummyi == nweight_a) {
            end_of_weight_a = 1;
            if (end_of_weight_a) break;
          } // for i
        }
      }
      else if (strcmp(mojiretsu,"weight_b_params") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nweight_b = atof(mojiretsu);
        weight_b = new double[ nweight_b ];
      }
      else if (strcmp(mojiretsu,"weight_b") == 0) {
        for (int i = 0; i < nweight_b; i++) {
          char sen2[MAXLINE],*moji;
          if (comm->me == 0) {
            moji = fgets(sen2,MAXLINE,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
          if (nsuji == 0) continue;
          while (nsuji < 2) {
            ns = strlen(sen2);
            if (comm->me == 0) {
              moji = fgets(&sen2[ns],MAXLINE-ns,fp);
              if (moji == NULL) {
                end_of_file = 1;
                fclose(fp);
              } else ns = strlen(sen2) + 1;
            }
            MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
            if (end_of_file) break;
            MPI_Bcast(&ns,1,MPI_INT,0,world);
            MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
            if ((moji = strchr(sen2,'#'))) *moji = '\0';
            nsuji = atom->count_words(sen2);
          }
          char **suji = new char*[ 2+1 ];
          nsuji = 0;
          suji[nsuji++] = strtok(sen2," \t\n\r\f");
          while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          weight_b[dummyi] = atof( suji[1] );
          delete [] suji,sen2,moji;
          if (dummyi == nweight_b) {
            end_of_weight_b = 1;
            if (end_of_weight_b) break;
          } // for i
        }
      }
      else if (strcmp(mojiretsu,"weight_c_params") == 0) {
        mojiretsu = strtok(NULL," \t\n\r\f");
        if (mojiretsu == NULL) return;
        nweight_c = atof(mojiretsu);
        weight_c = new double[ nweight_c ];
      }
      else if (strcmp(mojiretsu,"weight_c") == 0) {
        for (int i = 0; i < nweight_c; i++) {
          char sen2[MAXLINE],*moji;
          if (comm->me == 0) {
            moji = fgets(sen2,MAXLINE,fp);
            if (moji == NULL) {
              end_of_file = 1;
              fclose(fp);
            } else ns = strlen(sen2) + 1;
          }
          MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
          if (end_of_file) break;
          MPI_Bcast(&ns,1,MPI_INT,0,world);
          MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
          if ((moji = strchr(sen2,'#'))) *moji = '\0';
          nsuji = atom->count_words(sen2);
          if (nsuji == 0) continue;
          while (nsuji < 2) {
            ns = strlen(sen2);
            if (comm->me == 0) {
              moji = fgets(&sen2[ns],MAXLINE-ns,fp);
              if (moji == NULL) {
                end_of_file = 1;
                fclose(fp);
              } else ns = strlen(sen2) + 1;
            }
            MPI_Bcast(&end_of_file,1,MPI_INT,0,world);
            if (end_of_file) break;
            MPI_Bcast(&ns,1,MPI_INT,0,world);
            MPI_Bcast(sen2,ns,MPI_CHAR,0,world);
            if ((moji = strchr(sen2,'#'))) *moji = '\0';
            nsuji = atom->count_words(sen2);
          }
          char **suji = new char*[ 2+1 ];
          nsuji = 0;
          suji[nsuji++] = strtok(sen2," \t\n\r\f");
          while ((suji[nsuji++] = strtok(NULL," \t\n\r\f"))) continue;
          dummyi = atof( suji[0] );
          weight_c[dummyi] = atof( suji[1] );
          delete [] suji,sen2,moji;
          if (dummyi == nweight_c) {
            end_of_weight_c = 1;
            if (end_of_weight_c) break;
          } // for i
        }
      }
      mojiretsu = strtok(NULL," \t\n\r\f");
    } // while (mojiretsu)
    if(end_of_weight_c) break;
  } // while (1)


  // count nhidden
  nhidden_a = 0;
  for (int i = 1; i < nlayer_a-1; i++) {
    nhidden_a += layer_a[i] + 1;
  }
  nhidden_b = 0;
  for (int i = 1; i < nlayer_b-1; i++) {
    nhidden_b += layer_b[i] + 1;
  }
  nhidden_c = 0;
  for (int i = 1; i < nlayer_c-1; i++) {
    nhidden_c += layer_c[i] + 1;
  }

}

/* ---------------------------------------------------------------------- */

void PairNNP3v::setup_params()
{

  int i,j,k,m,n;

  // cutoff distance
  cutmax = rc;

}

/* ---------------------------------------------------------------------- */

void PairNNP3v::sf_g2(int i, int j, double rij,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double *(&g2), double ***(&g2_deriv), int num_g2, 
                     double *(&eta2), double *(&rs2) )
{
  double g2_const_1,g2_const_2,g2_const_x,g2_const_y,g2_const_z,g2_const_d;
  double fc1,fc2,sin1,inv_rij;

  fc1 = fc(rij);
  fc2 = 2.0 * fc1;
  inv_rij = 1.0 / rij;
  sin1 = 0.50 * (MY_PI/rc) * sin( (MY_PI*rij/rc) );

  for (int ig2 = 0; ig2 < num_g2; ig2++) {

    g2_const_1 = exp( -eta2[ig2] * pow(rij - rs2[ig2], 2.0) );
    g2[ig2] += g2_const_1 * fc1;
    //g2[ig2] += g2_const_1 * fc(rij);

    g2_const_d = -( eta2[ig2] * (rij - rs2[ig2]) * fc2 ) - sin1;
    //g2_const_d = - ( 2.0 * eta2[ig2] * (rij - rs2[ig2]) * fc(rij) )
    //             - ( 0.50 * (MY_PI/rc) * sin( (MY_PI*rij/rc) ) );

    g2_const_2 = g2_const_d * g2_const_1 * inv_rij;
    //g2_const_2 = g2_const_d * g2_const_1 / rij;

    g2_const_x = g2_const_2 * (x1 - x2);
    g2_const_y = g2_const_2 * (y1 - y2);
    g2_const_z = g2_const_2 * (z1 - z2);

//  (i,i)
    g2_deriv[i][0][ig2] += g2_const_x; //g2_const_1 * ((x1 - x2)/rij) * g2_const_d;
    g2_deriv[i][1][ig2] += g2_const_y; //g2_const_1 * ((y1 - y2)/rij) * g2_const_d;
    g2_deriv[i][2][ig2] += g2_const_z; //g2_const_1 * ((z1 - z2)/rij) * g2_const_d;

//  (i,j)
    g2_deriv[j][0][ig2] -= g2_const_x; //g2_const_1 * ((x1 - x2)/rij) * g2_const_d;
    g2_deriv[j][1][ig2] -= g2_const_y; //g2_const_1 * ((y1 - y2)/rij) * g2_const_d;
    g2_deriv[j][2][ig2] -= g2_const_z; //2_const_1 * ((z1 - z2)/rij) * g2_const_d;

  } // for ig2


}

/* ---------------------------------------------------------------------- */

double PairNNP3v::fc(double rij)
{
  return 0.5 * ( cos( ( MY_PI * rij )/rc ) + 1.0 );
}

/* ---------------------------------------------------------------------- */

double PairNNP3v::fc_deriv(double rij, double xi, double xj)
{
  return -( 0.5 * MY_PI * ( xi - xj ) * sin( ( MY_PI * rij )/rc ) )
          /( rc * rij ); 
}

/* ---------------------------------------------------------------------- */
void PairNNP3v::sf_g5(int i, int j, int k, double rij, double rik,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double x3, double y3, double z3,
                     double *(&g5), double ***(&g5_deriv), int num_g5,
                     double *(&eta5), double *(&lambda5), double *(&zeta5) )
{
  double rij_v[3], rik_v[3];
  double theta_ijk;
  double theta_ijk_dx,theta_ijk_dy,theta_ijk_dz;
  double theta_ijk_dx1,theta_ijk_dy1,theta_ijk_dz1;
  double theta_ijk_dx2,theta_ijk_dy2,theta_ijk_dz2;
  double g5_const_x,g5_const_y,g5_const_z;
  double fc_ij_ik,exp_ij_ik;
  double g5_const_1, g5_const_2;
  double rij_rik,rijrik,inv_rijrik,rij_p2,rik_p2,rij_p3,rik_p3,rikrij_p3,rijrik_p3,inv_rijrik_p3,inv_rikrij_p3;
  double theta_const_1,theta_const_2;
  double g5_const_3,g5_const_4,g5_const_5,g5_const_6,g5_const_7,g5_const_8;
  double g5_const_x0,g5_const_y0,g5_const_z0;
  double g5_const_x1,g5_const_y1,g5_const_z1;
  double g5_const_x2,g5_const_y2,g5_const_z2;
  double fc1,fc2,fc_dx12,fc_dy12,fc_dz12,fc_dx13,fc_dy13,fc_dz13;

  rij_v[0] = x2 - x1;
  rij_v[1] = y2 - y1;
  rij_v[2] = z2 - z1;

  rik_v[0] = x3 - x1;
  rik_v[1] = y3 - y1;
  rik_v[2] = z3 - z1;

  rij_rik = (   rij_v[0] * rik_v[0] 
              + rij_v[1] * rik_v[1]
              + rij_v[2] * rik_v[2] );

  rijrik = rij * rik;

  rij_p2 = pow(rij, 2.0);
  rik_p2 = pow(rik, 2.0);

  rij_p3 = pow(rij, 3.0);
  rik_p3 = pow(rik, 3.0);

  rijrik_p3 = rij * rik_p3;
  rikrij_p3 = rik * rij_p3;

  inv_rijrik = 1.0 / rijrik;

  inv_rijrik_p3 = 1.0 / rijrik_p3;
  inv_rikrij_p3 = 1.0 / rikrij_p3;

  theta_ijk = rij_rik * inv_rijrik;

  theta_const_1 = rij_rik * inv_rikrij_p3;
  theta_const_2 = rij_rik * inv_rijrik_p3;

  theta_ijk_dx = - ( rij_v[0] + rik_v[0] ) * inv_rijrik
                 + rij_v[0] * theta_const_1
                 + rik_v[0] * theta_const_2;
  theta_ijk_dy = - ( rij_v[1] + rik_v[1] ) * inv_rijrik
                 + rij_v[1] * theta_const_1
                 + rik_v[1] * theta_const_2;
  theta_ijk_dz = - ( rij_v[2] + rik_v[2] ) * inv_rijrik
                 + rij_v[2] * theta_const_1
                 + rik_v[2] * theta_const_2;

  theta_ijk_dx1 =   ( rik_v[0] * inv_rijrik ) 
                  - theta_const_1 * rij_v[0];
  theta_ijk_dy1 =   ( rik_v[1] * inv_rijrik ) 
                  - theta_const_1 * rij_v[1];
  theta_ijk_dz1 =   ( rik_v[2] * inv_rijrik ) 
                  - theta_const_1 * rij_v[2];

  theta_ijk_dx2 =   ( rij_v[0] * inv_rijrik ) 
                  - theta_const_2 * rik_v[0];
  theta_ijk_dy2 =   ( rij_v[1] * inv_rijrik ) 
                  - theta_const_2 * rik_v[1];
  theta_ijk_dz2 =   ( rij_v[2] * inv_rijrik ) 
                  - theta_const_2 * rik_v[2];


  fc1 = fc(rij);
  fc2 = fc(rik);
  
  fc_dx12 = fc_deriv(rij, x1, x2);
  fc_dy12 = fc_deriv(rij, y1, y2);
  fc_dz12 = fc_deriv(rij, z1, z2);

  fc_dx13 = fc_deriv(rik, x1, x3);
  fc_dy13 = fc_deriv(rik, y1, y3);
  fc_dz13 = fc_deriv(rik, z1, z3);

  g5_const_x =   fc_dx12 * fc2 + fc1 * fc_dx13;
  g5_const_y =   fc_dy12 * fc2 + fc1 * fc_dy13;
  g5_const_z =   fc_dz12 * fc2 + fc1 * fc_dz13;

  fc_ij_ik = fc1 * fc2;


  for (int ig5 = 0; ig5 < num_g5; ig5++) {

    exp_ij_ik = exp( -eta5[ig5] * ( rij_p2 + rik_p2 ) );
    //exp_ij_ik = exp( -eta5[ig5] * ( pow(rij, 2.0) + pow(rik, 2.0) ) );
 
    g5_const_8 = 1.0 + lambda5[ig5] * theta_ijk;
    g5_const_1 = pow( g5_const_8, zeta5[ig5] );
    g5_const_2 = pow( g5_const_8, zeta5[ig5] - 1.0 );
    g5_const_7 = g5_const_1 * exp_ij_ik;

    g5[ig5] += g5_const_7 * fc_ij_ik;

    g5_const_3 = exp_ij_ik * fc_ij_ik;
    g5_const_4 = g5_const_2 * g5_const_3;

    g5_const_x0 = g5_const_4 * theta_ijk_dx; 
    g5_const_y0 = g5_const_4 * theta_ijk_dy; 
    g5_const_z0 = g5_const_4 * theta_ijk_dz;

    g5_const_x1 = g5_const_4 * theta_ijk_dx1;
    g5_const_y1 = g5_const_4 * theta_ijk_dy1;
    g5_const_z1 = g5_const_4 * theta_ijk_dz1;

    g5_const_x2 = g5_const_4 * theta_ijk_dx2;
    g5_const_y2 = g5_const_4 * theta_ijk_dy2;
    g5_const_z2 = g5_const_4 * theta_ijk_dz2;

    g5_const_5 = zeta5[ig5] * lambda5[ig5];
    g5_const_6 = -2.0 * g5_const_1 * eta5[ig5] * g5_const_3;

//  (i,i)
    g5_deriv[i][0][ig5] +=   g5_const_5 * g5_const_x0
                           - g5_const_6 * ( rij_v[0] + rik_v[0] )
                           + g5_const_7 * g5_const_x;
    g5_deriv[i][1][ig5] +=   g5_const_5 * g5_const_y0
                           - g5_const_6 * ( rij_v[1] + rik_v[1] )
                           + g5_const_7 * g5_const_y;
    g5_deriv[i][2][ig5] +=   g5_const_5 * g5_const_z0
                           - g5_const_6 * ( rij_v[2] + rik_v[2] )
                           + g5_const_7 * g5_const_z;

//  (i,j)
    g5_deriv[j][0][ig5] +=   g5_const_5 * g5_const_x1
                           + g5_const_6 * rij_v[0] 
                           - g5_const_7 * fc_dx12 * fc2;
    g5_deriv[j][1][ig5] +=   g5_const_5 * g5_const_y1
                           + g5_const_6 * rij_v[1] 
                           - g5_const_7 * fc_dy12 * fc2;
    g5_deriv[j][2][ig5] +=   g5_const_5 * g5_const_z1
                           + g5_const_6 * rij_v[2]
                           - g5_const_7 * fc_dz12 * fc2;

//  (i,k)
    g5_deriv[k][0][ig5] +=   g5_const_5 * g5_const_x2
                           + g5_const_6 * rik_v[0]
                           - g5_const_7 * fc_dx13 * fc1;

    g5_deriv[k][1][ig5] +=   g5_const_5 * g5_const_y2
                           + g5_const_6 * rik_v[1]
                           - g5_const_7 * fc_dy13 * fc1;

    g5_deriv[k][2][ig5] +=   g5_const_5 * g5_const_z2
                           + g5_const_6 * rik_v[2]
                           - g5_const_7 * fc_dz13 * fc1;

  } // for ig5


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_e_a(double *(&g2_a_a),  int num_g2_a_a,
                          double *(&g2_a_b),  int num_g2_a_b,
                          double *(&g2_a_c),  int num_g2_a_c,
                          double *(&g5_a_aa), int num_g5_a_aa, 
                          double *(&g5_a_ab), int num_g5_a_ab, 
                          double *(&g5_a_ac), int num_g5_a_ac, 
                          double *(&g5_a_bb), int num_g5_a_bb, 
                          double *(&g5_a_bc), int num_g5_a_bc, 
                          double *(&g5_a_cc), int num_g5_a_cc, 
                          double *(&hidden_a), double &eng )
// eng: energy
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  //g5_a_aa
  for (ig5 = 0; ig5 < num_g5_a_aa; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_aa[ig5] );
    g5_a_aa[ig5] = g5_a_aa[ig5] * g5_const_3;
  } // ig5
  //g5_a_ab
  for (ig5 = 0; ig5 < num_g5_a_ab; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_ab[ig5] );
    g5_a_ab[ig5] = g5_a_ab[ig5] * g5_const_3;
  } // ig5
  //g5_a_ac
  for (ig5 = 0; ig5 < num_g5_a_ac; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_ac[ig5] );
    g5_a_ac[ig5] = g5_a_ac[ig5] * g5_const_3;
  } // ig5
  //g5_a_bb
  for (ig5 = 0; ig5 < num_g5_a_bb; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_bb[ig5] );
    g5_a_bb[ig5] = g5_a_bb[ig5] * g5_const_3;
  } // ig5
  //g5_a_bc
  for (ig5 = 0; ig5 < num_g5_a_bc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_bc[ig5] );
    g5_a_bc[ig5] = g5_a_bc[ig5] * g5_const_3;
  } // ig5
  //g5_a_cc
  for (ig5 = 0; ig5 < num_g5_a_cc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_a_cc[ig5] );
    g5_a_cc[ig5] = g5_a_cc[ig5] * g5_const_3;
  } // ig5


  /*------------------------
     preparation of input
  ------------------------*/
  // g2_a_* & g5_a_**
  input_a[0] = 1.0;
  // g2_a_a
  for (j = 0; j < num_g2_a_a; j++) {
    input_a[j+1] = ( 2.0*( g2_a_a[j] - g2_a_a_min[j] ) )/
                   ( g2_a_a_max[j] - g2_a_a_min[j] ) - 1.0;
  }
  // g2_a_b
  dummyi = num_g2_a_a;
  for (j = 0; j < num_g2_a_b; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g2_a_b[j] - g2_a_b_min[j] ) )/
                   ( g2_a_b_max[j] - g2_a_b_min[j] ) - 1.0;
  }
  // g2_a_c
  dummyi += num_g2_a_b;
  for (j = 0; j < num_g2_a_c; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g2_a_c[j] - g2_a_c_min[j] ) )/
                   ( g2_a_c_max[j] - g2_a_c_min[j] ) - 1.0;
  }
  // g5_a_aa
  dummyi += num_g2_a_c;
  for (j = 0; j < num_g5_a_aa; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_aa[j] - g5_a_aa_min[j] ) )/
                   ( g5_a_aa_max[j] - g5_a_aa_min[j] ) - 1.0;
  }
  // g5_a_ab
  dummyi += num_g5_a_aa;
  for (j = 0; j < num_g5_a_ab; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_ab[j] - g5_a_ab_min[j] ) )/
                   ( g5_a_ab_max[j] - g5_a_ab_min[j] ) - 1.0;
  }
  // g5_a_ac
  dummyi += num_g5_a_ab;
  for (j = 0; j < num_g5_a_ac; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_ac[j] - g5_a_ac_min[j] ) )/
                   ( g5_a_ac_max[j] - g5_a_ac_min[j] ) - 1.0;
  }
  // g5_a_bb
  dummyi += num_g5_a_ac;
  for (j = 0; j < num_g5_a_bb; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_bb[j] - g5_a_bb_min[j] ) )/
                   ( g5_a_bb_max[j] - g5_a_bb_min[j] ) - 1.0;
  }
  // g5_a_bc
  dummyi += num_g5_a_bb;
  for (j = 0; j < num_g5_a_bc; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_bc[j] - g5_a_bc_min[j] ) )/
                   ( g5_a_bc_max[j] - g5_a_bc_min[j] ) - 1.0;
  }
  // g5_a_cc
  dummyi += num_g5_a_bc;
  for (j = 0; j < num_g5_a_cc; j++) {
    k = dummyi + j;
    input_a[k+1] = ( 2.0*( g5_a_cc[j] - g5_a_cc_min[j] ) )/
                   ( g5_a_cc_max[j] - g5_a_cc_min[j] ) - 1.0;
  }


  // initialization
  for (n = 0; n < nhidden_a; n++) hidden_a[n] = 0.0;


  /*---------------
     energy loop
  ---------------*/
  if (nlayer_a == 4) {
    int l1,l2,h,h_pre,w,w_pre;
    h_pre = 0; w_pre = 0;
    /*------------------
       1st hidden layer
    --------------------*/
    hidden_a[0] = 1.0;
    for (l2 = 1; l2 <= layer_a[1]; l2++) {
      h = l2;
      for (l1 = 0; l1 <= layer_a[0]; l1++) {
        w = ( l2 - 1 )*( layer_a[0] + 1 ) + ( l1 );
        hidden_a[h] += weight_a[w] * input_a[l1];
      } // for l1
    } // for l2
    for (l2 = 1; l2 <= layer_a[1]; l2++) {
      h = l2;
      hidden_a[h] = tanh( hidden_a[h] );
    } // for l2
    /*-------------------------------------------
       2nd hidden layer
    ---------------------------------------------*/
    w_pre = ( layer_a[0] + 1 ) * layer_a[1];
    h_pre = layer_a[1] + 1;
    hidden_a[h_pre] = 1.0;
    for (l2 = 1; l2 <= layer_a[2]; l2++) {
      h = h_pre + l2;
      for (l1 = 0; l1 <= layer_a[1]; l1++) {
        w = w_pre + ( l2 - 1 )*( layer_a[1] + 1 ) + l1 ;
        hidden_a[h] += weight_a[w] * hidden_a[l1];
      } // for l1
    } // for l2
    for (l2 = 1; l2 <= layer_a[2]; l2++) {
      h = h_pre + l2;
      hidden_a[h] = tanh( hidden_a[h] );
    } // for l2
    /*---------------
       atomic energy
    -----------------*/
    w_pre =   ( layer_a[0] + 1 ) * layer_a[1] 
            + ( layer_a[ nlayer_a - 3 ] + 1 ) * layer_a[ nlayer_a - 2 ];
    eng = 0.0;
    for (l1 = 0; l1 <= layer_a[ nlayer_a - 2 ]; l1++) {
      w = w_pre + l1;
      h = ( layer_a[1] + 1 ) + l1;
      eng += weight_a[w] * hidden_a[h];
    } // for l1


  } // if nlayer_a == 4


  /*------------------------
     deallocate variables
  ------------------------*/ 

  //for (n = 0; n < nhidden_a; n++) hidden_a[n] = 0.0;

  for (n = 0; n < num_g2_a_a; n++) g2_a_a[n] = 0.0;
  for (n = 0; n < num_g2_a_b; n++) g2_a_b[n] = 0.0;
  for (n = 0; n < num_g2_a_c; n++) g2_a_c[n] = 0.0;
  for (n = 0; n < num_g5_a_aa; n++) g5_a_aa[n] = 0.0;
  for (n = 0; n < num_g5_a_ab; n++) g5_a_ab[n] = 0.0;
  for (n = 0; n < num_g5_a_ac; n++) g5_a_ac[n] = 0.0;
  for (n = 0; n < num_g5_a_bb; n++) g5_a_bb[n] = 0.0;
  for (n = 0; n < num_g5_a_bc; n++) g5_a_bc[n] = 0.0;
  for (n = 0; n < num_g5_a_cc; n++) g5_a_cc[n] = 0.0;


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_f_a(double ***(&g2_deriv_a_a),  int num_g2_a_a,
                          double ***(&g2_deriv_a_b),  int num_g2_a_b,
                          double ***(&g2_deriv_a_c),  int num_g2_a_c,
                          double ***(&g5_deriv_a_aa), int num_g5_a_aa, 
                          double ***(&g5_deriv_a_ab), int num_g5_a_ab, 
                          double ***(&g5_deriv_a_ac), int num_g5_a_ac, 
                          double ***(&g5_deriv_a_bb), int num_g5_a_bb, 
                          double ***(&g5_deriv_a_bc), int num_g5_a_bc, 
                          double ***(&g5_deriv_a_cc), int num_g5_a_cc, 
                          int ii, double *(&hidden_a), double &fx, double &fy, double &fz )
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  fx = 0.0; fy = 0.0; fz = 0.0;

    for (ig5 = 0; ig5 < num_g5_a_aa; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_aa[ig5] );
      g5_deriv_a_aa[ii][0][ig5] = g5_deriv_a_aa[ii][0][ig5] * g5_const_3;
      g5_deriv_a_aa[ii][1][ig5] = g5_deriv_a_aa[ii][1][ig5] * g5_const_3;
      g5_deriv_a_aa[ii][2][ig5] = g5_deriv_a_aa[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_a_ab; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_ab[ig5] );
      g5_deriv_a_ab[ii][0][ig5] = g5_deriv_a_ab[ii][0][ig5] * g5_const_3;
      g5_deriv_a_ab[ii][1][ig5] = g5_deriv_a_ab[ii][1][ig5] * g5_const_3;
      g5_deriv_a_ab[ii][2][ig5] = g5_deriv_a_ab[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_a_ac; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_ac[ig5] );
      g5_deriv_a_ac[ii][0][ig5] = g5_deriv_a_ac[ii][0][ig5] * g5_const_3;
      g5_deriv_a_ac[ii][1][ig5] = g5_deriv_a_ac[ii][1][ig5] * g5_const_3;
      g5_deriv_a_ac[ii][2][ig5] = g5_deriv_a_ac[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_a_bb; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_bb[ig5] );
      g5_deriv_a_bb[ii][0][ig5] = g5_deriv_a_bb[ii][0][ig5] * g5_const_3;
      g5_deriv_a_bb[ii][1][ig5] = g5_deriv_a_bb[ii][1][ig5] * g5_const_3;
      g5_deriv_a_bb[ii][2][ig5] = g5_deriv_a_bb[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_a_bc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_bc[ig5] );
      g5_deriv_a_bc[ii][0][ig5] = g5_deriv_a_bc[ii][0][ig5] * g5_const_3;
      g5_deriv_a_bc[ii][1][ig5] = g5_deriv_a_bc[ii][1][ig5] * g5_const_3;
      g5_deriv_a_bc[ii][2][ig5] = g5_deriv_a_bc[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_a_cc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_a_cc[ig5] );
      g5_deriv_a_cc[ii][0][ig5] = g5_deriv_a_cc[ii][0][ig5] * g5_const_3;
      g5_deriv_a_cc[ii][1][ig5] = g5_deriv_a_cc[ii][1][ig5] * g5_const_3;
      g5_deriv_a_cc[ii][2][ig5] = g5_deriv_a_cc[ii][2][ig5] * g5_const_3;
    } // ig5


  /*------------------------
     preparation of input
  ------------------------*/
  // g2_deriv_a_* & g5_deriv_a_**
    // g2_deriv_a_a
    for (k = 0; k < num_g2_a_a; k++) {
      dGdx_a[k] = 2.0*g2_deriv_a_a[ii][0][k]/( g2_a_a_max[k] - g2_a_a_min[k] );
      dGdy_a[k] = 2.0*g2_deriv_a_a[ii][1][k]/( g2_a_a_max[k] - g2_a_a_min[k] );
      dGdz_a[k] = 2.0*g2_deriv_a_a[ii][2][k]/( g2_a_a_max[k] - g2_a_a_min[k] );
    } // for k
    // g2_deriv_a_b
    dummyi = num_g2_a_a;
    for (k = 0; k < num_g2_a_b; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g2_deriv_a_b[ii][0][k]/( g2_a_b_max[k] - g2_a_b_min[k] );
      dGdy_a[m] = 2.0*g2_deriv_a_b[ii][1][k]/( g2_a_b_max[k] - g2_a_b_min[k] );
      dGdz_a[m] = 2.0*g2_deriv_a_b[ii][2][k]/( g2_a_b_max[k] - g2_a_b_min[k] );
    } // for k
    // g2_deriv_a_c
    dummyi += num_g2_a_b;
    for (k = 0; k < num_g2_a_c; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g2_deriv_a_c[ii][0][k]/( g2_a_c_max[k] - g2_a_c_min[k] );
      dGdy_a[m] = 2.0*g2_deriv_a_c[ii][1][k]/( g2_a_c_max[k] - g2_a_c_min[k] );
      dGdz_a[m] = 2.0*g2_deriv_a_c[ii][2][k]/( g2_a_c_max[k] - g2_a_c_min[k] );
    } // for k
    // g5_deriv_a_aa
    dummyi += num_g2_a_c;
    for (k = 0; k < num_g5_a_aa; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_aa[ii][0][k]/( g5_a_aa_max[k] - g5_a_aa_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_aa[ii][1][k]/( g5_a_aa_max[k] - g5_a_aa_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_aa[ii][2][k]/( g5_a_aa_max[k] - g5_a_aa_min[k] );
    } // for k
    // g5_deriv_a_ab
    dummyi += num_g5_a_aa;
    for (k = 0; k < num_g5_a_ab; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_ab[ii][0][k]/( g5_a_ab_max[k] - g5_a_ab_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_ab[ii][1][k]/( g5_a_ab_max[k] - g5_a_ab_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_ab[ii][2][k]/( g5_a_ab_max[k] - g5_a_ab_min[k] );
    } // for k
    // g5_deriv_a_ac
    dummyi += num_g5_a_ab;
    for (k = 0; k < num_g5_a_ac; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_ac[ii][0][k]/( g5_a_ac_max[k] - g5_a_ac_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_ac[ii][1][k]/( g5_a_ac_max[k] - g5_a_ac_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_ac[ii][2][k]/( g5_a_ac_max[k] - g5_a_ac_min[k] );
    } // for k
    // g5_deriv_a_bb
    dummyi += num_g5_a_ac;
    for (k = 0; k < num_g5_a_bb; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_bb[ii][0][k]/( g5_a_bb_max[k] - g5_a_bb_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_bb[ii][1][k]/( g5_a_bb_max[k] - g5_a_bb_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_bb[ii][2][k]/( g5_a_bb_max[k] - g5_a_bb_min[k] );
    } // for k
    // g5_deriv_a_bc
    dummyi += num_g5_a_bb;
    for (k = 0; k < num_g5_a_bc; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_bc[ii][0][k]/( g5_a_bc_max[k] - g5_a_bc_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_bc[ii][1][k]/( g5_a_bc_max[k] - g5_a_bc_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_bc[ii][2][k]/( g5_a_bc_max[k] - g5_a_bc_min[k] );
    } // for k
    // g5_deriv_a_cc
    dummyi += num_g5_a_bc;
    for (k = 0; k < num_g5_a_cc; k++) {
      m = dummyi + k;
      dGdx_a[m] = 2.0*g5_deriv_a_cc[ii][0][k]/( g5_a_cc_max[k] - g5_a_cc_min[k] );
      dGdy_a[m] = 2.0*g5_deriv_a_cc[ii][1][k]/( g5_a_cc_max[k] - g5_a_cc_min[k] );
      dGdz_a[m] = 2.0*g5_deriv_a_cc[ii][2][k]/( g5_a_cc_max[k] - g5_a_cc_min[k] );
    } // for k


  /*--------------
     force loop
  --------------*/
  if (nlayer_a == 4) {
    int w1,w2,w2_pre,w3,w3_pre,h1,h2,h2_pre;
    double dummy,dummy1,dummy2;
    w2_pre = ( 1 + layer_a[0] ) * layer_a[1];
    w3_pre = w2_pre + ( 1 + layer_a[1] ) * layer_a[2];
    h2_pre = 1 + layer_a[1];
    for (int l0 = 1; l0 <= layer_a[0]; l0++){
    for (int l1 = 1; l1 <= layer_a[1]; l1++){
      w1 = ( 1 + layer_a[0] )*( l1 - 1 ) + l0;
      h1 = l1;
      dummy1 = weight_a[w1] * ( 1.0 - pow(hidden_a[h1], 2.0) );
    for (int l2 = 1; l2 <= layer_a[2]; l2++){
      w2 = w2_pre + ( 1 + layer_a[1] )*( l2 - 1 ) + l1;
      h2 = h2_pre + l2;
      dummy2 = weight_a[w2] * ( 1.0 - pow(hidden_a[h2], 2.0) );
    for (int l3 = 1; l3 <= layer_a[3]; l3++){
      w3 = w3_pre + ( 1 + layer_a[2] )*( l3 - 1 ) + l2;
      dummy = weight_a[w3] * dummy1 * dummy2;
      //dummy =   weight_a[w1] * ( 1.0 - pow(hidden_a[h1], 2.0) ) 
      //        * weight_a[w2] * ( 1.0 - pow(hidden_a[h2], 2.0) )
      //        * weight_a[w3];
      fx += -dGdx_a[l0-1] * dummy;
      fy += -dGdy_a[l0-1] * dummy;
      fz += -dGdz_a[l0-1] * dummy;
    } // for l4
    } // for l3
    } // for l2
    } // for l1

  } // if nlayer_a == 4


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_e_b(double *(&g2_b_a),  int num_g2_b_a,
                          double *(&g2_b_b),  int num_g2_b_b,
                          double *(&g2_b_c),  int num_g2_b_c,
                          double *(&g5_b_aa), int num_g5_b_aa, 
                          double *(&g5_b_ab), int num_g5_b_ab, 
                          double *(&g5_b_ac), int num_g5_b_ac, 
                          double *(&g5_b_bb), int num_g5_b_bb, 
                          double *(&g5_b_bc), int num_g5_b_bc, 
                          double *(&g5_b_cc), int num_g5_b_cc, 
                          double *(&hidden_b), double &eng )
// eng: energy
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  //g5_b_aa
  for (ig5 = 0; ig5 < num_g5_b_aa; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_aa[ig5] );
    g5_b_aa[ig5] = g5_b_aa[ig5] * g5_const_3;
  } // ig5
  //g5_b_ab
  for (ig5 = 0; ig5 < num_g5_b_ab; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_ab[ig5] );
    g5_b_ab[ig5] = g5_b_ab[ig5] * g5_const_3;
  } // ig5
  //g5_b_ac
  for (ig5 = 0; ig5 < num_g5_b_ac; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_ac[ig5] );
    g5_b_ac[ig5] = g5_b_ac[ig5] * g5_const_3;
  } // ig5
  //g5_b_bb
  for (ig5 = 0; ig5 < num_g5_b_bb; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_bb[ig5] );
    g5_b_bb[ig5] = g5_b_bb[ig5] * g5_const_3;
  } // ig5
  //g5_b_bc
  for (ig5 = 0; ig5 < num_g5_b_bc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_bc[ig5] );
    g5_b_bc[ig5] = g5_b_bc[ig5] * g5_const_3;
  } // ig5
  //g5_b_cc
  for (ig5 = 0; ig5 < num_g5_b_cc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_b_cc[ig5] );
    g5_b_cc[ig5] = g5_b_cc[ig5] * g5_const_3;
  } // ig5


  /*------------------------
     preparation of input
  ------------------------*/
  // g2_b_* & g5_b_**
  input_b[0] = 1.0;
  // g2_b_a
  for (j = 0; j < num_g2_b_a; j++) {
    input_b[j+1] = ( 2.0*( g2_b_a[j] - g2_b_a_min[j] ) )/
                   ( g2_b_a_max[j] - g2_b_a_min[j] ) - 1.0;
  }
  // g2_b_b
  dummyi = num_g2_b_a;
  for (j = 0; j < num_g2_b_b; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g2_b_b[j] - g2_b_b_min[j] ) )/
                   ( g2_b_b_max[j] - g2_b_b_min[j] ) - 1.0;
  }
  // g2_b_c
  dummyi += num_g2_b_b;
  for (j = 0; j < num_g2_b_c; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g2_b_c[j] - g2_b_c_min[j] ) )/
                   ( g2_b_c_max[j] - g2_b_c_min[j] ) - 1.0;
  }
  // g5_b_aa
  dummyi += num_g2_b_c;
  for (j = 0; j < num_g5_b_aa; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_aa[j] - g5_b_aa_min[j] ) )/
                   ( g5_b_aa_max[j] - g5_b_aa_min[j] ) - 1.0;
  }
  // g5_b_ab
  dummyi += num_g5_b_aa;
  for (j = 0; j < num_g5_b_ab; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_ab[j] - g5_b_ab_min[j] ) )/
                   ( g5_b_ab_max[j] - g5_b_ab_min[j] ) - 1.0;
  }
  // g5_b_ac
  dummyi += num_g5_b_ab;
  for (j = 0; j < num_g5_b_ac; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_ac[j] - g5_b_ac_min[j] ) )/
                   ( g5_b_ac_max[j] - g5_b_ac_min[j] ) - 1.0;
  }
  // g5_b_bb
  dummyi += num_g5_b_ac;
  for (j = 0; j < num_g5_b_bb; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_bb[j] - g5_b_bb_min[j] ) )/
                   ( g5_b_bb_max[j] - g5_b_bb_min[j] ) - 1.0;
  }
  // g5_b_bc
  dummyi += num_g5_b_bb;
  for (j = 0; j < num_g5_b_bc; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_bc[j] - g5_b_bc_min[j] ) )/
                   ( g5_b_bc_max[j] - g5_b_bc_min[j] ) - 1.0;
  }
  // g5_b_cc
  dummyi += num_g5_b_bc;
  for (j = 0; j < num_g5_b_cc; j++) {
    k = dummyi + j;
    input_b[k+1] = ( 2.0*( g5_b_cc[j] - g5_b_cc_min[j] ) )/
                   ( g5_b_cc_max[j] - g5_b_cc_min[j] ) - 1.0;
  }


  // initialization
  for (n = 0; n < nhidden_b; n++) hidden_b[n] = 0.0;


  /*---------------
     energy loop
  ---------------*/
  if (nlayer_b == 4) {
    int l1,l2,h,h_pre,w,w_pre;
    h_pre = 0; w_pre = 0;
    /*------------------
       1st hidden layer
    --------------------*/
    hidden_b[0] = 1.0;
    for (l2 = 1; l2 <= layer_b[1]; l2++) {
      h = l2;
      for (l1 = 0; l1 <= layer_b[0]; l1++) {
        w = ( l2 - 1 )*( layer_b[0] + 1 ) + ( l1 );
        hidden_b[h] += weight_b[w] * input_b[l1];
      } // for l1
    } // for l2
    for (l2 = 1; l2 <= layer_b[1]; l2++) {
      h = l2;
      hidden_b[h] = tanh( hidden_b[h] );
    } // for l2
    /*-------------------------------------------
       2nd hidden layer
    ---------------------------------------------*/
    w_pre = ( layer_b[0] + 1 ) * layer_b[1];
    h_pre = layer_b[1] + 1;
    hidden_b[h_pre] = 1.0;
    for (l2 = 1; l2 <= layer_b[2]; l2++) {
      h = h_pre + l2;
      for (l1 = 0; l1 <= layer_b[1]; l1++) {
        w = w_pre + ( l2 - 1 )*( layer_b[1] + 1 ) + l1 ;
        hidden_b[h] += weight_b[w] * hidden_b[l1];
      } // for l1
    } // for l2
    for (l2 = 1; l2 <= layer_b[2]; l2++) {
      h = h_pre + l2;
      hidden_b[h] = tanh( hidden_b[h] );
    } // for l2
    /*---------------
       atomic energy
    -----------------*/
    w_pre =   ( layer_b[0] + 1 ) * layer_b[1] 
            + ( layer_b[ nlayer_b - 3 ] + 1 ) * layer_b[ nlayer_b - 2 ];
    eng = 0.0;
    for (l1 = 0; l1 <= layer_b[ nlayer_b - 2 ]; l1++) {
      w = w_pre + l1;
      h = layer_b[1] + 1 + l1;
      eng += weight_b[w] * hidden_b[h];
    } // for l1


  } // if nlayer_b == 4

  /*------------------------
     deallocate variables
  ------------------------*/ 
  //for (n = 0; n < nhidden_b; n++) hidden_b[n] = 0.0;

  for (n = 0; n < num_g2_b_a; n++)  g2_b_a[n] = 0.0;
  for (n = 0; n < num_g2_b_b; n++)  g2_b_b[n] = 0.0;
  for (n = 0; n < num_g2_b_c; n++)  g2_b_c[n] = 0.0;
  for (n = 0; n < num_g5_b_aa; n++) g5_b_aa[n] = 0.0;
  for (n = 0; n < num_g5_b_ab; n++) g5_b_ab[n] = 0.0;
  for (n = 0; n < num_g5_b_ac; n++) g5_b_ac[n] = 0.0;
  for (n = 0; n < num_g5_b_bb; n++) g5_b_bb[n] = 0.0;
  for (n = 0; n < num_g5_b_bc; n++) g5_b_bc[n] = 0.0;
  for (n = 0; n < num_g5_b_cc; n++) g5_b_cc[n] = 0.0;


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_f_b(double ***(&g2_deriv_b_a),  int num_g2_b_a,
                          double ***(&g2_deriv_b_b),  int num_g2_b_b,
                          double ***(&g2_deriv_b_c),  int num_g2_b_c,
                          double ***(&g5_deriv_b_aa), int num_g5_b_aa, 
                          double ***(&g5_deriv_b_ab), int num_g5_b_ab, 
                          double ***(&g5_deriv_b_ac), int num_g5_b_ac, 
                          double ***(&g5_deriv_b_bb), int num_g5_b_bb, 
                          double ***(&g5_deriv_b_bc), int num_g5_b_bc, 
                          double ***(&g5_deriv_b_cc), int num_g5_b_cc, 
                          int ii, double *(&hidden_b), double &fx, double &fy, double &fz )
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  fx = 0.0; fy = 0.0; fz = 0.0;

    for (ig5 = 0; ig5 < num_g5_b_aa; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_aa[ig5] );
      g5_deriv_b_aa[ii][0][ig5] = g5_deriv_b_aa[ii][0][ig5] * g5_const_3;
      g5_deriv_b_aa[ii][1][ig5] = g5_deriv_b_aa[ii][1][ig5] * g5_const_3;
      g5_deriv_b_aa[ii][2][ig5] = g5_deriv_b_aa[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_b_ab; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_ab[ig5] );
      g5_deriv_b_ab[ii][0][ig5] = g5_deriv_b_ab[ii][0][ig5] * g5_const_3;
      g5_deriv_b_ab[ii][1][ig5] = g5_deriv_b_ab[ii][1][ig5] * g5_const_3;
      g5_deriv_b_ab[ii][2][ig5] = g5_deriv_b_ab[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_b_ac; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_ac[ig5] );
      g5_deriv_b_ac[ii][0][ig5] = g5_deriv_b_ac[ii][0][ig5] * g5_const_3;
      g5_deriv_b_ac[ii][1][ig5] = g5_deriv_b_ac[ii][1][ig5] * g5_const_3;
      g5_deriv_b_ac[ii][2][ig5] = g5_deriv_b_ac[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_b_bb; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_bb[ig5] );
      g5_deriv_b_bb[ii][0][ig5] = g5_deriv_b_bb[ii][0][ig5] * g5_const_3;
      g5_deriv_b_bb[ii][1][ig5] = g5_deriv_b_bb[ii][1][ig5] * g5_const_3;
      g5_deriv_b_bb[ii][2][ig5] = g5_deriv_b_bb[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_b_bc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_bc[ig5] );
      g5_deriv_b_bc[ii][0][ig5] = g5_deriv_b_bc[ii][0][ig5] * g5_const_3;
      g5_deriv_b_bc[ii][1][ig5] = g5_deriv_b_bc[ii][1][ig5] * g5_const_3;
      g5_deriv_b_bc[ii][2][ig5] = g5_deriv_b_bc[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_b_cc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_b_cc[ig5] );
      g5_deriv_b_cc[ii][0][ig5] = g5_deriv_b_cc[ii][0][ig5] * g5_const_3;
      g5_deriv_b_cc[ii][1][ig5] = g5_deriv_b_cc[ii][1][ig5] * g5_const_3;
      g5_deriv_b_cc[ii][2][ig5] = g5_deriv_b_cc[ii][2][ig5] * g5_const_3;
    } // ig5


  /*------------------------
     preparation of input
  ------------------------*/
  // g2_deriv_b_* & g5_deriv_b_**
    // g2_deriv_b_a
    for (k = 0; k < num_g2_b_a; k++) {
      dGdx_b[k] = 2.0*g2_deriv_b_a[ii][0][k]/( g2_b_a_max[k] - g2_b_a_min[k] );
      dGdy_b[k] = 2.0*g2_deriv_b_a[ii][1][k]/( g2_b_a_max[k] - g2_b_a_min[k] );
      dGdz_b[k] = 2.0*g2_deriv_b_a[ii][2][k]/( g2_b_a_max[k] - g2_b_a_min[k] );
    } // for k
    // g2_deriv_b_b
    dummyi = num_g2_b_a;
    for (k = 0; k < num_g2_b_b; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g2_deriv_b_b[ii][0][k]/( g2_b_b_max[k] - g2_b_b_min[k] );
      dGdy_b[m] = 2.0*g2_deriv_b_b[ii][1][k]/( g2_b_b_max[k] - g2_b_b_min[k] );
      dGdz_b[m] = 2.0*g2_deriv_b_b[ii][2][k]/( g2_b_b_max[k] - g2_b_b_min[k] );
    } // for k
    // g2_deriv_b_c
    dummyi += num_g2_b_b;
    for (int k = 0; k < num_g2_b_c; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g2_deriv_b_c[ii][0][k]/( g2_b_c_max[k] - g2_b_c_min[k] );
      dGdy_b[m] = 2.0*g2_deriv_b_c[ii][1][k]/( g2_b_c_max[k] - g2_b_c_min[k] );
      dGdz_b[m] = 2.0*g2_deriv_b_c[ii][2][k]/( g2_b_c_max[k] - g2_b_c_min[k] );
    } // for k
    // g5_deriv_b_aa
    dummyi += num_g2_b_c;
    for (int k = 0; k < num_g5_b_aa; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_aa[ii][0][k]/( g5_b_aa_max[k] - g5_b_aa_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_aa[ii][1][k]/( g5_b_aa_max[k] - g5_b_aa_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_aa[ii][2][k]/( g5_b_aa_max[k] - g5_b_aa_min[k] );
    } // for k
    // g5_deriv_b_ab
    dummyi += num_g5_b_aa;
    for (int k = 0; k < num_g5_b_ab; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_ab[ii][0][k]/( g5_b_ab_max[k] - g5_b_ab_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_ab[ii][1][k]/( g5_b_ab_max[k] - g5_b_ab_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_ab[ii][2][k]/( g5_b_ab_max[k] - g5_b_ab_min[k] );
    } // for k
    // g5_deriv_b_ac
    dummyi += num_g5_b_ab;
    for (int k = 0; k < num_g5_b_ac; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_ac[ii][0][k]/( g5_b_ac_max[k] - g5_b_ac_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_ac[ii][1][k]/( g5_b_ac_max[k] - g5_b_ac_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_ac[ii][2][k]/( g5_b_ac_max[k] - g5_b_ac_min[k] );
    } // for k
    // g5_deriv_b_bb
    dummyi += num_g5_b_ac;
    for (int k = 0; k < num_g5_b_bb; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_bb[ii][0][k]/( g5_b_bb_max[k] - g5_b_bb_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_bb[ii][1][k]/( g5_b_bb_max[k] - g5_b_bb_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_bb[ii][2][k]/( g5_b_bb_max[k] - g5_b_bb_min[k] );
    } // for k
    // g5_deriv_b_bc
    dummyi += num_g5_b_bb;
    for (int k = 0; k < num_g5_b_bc; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_bc[ii][0][k]/( g5_b_bc_max[k] - g5_b_bc_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_bc[ii][1][k]/( g5_b_bc_max[k] - g5_b_bc_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_bc[ii][2][k]/( g5_b_bc_max[k] - g5_b_bc_min[k] );
    } // for k
    // g5_deriv_b_cc
    dummyi += num_g5_b_bc;
    for (int k = 0; k < num_g5_b_cc; k++) {
      m = dummyi + k;
      dGdx_b[m] = 2.0*g5_deriv_b_cc[ii][0][k]/( g5_b_cc_max[k] - g5_b_cc_min[k] );
      dGdy_b[m] = 2.0*g5_deriv_b_cc[ii][1][k]/( g5_b_cc_max[k] - g5_b_cc_min[k] );
      dGdz_b[m] = 2.0*g5_deriv_b_cc[ii][2][k]/( g5_b_cc_max[k] - g5_b_cc_min[k] );
    } // for k


  /*--------------
     force loop
  --------------*/
  if (nlayer_b == 4) {
    int w1,w2,w2_pre,w3,w3_pre;
    int h1,h2,h2_pre;
    double dummy,dummy1,dummy2;
    w2_pre = ( 1 + layer_b[0] ) * layer_b[1];
    w3_pre = w2_pre + ( 1 + layer_b[1] ) * layer_b[2];
    h2_pre = 1 + layer_b[1];
    for (int l0 = 1; l0 <= layer_b[0]; l0++){
    for (int l1 = 1; l1 <= layer_b[1]; l1++){
      w1 = ( 1 + layer_b[0] )*( l1 - 1 ) + l0;
      h1 = l1;
      dummy1 = weight_b[w1] * ( 1.0 - pow(hidden_b[h1], 2.0) );
    for (int l2 = 1; l2 <= layer_b[2]; l2++){
      w2 = w2_pre + ( 1 + layer_b[1] )*( l2 - 1 ) + l1;
      h2 = h2_pre + l2;
      dummy2 = weight_b[w2] * ( 1.0 - pow(hidden_b[h2], 2.0) );
    for (int l3 = 1; l3 <= layer_b[3]; l3++){
      w3 = w3_pre + ( 1 + layer_b[2] )*( l3 - 1 ) + l2;
      dummy = weight_b[w3] * dummy1 * dummy2;
      //dummy =   weight_b[w1] * ( 1.0 - pow(hidden_b[h1], 2.0) ) 
      //        * weight_b[w2] * ( 1.0 - pow(hidden_b[h2], 2.0) )
      //        * weight_b[w3];
      fx += -dGdx_b[l0-1] * dummy;
      fy += -dGdy_b[l0-1] * dummy;
      fz += -dGdz_b[l0-1] * dummy;
    } // for l4
    } // for l3
    } // for l2
    } // for l1

  } // if nlayer_b == 4


}

/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_e_c(double *(&g2_c_a),  int num_g2_c_a,
                          double *(&g2_c_b),  int num_g2_c_b,
                          double *(&g2_c_c),  int num_g2_c_c,
                          double *(&g5_c_aa), int num_g5_c_aa, 
                          double *(&g5_c_ab), int num_g5_c_ab, 
                          double *(&g5_c_ac), int num_g5_c_ac, 
                          double *(&g5_c_bb), int num_g5_c_bb, 
                          double *(&g5_c_bc), int num_g5_c_bc, 
                          double *(&g5_c_cc), int num_g5_c_cc, 
                          double *(&hidden_c), double &eng )
// eng: energy
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  //g5_c_aa
  for (ig5 = 0; ig5 < num_g5_c_aa; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_aa[ig5] );
    g5_c_aa[ig5] = g5_c_aa[ig5] * g5_const_3;
  } // ig5
  //g5_c_ab
  for (ig5 = 0; ig5 < num_g5_c_ab; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_ab[ig5] );
    g5_c_ab[ig5] = g5_c_ab[ig5] * g5_const_3;
  } // ig5
  //g5_c_ac
  for (ig5 = 0; ig5 < num_g5_c_ac; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_ac[ig5] );
    g5_c_ac[ig5] = g5_c_ac[ig5] * g5_const_3;
  } // ig5
  //g5_c_bb
  for (ig5 = 0; ig5 < num_g5_c_bb; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_bb[ig5] );
    g5_c_bb[ig5] = g5_c_bb[ig5] * g5_const_3;
  } // ig5
  //g5_c_bc
  for (ig5 = 0; ig5 < num_g5_c_bc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_bc[ig5] );
    g5_c_bc[ig5] = g5_c_bc[ig5] * g5_const_3;
  } // ig5
  //g5_c_cc
  for (ig5 = 0; ig5 < num_g5_c_cc; ig5++) {
    g5_const_3 = pow( 2.0, 1.0 - zeta5_c_cc[ig5] );
    g5_c_cc[ig5] = g5_c_cc[ig5] * g5_const_3;
  } // ig5

 
  /*------------------------
     preparation of input
  ------------------------*/
  // g2_c_* & g5_c_**
  input_c[0] = 1.0;
  // g2_c_a
  for (j = 0; j < num_g2_c_a; j++) {
    input_c[j+1] = ( 2.0*( g2_c_a[j] - g2_c_a_min[j] ) )/
                   ( g2_c_a_max[j] - g2_c_a_min[j] ) - 1.0;
  }
  // g2_c_b
  dummyi = num_g2_c_a;
  for (j = 0; j < num_g2_c_b; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g2_c_b[j] - g2_c_b_min[j] ) )/
                   ( g2_c_b_max[j] - g2_c_b_min[j] ) - 1.0;
  }
  // g2_c_c
  dummyi += num_g2_c_b;
  for (j = 0; j < num_g2_c_c; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g2_c_c[j] - g2_c_c_min[j] ) )/
                   ( g2_c_c_max[j] - g2_c_c_min[j] ) - 1.0;
  }
  // g5_c_aa
  dummyi += num_g2_c_c;
  for (j = 0; j < num_g5_c_aa; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_aa[j] - g5_c_aa_min[j] ) )/
                   ( g5_c_aa_max[j] - g5_c_aa_min[j] ) - 1.0;
  }
  // g5_c_ab
  dummyi += num_g5_c_aa;
  for (j = 0; j < num_g5_c_ab; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_ab[j] - g5_c_ab_min[j] ) )/
                   ( g5_c_ab_max[j] - g5_c_ab_min[j] ) - 1.0;
  }
  // g5_c_ac
  dummyi += num_g5_c_ab;
  for (j = 0; j < num_g5_c_ac; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_ac[j] - g5_c_ac_min[j] ) )/
                   ( g5_c_ac_max[j] - g5_c_ac_min[j] ) - 1.0;
  }
  // g5_c_bb
  dummyi += num_g5_c_ac;
  for (j = 0; j < num_g5_c_bb; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_bb[j] - g5_c_bb_min[j] ) )/
                   ( g5_c_bb_max[j] - g5_c_bb_min[j] ) - 1.0;
  }
  // g5_c_bc
  dummyi += num_g5_c_bb;
  for (j = 0; j < num_g5_c_bc; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_bc[j] - g5_c_bc_min[j] ) )/
                   ( g5_c_bc_max[j] - g5_c_bc_min[j] ) - 1.0;
  }
  // g5_c_cc
  dummyi += num_g5_c_bc;
  for (j = 0; j < num_g5_c_cc; j++) {
    k = dummyi + j;
    input_c[k+1] = ( 2.0*( g5_c_cc[j] - g5_c_cc_min[j] ) )/
                   ( g5_c_cc_max[j] - g5_c_cc_min[j] ) - 1.0;
  }


  // initialization
  for (n = 0; n < nhidden_c; n++) hidden_c[n] = 0.0;
  

  /*---------------
     energy loop
  ---------------*/
  if (nlayer_c == 4) {
    int h,h_pre,w,w_pre;
    h_pre = 0; w_pre = 0;
    /*------------------
       1st hidden layer
    --------------------*/
    hidden_c[0] = 1.0;
    for (int l2 = 1; l2 <= layer_c[1]; l2++) {
      h = l2;
      for (int l1 = 0; l1 <= layer_c[0]; l1++) {
        w = ( l2 - 1 )*( layer_c[0] + 1 ) + ( l1 );
        hidden_c[h] += weight_c[w] * input_c[l1];
      } // for l1
    } // for l2
    for (int l2 = 1; l2 <= layer_c[1]; l2++) {
      h = l2;
      hidden_c[h] = tanh( hidden_c[h] );
    } // for l2
    /*-------------------------------------------
       2nd hidden layer
    ---------------------------------------------*/
    w_pre = ( layer_c[0] + 1 ) * layer_c[1];
    h_pre = layer_c[1] + 1;
    hidden_c[h_pre] = 1.0;
    for (int l2 = 1; l2 <= layer_c[2]; l2++) {
      h = h_pre + l2;
      for (int l1 = 0; l1 <= layer_c[1]; l1++) {
        w = w_pre + ( l2 - 1 )*( layer_c[1] + 1 ) + l1 ;
        hidden_c[h] += weight_c[w] * hidden_c[l1];
      } // for l1
    } // for l2
    for (int l2 = 1; l2 <= layer_c[2]; l2++) {
      h = h_pre + l2;
      hidden_c[h] = tanh( hidden_c[h] );
    } // for l2
    /*---------------
       atomic energy
    -----------------*/
    w_pre =   ( layer_c[0] + 1 ) * layer_c[1] 
            + ( layer_c[ nlayer_c - 3 ] + 1 ) * layer_c[ nlayer_c - 2 ];
    eng = 0.0;
    for (int l1 = 0; l1 <= layer_c[ nlayer_c - 2 ]; l1++) {
      w = w_pre + l1;
      h = layer_c[1] + 1 + l1;
      eng += weight_c[w] * hidden_c[h];
    } // for l1


  } // if nlayer_c == 4

  /*------------------------
     deallocate variables
  ------------------------*/ 
  //for (int n = 0; n < nhidden_c; n++) hidden_c[n] = 0.0;

  for (int n = 0; n < num_g2_c_a; n++) g2_c_a[n] = 0.0;
  for (int n = 0; n < num_g2_c_b; n++) g2_c_b[n] = 0.0;
  for (int n = 0; n < num_g2_c_c; n++) g2_c_c[n] = 0.0;
  for (int n = 0; n < num_g5_c_aa; n++) g5_c_aa[n] = 0.0;
  for (int n = 0; n < num_g5_c_ab; n++) g5_c_ab[n] = 0.0;
  for (int n = 0; n < num_g5_c_ac; n++) g5_c_ac[n] = 0.0;
  for (int n = 0; n < num_g5_c_bb; n++) g5_c_bb[n] = 0.0;
  for (int n = 0; n < num_g5_c_bc; n++) g5_c_bc[n] = 0.0;
  for (int n = 0; n < num_g5_c_cc; n++) g5_c_cc[n] = 0.0;


}
/* ---------------------------------------------------------------------- */
void PairNNP3v::hdnnp_f_c(double ***(&g2_deriv_c_a),  int num_g2_c_a,
                          double ***(&g2_deriv_c_b),  int num_g2_c_b,
                          double ***(&g2_deriv_c_c),  int num_g2_c_c,
                          double ***(&g5_deriv_c_aa), int num_g5_c_aa, 
                          double ***(&g5_deriv_c_ab), int num_g5_c_ab, 
                          double ***(&g5_deriv_c_ac), int num_g5_c_ac, 
                          double ***(&g5_deriv_c_bb), int num_g5_c_bb, 
                          double ***(&g5_deriv_c_bc), int num_g5_c_bc, 
                          double ***(&g5_deriv_c_cc), int num_g5_c_cc, 
                          int ii, double *(&hidden_c), double &fx, double &fy, double &fz )
{
  int i,j,k,m,n,ig5,dummyi;
  double g5_const_3;

  fx = 0.0; fy = 0.0; fz =0.0;

    for (ig5 = 0; ig5 < num_g5_c_aa; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_aa[ig5] );
      g5_deriv_c_aa[ii][0][ig5] = g5_deriv_c_aa[ii][0][ig5] * g5_const_3;
      g5_deriv_c_aa[ii][1][ig5] = g5_deriv_c_aa[ii][1][ig5] * g5_const_3;
      g5_deriv_c_aa[ii][2][ig5] = g5_deriv_c_aa[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_c_ab; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_ab[ig5] );
      g5_deriv_c_ab[ii][0][ig5] = g5_deriv_c_ab[ii][0][ig5] * g5_const_3;
      g5_deriv_c_ab[ii][1][ig5] = g5_deriv_c_ab[ii][1][ig5] * g5_const_3;
      g5_deriv_c_ab[ii][2][ig5] = g5_deriv_c_ab[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_c_ac; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_ac[ig5] );
      g5_deriv_c_ac[ii][0][ig5] = g5_deriv_c_ac[ii][0][ig5] * g5_const_3;
      g5_deriv_c_ac[ii][1][ig5] = g5_deriv_c_ac[ii][1][ig5] * g5_const_3;
      g5_deriv_c_ac[ii][2][ig5] = g5_deriv_c_ac[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_c_bb; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_bb[ig5] );
      g5_deriv_c_bb[ii][0][ig5] = g5_deriv_c_bb[ii][0][ig5] * g5_const_3;
      g5_deriv_c_bb[ii][1][ig5] = g5_deriv_c_bb[ii][1][ig5] * g5_const_3;
      g5_deriv_c_bb[ii][2][ig5] = g5_deriv_c_bb[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_c_bc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_bc[ig5] );
      g5_deriv_c_bc[ii][0][ig5] = g5_deriv_c_bc[ii][0][ig5] * g5_const_3;
      g5_deriv_c_bc[ii][1][ig5] = g5_deriv_c_bc[ii][1][ig5] * g5_const_3;
      g5_deriv_c_bc[ii][2][ig5] = g5_deriv_c_bc[ii][2][ig5] * g5_const_3;
    } // ig5
    for (ig5 = 0; ig5 < num_g5_c_cc; ig5++) {
      g5_const_3 = pow( 2.0, 1.0 - zeta5_c_cc[ig5] );
      g5_deriv_c_cc[ii][0][ig5] = g5_deriv_c_cc[ii][0][ig5] * g5_const_3;
      g5_deriv_c_cc[ii][1][ig5] = g5_deriv_c_cc[ii][1][ig5] * g5_const_3;
      g5_deriv_c_cc[ii][2][ig5] = g5_deriv_c_cc[ii][2][ig5] * g5_const_3;
    } // ig5

 
  /*------------------------
     preparation of input
  ------------------------*/
  // g2_deriv_c_* & g5_deriv_c_**
    // g2_deriv_c_a
    for (k = 0; k < num_g2_c_a; k++) {
      dGdx_c[k] = 2.0*g2_deriv_c_a[ii][0][k]/( g2_c_a_max[k] - g2_c_a_min[k] );
      dGdy_c[k] = 2.0*g2_deriv_c_a[ii][1][k]/( g2_c_a_max[k] - g2_c_a_min[k] );
      dGdz_c[k] = 2.0*g2_deriv_c_a[ii][2][k]/( g2_c_a_max[k] - g2_c_a_min[k] );
    } // for k
    // g2_deriv_c_b
    dummyi = num_g2_c_a;
    for (int k = 0; k < num_g2_c_b; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g2_deriv_c_b[ii][0][k]/( g2_c_b_max[k] - g2_c_b_min[k] );
      dGdy_c[m] = 2.0*g2_deriv_c_b[ii][1][k]/( g2_c_b_max[k] - g2_c_b_min[k] );
      dGdz_c[m] = 2.0*g2_deriv_c_b[ii][2][k]/( g2_c_b_max[k] - g2_c_b_min[k] );
    } // for k
    // g2_deriv_c_c
    dummyi += num_g2_c_b;
    for (int k = 0; k < num_g2_c_c; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g2_deriv_c_c[ii][0][k]/( g2_c_c_max[k] - g2_c_c_min[k] );
      dGdy_c[m] = 2.0*g2_deriv_c_c[ii][1][k]/( g2_c_c_max[k] - g2_c_c_min[k] );
      dGdz_c[m] = 2.0*g2_deriv_c_c[ii][2][k]/( g2_c_c_max[k] - g2_c_c_min[k] );
    } // for k
    // g5_deriv_c_aa
    dummyi += num_g2_c_c;
    for (int k = 0; k < num_g5_c_aa; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_aa[ii][0][k]/( g5_c_aa_max[k] - g5_c_aa_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_aa[ii][1][k]/( g5_c_aa_max[k] - g5_c_aa_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_aa[ii][2][k]/( g5_c_aa_max[k] - g5_c_aa_min[k] );
    } // for k
    // g5_deriv_c_ab
    dummyi += num_g5_c_aa;
    for (int k = 0; k < num_g5_c_ab; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_ab[ii][0][k]/( g5_c_ab_max[k] - g5_c_ab_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_ab[ii][1][k]/( g5_c_ab_max[k] - g5_c_ab_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_ab[ii][2][k]/( g5_c_ab_max[k] - g5_c_ab_min[k] );
    } // for k
    // g5_deriv_c_ac
    dummyi += num_g5_c_ab;
    for (int k = 0; k < num_g5_c_ac; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_ac[ii][0][k]/( g5_c_ac_max[k] - g5_c_ac_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_ac[ii][1][k]/( g5_c_ac_max[k] - g5_c_ac_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_ac[ii][2][k]/( g5_c_ac_max[k] - g5_c_ac_min[k] );
    } // for k
    // g5_deriv_c_bb
    dummyi += num_g5_c_ac;
    for (int k = 0; k < num_g5_c_bb; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_bb[ii][0][k]/( g5_c_bb_max[k] - g5_c_bb_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_bb[ii][1][k]/( g5_c_bb_max[k] - g5_c_bb_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_bb[ii][2][k]/( g5_c_bb_max[k] - g5_c_bb_min[k] );
    } // for k
    // g5_deriv_c_bc
    dummyi += num_g5_c_bb;
    for (int k = 0; k < num_g5_c_bc; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_bc[ii][0][k]/( g5_c_bc_max[k] - g5_c_bc_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_bc[ii][1][k]/( g5_c_bc_max[k] - g5_c_bc_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_bc[ii][2][k]/( g5_c_bc_max[k] - g5_c_bc_min[k] );
    } // for k
    // g5_deriv_c_cc
    dummyi += num_g5_c_bc;
    for (int k = 0; k < num_g5_c_cc; k++) {
      m = dummyi + k;
      dGdx_c[m] = 2.0*g5_deriv_c_cc[ii][0][k]/( g5_c_cc_max[k] - g5_c_cc_min[k] );
      dGdy_c[m] = 2.0*g5_deriv_c_cc[ii][1][k]/( g5_c_cc_max[k] - g5_c_cc_min[k] );
      dGdz_c[m] = 2.0*g5_deriv_c_cc[ii][2][k]/( g5_c_cc_max[k] - g5_c_cc_min[k] );
    } // for k


  /*--------------
     force loop
  --------------*/
  if (nlayer_c == 4) {
    int w1,w2,w2_pre,w3,w3_pre;
    int h1,h2,h2_pre;
    double dummy,dummy1,dummy2;
    w2_pre = ( 1 + layer_c[0] ) * layer_c[1];
    w3_pre = w2_pre + ( 1 + layer_c[1] ) * layer_c[2];
    h2_pre = 1 + layer_c[1];
    for (int l0 = 1; l0 <= layer_c[0]; l0++){
    for (int l1 = 1; l1 <= layer_c[1]; l1++){
      w1 = ( 1 + layer_c[0] )*( l1 - 1 ) + l0;
      h1 = l1;
      dummy1 = weight_c[w1] * ( 1.0 - pow(hidden_c[h1], 2.0) );
    for (int l2 = 1; l2 <= layer_c[2]; l2++){
      w2 = w2_pre + ( 1 + layer_c[1] )*( l2 - 1 ) + l1;
      h2 = h2_pre + l2;
      dummy2 = weight_c[w2] * ( 1.0 - pow(hidden_c[h2], 2.0) );
    for (int l3 = 1; l3 <= layer_c[3]; l3++){
      w3 = w3_pre + ( 1 + layer_c[2] )*( l3 - 1 ) + l2;
      dummy = weight_c[w3] * dummy1 * dummy2;
      //dummy =   weight_c[w1] * ( 1.0 - pow(hidden_c[h1], 2.0) ) 
      //        * weight_c[w2] * ( 1.0 - pow(hidden_c[h2], 2.0) )
      //        * weight_c[w3];
      fx += -dGdx_c[l0-1] * dummy;
      fy += -dGdy_c[l0-1] * dummy;
      fz += -dGdz_c[l0-1] * dummy;
    } // for l4
    } // for l3
    } // for l2
    } // for l1

  } // if nlayer_c == 4


}

/* ---------------------------------------------------------------------- */

void PairNNP3v::arraynnp()
{
//printf("\n arraynnp daze !!");
  int n;

  memory->create(hidden_a,nhidden_a,"pair:nnp3:hidden_a");
  for (n = 0; n < nhidden_a; n++) hidden_a[n] = 0.0;
  memory->create(hidden_b,nhidden_b,"pair:nnp3:hidden_b");
  for (n = 0; n < nhidden_b; n++) hidden_b[n] = 0.0;
  memory->create(hidden_c,nhidden_c,"pair:nnp3:hidden_c");
  for (n = 0; n < nhidden_c; n++) hidden_c[n] = 0.0;

  // input
  memory->create(input_a,layer_a[0]+1,"pair:nnp3:input_a");
  memory->create(input_b,layer_b[0]+1,"pair:nnp3:input_b");
  memory->create(input_c,layer_c[0]+1,"pair:nnp3:input_c");

  // dGdx, dGdy, dGdz
  memory->create(dGdx_a,layer_a[0],"pair:nnp3:dGdx_a");
  memory->create(dGdy_a,layer_a[0],"pair:nnp3:dGdy_a");
  memory->create(dGdz_a,layer_a[0],"pair:nnp3:dGdz_a");

  memory->create(dGdx_b,layer_b[0],"pair:nnp3:dGdx_b");
  memory->create(dGdy_b,layer_b[0],"pair:nnp3:dGdy_b");
  memory->create(dGdz_b,layer_b[0],"pair:nnp3:dGdz_b");

  memory->create(dGdx_c,layer_c[0],"pair:nnp3:dGdx_c");
  memory->create(dGdy_c,layer_c[0],"pair:nnp3:dGdy_c");
  memory->create(dGdz_c,layer_c[0],"pair:nnp3:dGdz_c");

  memory->create(g2_a_a,g_a[0],"pair:nnp3:g2_a_a");
  for (n = 0; n < g_a[0]; n++) g2_a_a[n] = 0.0;
  memory->create(g2_a_b,g_a[1],"pair:nnp3:g2_a_b");
  for (n = 0; n < g_a[1]; n++) g2_a_b[n] = 0.0;
  memory->create(g2_a_c,g_a[2],"pair:nnp3:g2_a_c");
  for (n = 0; n < g_a[2]; n++) g2_a_c[n] = 0.0;
  memory->create(g5_a_aa,g_a[3],"pair:nnp3:g5_a_aa");
  for (n = 0; n < g_a[3]; n++) g5_a_aa[n] = 0.0;
  memory->create(g5_a_ab,g_a[4],"pair:nnp3:g5_a_ab");
  for (n = 0; n < g_a[4]; n++) g5_a_ab[n] = 0.0;
  memory->create(g5_a_ac,g_a[5],"pair:nnp3:g5_a_ac");
  for (n = 0; n < g_a[5]; n++) g5_a_ac[n] = 0.0;
  memory->create(g5_a_bb,g_a[6],"pair:nnp3:g5_a_bb");
  for (n = 0; n < g_a[6]; n++) g5_a_bb[n] = 0.0;
  memory->create(g5_a_bc,g_a[7],"pair:nnp3:g5_a_bc");
  for (n = 0; n < g_a[7]; n++) g5_a_bc[n] = 0.0;
  memory->create(g5_a_cc,g_a[8],"pair:nnp3:g5_a_cc");
  for (n = 0; n < g_a[8]; n++) g5_a_cc[n] = 0.0;

  memory->create(g2_b_a,g_b[0],"pair:nnp3:g2_b_a");
  for (n = 0; n < g_b[0]; n++) g2_b_a[n] = 0.0;
  memory->create(g2_b_b,g_b[1],"pair:nnp3:g2_b_b");
  for (n = 0; n < g_b[1]; n++) g2_b_b[n] = 0.0;
  memory->create(g2_b_c,g_b[2],"pair:nnp3:g2_b_c");
  for (n = 0; n < g_b[2]; n++) g2_b_c[n] = 0.0;
  memory->create(g5_b_aa,g_b[3],"pair:nnp3:g5_b_aa");
  for (n = 0; n < g_b[3]; n++) g5_b_aa[n] = 0.0;
  memory->create(g5_b_ab,g_b[4],"pair:nnp3:g5_b_ab");
  for (n = 0; n < g_b[4]; n++) g5_b_ab[n] = 0.0;
  memory->create(g5_b_ac,g_b[5],"pair:nnp3:g5_b_ac");
  for (n = 0; n < g_b[5]; n++) g5_b_ac[n] = 0.0;
  memory->create(g5_b_bb,g_b[6],"pair:nnp3:g5_b_bb");
  for (n = 0; n < g_b[6]; n++) g5_b_bb[n] = 0.0;
  memory->create(g5_b_bc,g_b[7],"pair:nnp3:g5_b_bc");
  for (n = 0; n < g_b[7]; n++) g5_b_bc[n] = 0.0;
  memory->create(g5_b_cc,g_b[8],"pair:nnp3:g5_b_cc");
  for (n = 0; n < g_b[8]; n++) g5_b_cc[n] = 0.0;

  memory->create(g2_c_a,g_c[0],"pair:nnp3:g2_c_a");
  for (n = 0; n < g_c[0]; n++) g2_c_a[n] = 0.0;
  memory->create(g2_c_b,g_c[1],"pair:nnp3:g2_c_b");
  for (n = 0; n < g_c[1]; n++) g2_c_b[n] = 0.0;
  memory->create(g2_c_c,g_c[2],"pair:nnp3:g2_c_c");
  for (n = 0; n < g_c[2]; n++) g2_c_c[n] = 0.0;
  memory->create(g5_c_aa,g_c[3],"pair:nnp3:g5_c_aa");
  for (n = 0; n < g_c[3]; n++) g5_c_aa[n] = 0.0;
  memory->create(g5_c_ab,g_c[4],"pair:nnp3:g5_c_ab");
  for (n = 0; n < g_c[4]; n++) g5_c_ab[n] = 0.0;
  memory->create(g5_c_ac,g_c[5],"pair:nnp3:g5_c_ac");
  for (n = 0; n < g_c[5]; n++) g5_c_ac[n] = 0.0;
  memory->create(g5_c_bb,g_c[6],"pair:nnp3:g5_c_bb");
  for (n = 0; n < g_c[6]; n++) g5_c_bb[n] = 0.0;
  memory->create(g5_c_bc,g_c[7],"pair:nnp3:g5_c_bc");
  for (n = 0; n < g_c[7]; n++) g5_c_bc[n] = 0.0;
  memory->create(g5_c_cc,g_c[8],"pair:nnp3:g5_c_cc");
  for (n = 0; n < g_c[8]; n++) g5_c_cc[n] = 0.0;



}
