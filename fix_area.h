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
#ifdef FIX_CLASS

FixStyle(area,FixArea)

#else

#ifndef LMP_FIX_AREA_H
#define LMP_FIX_AREA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixArea : public Fix {
 public:
  FixArea(class LAMMPS *, int, char **);
  ~FixArea() override;
  int setmask();
  virtual void post_force(int);

 protected:
 
  double A0;
  double Xi;
  int n ;
  int n_mol; 
  double *currareaproc;
  double *curarea;
};

}

#endif
#endif