/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nomem,FixNomem);
// clang-format on
#else

#ifndef LMP_FIX_NOMEM_H
#define LMP_FIX_NOMEM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNomem: public Fix {
 public:
  FixNomem(class LAMMPS *, int, char **);
 ~FixNomem() override;
  int setmask();
  virtual void post_force(int);
 private:
  
  
  int n_mol; 
  int nrepeat;
  double f_m ;

  double **pol_tmp , **pol ;
  double polar ;
  

};

}    // namespace LAMMPS_NS

#endif
#endif
