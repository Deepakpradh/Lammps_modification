// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "fix_molactive.h"
#include "domain.h"
#include "region.h"
#include "arg_info.h"
#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "group.h"
#include "modify.h"
#include "update.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */

FixMolactive::FixMolactive(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg),molvel(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix ave/atom", error);
  //if (atom->molecular ==1)error->all(FLERR ,"molecular system cannot use fix molactive");

  f_m = utils::numeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);

  time_depend = 1;

  // expand args if any have wildcard character "*"
  // this can reset nvalues
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;

  int id_mol = 0 ;
for (int i = 0; i < atom->nlocal; i++) {
  if (mask[i] & groupbit)
{
  
    // molecule id starts from 1. max(id_mol) equals to the number of molecules in the system
   id_mol = std::max(id_mol, molecule[i]);

  
}
}


MPI_Allreduce(&id_mol, &n_mol, 1, MPI_LMP_TAGINT, MPI_MAX, world); 







memory->grow(molvel, n_mol , 2 *nrepeat, "fix_molactive::v_mol");


for (int i = 0; i < n_mol; i++){
    for (int j = 0; j < 2*nrepeat; j++) {
      molvel[i][j] = 0.0;

  }
}


fprintf( screen,"constructor is called \n");

}
FixMolactive::~FixMolactive()
{
memory ->destroy(molvel);
}



int FixMolactive::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}



void FixMolactive::post_force(int /*vflag*/)
{
   



   int nlocal = atom->nlocal;
   tagint *tag = atom->tag;
   tagint *molecule = atom->molecule;
   double **x = atom->x;
   double **f = atom->f;
   double **v = atom->v;
   int *molatom = atom->molatom;
    int *mask = atom->mask;

   
   double T[2];
   imageint *image = atom->image;

  memory->create(v_mol, n_mol , 2, "fix_molactive::v_mol");
  memory->create(v_mol_tmp, n_mol , 2, "fix_molactive::v_mol_tmp");
  memset(*v_mol_tmp, 0, sizeof(double) * (n_mol ) * 2);

  for(int i=0 ;i< nlocal;++i){
    if (mask[i] & groupbit){
  
    v_mol_tmp[molecule[i]-1][0]+=v[i][0];
    v_mol_tmp[molecule[i]-1][1]+=v[i][1];
  
  }
  }


MPI_Allreduce(*v_mol_tmp, *v_mol, (n_mol) * 2, MPI_DOUBLE, MPI_SUM, world);

bigint ntimestep = update->ntimestep;
int ti = (update->ntimestep) % nrepeat;





for(int j=0 ; j< n_mol; ++j){
    molvel[j][ti]=v_mol[j][0];
    molvel[j][nrepeat+ti]=v_mol[j][1];
}



if(ntimestep<100) return ;
else{

  
  memory->create(pol_tmp , n_mol , 2, "fix_motile::avg_molv");
  

  memset(*pol_tmp, 0, sizeof(double) * (n_mol ) * 2);

   //calculate the average molecular volume and polarity over all atoms in the system
   for ( int l=0 ; l< nlocal;l++){
        if (mask[l] & groupbit){
   
          double unwrap[3];  
         if (molatom[l]==0){
          
      
          domain->unmap(x[l],image[l],unwrap);
          
         /* fprintf( screen,"unwrap pol 1 %d : %f %f\n",molatom[l], unwrap[0],unwrap[1]);*/
          
          pol_tmp[molecule[l]-1][0] +=unwrap[0];
          pol_tmp[molecule[l]-1][1] +=unwrap[1];
         }
         if(molatom[l]==19){
         
         
         domain->unmap(x[l],image[l],unwrap);
         
        /* fprintf( screen,"unwrap pol 20 %d : %f %f\n",molatom[l], unwrap[0],unwrap[1]); */
         pol_tmp[molecule[l]-1][0]-=unwrap[0];
         pol_tmp[molecule[l]-1][1]-=unwrap[1];

      } 

 }
   }

   

 
 memory->create(pol, n_mol, 2, "fix_motile::avg_molv_tmp");
 MPI_Allreduce(*pol_tmp,*pol,n_mol * 2, MPI_DOUBLE,MPI_SUM,world);
 memory->destroy(pol_tmp);
for(int i =0 ; i < n_mol ; i++){
  polar= sqrt(pol[i][0]*pol[i][0] + pol[i][1]*pol[i][1] ) ;
  pol[i][0]/=polar ;
  pol[i][1]/=polar ;
}




memory->create(v_molavg, n_mol , 2, "fix_tgnh_drude::v_mol");
memset(*v_molavg, 0, sizeof(double) * (n_mol ) * 2);

for(int k=0 ;k< n_mol ;++k){
    for(int l=0 ;l <nrepeat ; ++l){
        v_molavg[k][0]+=molvel[k][l]/nrepeat ;
        
    }

    for(int m =nrepeat ;m < 2 *nrepeat ; ++m){
        v_molavg[k][1]+=molvel[k][m]/nrepeat ;
        
    }
}
/*
fprintf(screen,"Timestep %ld  ",ntimestep) ;
for(int i =0 ; i < n_mol ; i++){
   fprintf( screen,"Molecule vel molactive  %d : %f %f\n",i ,v_molavg[i][0],v_molavg[i][1]);
 }
 */







memory->create(forcesign,n_mol, "fix_motile::forcesign");

for(int k =0 ;  k < n_mol ; k++) {
   if(v_molavg[k][0] * pol[k][0] + v_molavg[k][1]*pol[k][1] > 0 ){
    forcesign[k] = 1 ;
   }
   else{
     forcesign[k]= -1;
   }
   }

 for(int i=0 ; i < nlocal ; i++ ){
  if (mask[i] & groupbit){
    
   
        f[i][0]+=f_m * pol[molecule[i]-1][0] * forcesign[molecule[i]-1];
        f[i][1]+= f_m * pol[molecule[i]-1][1] *forcesign[molecule[i]-1];
    } 
  

 }

memory->destroy(forcesign); 
memory->destroy(pol); 
memory->destroy(v_mol_tmp);
memory->destroy(v_mol);
memory->destroy(v_molavg);




}


}





