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


#include "fix_nomem.h"
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

FixNomem::FixNomem(LAMMPS *lmp, int narg, char **arg) :
 Fix(lmp, narg, arg),pol_tmp(NULL),pol(NULL)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix ave/atom", error);

  f_m = utils::numeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);

  time_depend = 1;

  // expand args if any have wildcard character "*"
  // this can reset nvalues
  


}
FixNomem::~FixNomem()
{

memory ->destroy(pol_tmp);

memory ->destroy(pol);

}


int FixNomem::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}




 

 
void FixNomem::post_force(int /*vflag*/)
{

   int nlocal = atom->nlocal;
   tagint *tag = atom->tag;
   tagint *molecule = atom->molecule;
   double **x = atom->x;
   double **f = atom->f;
   int *molatom = atom->molatom;
   int id_mol = 0 ;
   imageint *image = atom->image;
   int *mask = atom->mask;

 


    for (int i = 0; i < atom->nlocal; i++){
    
    // molecule id starts from 1. max(id_mol) equals to the number of molecules in the system
    id_mol = std::max(id_mol, molecule[i]);

  }
    
  MPI_Allreduce(&id_mol, &n_mol, 1, MPI_LMP_TAGINT, MPI_MAX, world); 






  memory->destroy(pol_tmp);
  memory->create(pol_tmp , n_mol , 2, "fix_nomem::pol_tmp");
  

  memset(*pol_tmp, 0, sizeof(double) * (n_mol ) * 2);

   //calculate the average molecular volume and polarity over all atoms in the system
   for ( int l=0 ; l< nlocal;l++){
   
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
   

/*for(int i =0 ; i < n_mol ; i++){
   fprintf( screen,"Molecule poltmp  %d : %f %f\n",i+1,pol_tmp[i][0],pol_tmp[i][1]);
 }*/
 
 memory->destroy(pol);
 memory->create(pol, n_mol, 2, "fix_nomem::pol");
 MPI_Allreduce(*pol_tmp,*pol,n_mol * 2, MPI_DOUBLE,MPI_SUM,world);
 


  /*for(int i =0 ; i < n_mol ; i++){
  fprintf( screen,"Molecule pol %d : %f %f\n",i+1,pol[i][0],pol[i][1]);
 }*/
  
for(int i =0 ; i < n_mol ; i++){
  polar= sqrt(pol[i][0]*pol[i][0] + pol[i][1]*pol[i][1] ) ;
  pol[i][0]/=polar ;
  pol[i][1]/=polar ;
} 
 
 /*for(int i =0 ; i < n_mol ; i++){
   printf("Molecule pol %d : %f %f\n",i+1,pol[i][0],pol[i][1]);
 }*/

 /* for(int i =0 ; i < n_mol ; i++){
   fprintf(screen ,"Moleculeforce sgn  %d : %d\n",i+1,forcesign[i]);
 } */

 

 for(int i=0 ; i < nlocal ; i++ ){
  if (mask[i] & groupbit){
    
   
        f[i][0]+=f_m * pol[molecule[i]-1][0];
        f[i][1]+= f_m * pol[molecule[i]-1][1];
    } 

 }

    
        

 


 }
   
  
        
        
          

       
       
        

       


       
      
       
 



  
