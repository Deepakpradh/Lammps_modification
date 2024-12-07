/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory. **/


#include "fix_area.h"
#include <mpi.h>
#include <cstring>

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom_vec.h"
#include "molecule.h"
#include "neighbor.h"
#include "group.h"
#include <cmath>
#include <map>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;


FixArea::FixArea(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),curarea(NULL),currareaproc(NULL)

{
  if (strcmp(style,"area") != 0 && narg < 5)
    error->all(FLERR,"Illegal fix area command: not enough args");
  A0 = utils::numeric(FLERR,arg[3],false,lmp);
  Xi =utils::numeric(FLERR,arg[4],false,lmp);
  n =utils::numeric(FLERR,arg[5],false,lmp); //number of atoms per molecule 

 

 
 
 


}

FixArea::~FixArea()
{
  memory->destroy(currareaproc);
  memory->destroy(curarea) ;
  
}




int FixArea::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  
  return mask;
}

void FixArea::post_force(int /*vflag*/)
{

  double fx, fy ;

  int xbox,ybox,zbox;

  
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int id_mol = 0 ;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double *h = domain->h;
  double **x = atom->x;
  double **f = atom->f;
  int *molatom = atom->molatom;
  int molecular = atom->molecular;
  imageint *image = atom->image;
  int *type = atom->type;

  int nlocal= atom->nlocal;
  int *mask = atom->mask;
  tagint **bond_atom = atom->bond_atom;
  tagint tagprev;
  int *num_bond = atom->num_bond;
  int *molindex = atom->molindex;
  Molecule **onemols = atom->avec->onemols;
  int atom1, atom2,atom3 , yatom ,atom4,atom5;
  int jmol,jatom,nb,imol;  



  
 /*  fprintf(screen, "TIMESTEP= %ld\n", update->ntimestep);
  fprintf(screen,"molecular = %d \n", molecular);*/


  for (int i = 0; i < atom->nlocal; i++) {
    // molecule id starts from 1. max(id_mol) equals to the number of molecules in the system
    id_mol = std::max(id_mol, molecule[i]);

  }
  MPI_Allreduce(&id_mol, &n_mol, 1, MPI_LMP_TAGINT, MPI_MAX, world);





// zero out local per-mol areas

memory->destroy(currareaproc);
memory->create(currareaproc,n_mol,"FixArea::currareaproc");


for (int i = 0; i < n_mol; i++) {
currareaproc[i] = 0.0;
}


double unwrap3[3];
double unwrap4[3];



// compute current per-mol areas
for (int i = 0; i < nlocal; i++) {		
		if (mask[i] & groupbit){
    atom1=i ;

    if (molecular == 1) {
	nb = num_bond[atom1];
    }
				
	
	else {
		  if (molindex[atom1] < 0) continue;
		   jmol = molindex[atom1];
		   jatom = molatom[atom1];
       tagprev = tag[atom1] - jatom - 1;

			 nb = onemols[jmol]->num_bond[jatom];
			}

      
    for (int j = 0; j < nb; j++) {
				if (molecular == 1) {
					atom2 = atom->map(bond_atom[atom1][j]);

        }
        else  {
              
			      atom2 = atom->map(onemols[jmol]->bond_atom[jatom][j]+tagprev);
           /*fprintf(screen, "jatj =  %d  \n",onemols[jmol]->bond_atom[jatom][j]);*/
        }

    

      /*if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;	*/
      
      domain->unmap(x[atom1],image[atom1] ,unwrap3);
      domain->unmap(x[atom2] ,image[atom2] ,unwrap4 );

    /** fprintf(screen, "atom1=  %f %f  \n",unwrap3[0],unwrap3[1] );
      fprintf(screen, "atom2=  %f %f  \n",unwrap4[0],unwrap4[1] );  */




// get the molecule index of the atom // to see which cell the atom is connected to
int index = molecule[atom1]-1;
if (index < 0) continue;

// calculate the cross products locally over procs




currareaproc[index] += ((unwrap3[0])*(unwrap4[1]) - (unwrap4[0])*(unwrap3[1]))/2.0;



}

    }
}   
  


// sum over the procs and distribute to each proc
memory->destroy(curarea);
memory->create(curarea,n_mol,"FixArea::curarea");

MPI_Allreduce(currareaproc,curarea,n_mol,MPI_DOUBLE,MPI_SUM,world);
MPI_Barrier(world);

/*fprintf(screen, "TIMESTEP= %d\n", update->ntimestep);

for (int i = 0; i < n_mol; i++) {
     fprintf(screen, "area =  %d %f \n",i, curarea[i]);

} */

/*fprintf(screen, "TIMESTEP= %d\n", update->ntimestep);*/

   double unwrap1[3];
    double unwrap2[3];
    for (int i = 0; i < nlocal; i++) {		
		if (mask[i] & groupbit) {
      

      for (int j = 0; j < nb; j++){
      
      

      




       atom3 =i ;
       jmol = molindex[atom3];
		   jatom = molatom[atom3];
       tagprev = tag[atom3] - jatom - 1;
       /*fprintf(screen, "atom_no= %d \n", tag[atom3]);*/
       atom4 =atom->map(onemols[jmol]->angle_atom1[jatom][j]+tagprev) ;
       atom5 =atom->map(onemols[jmol]->angle_atom3[jatom][j]+tagprev) ;


    

     /*fprintf(screen, "atom4 = %d \n", tag[atom4] );
       fprintf(screen,"atom5 = %d \n" , tag[atom5]);*/

       domain->unmap(x[atom4],image[atom4] ,unwrap1);
      /* fprintf(screen,"atom4 unwrap = %d  %f %f \n" ,atom4 , unwrap1[0],unwrap1[1]);*/


       domain->unmap(x[atom5] ,image[atom5] ,unwrap2);
       /*fprintf(screen,"atom5 unwrap  = %d  %f %f \n" ,atom5 , unwrap2[0],unwrap2[1]); */
       fx=  0.5*Xi *(1-(curarea[molecule[i]-1]/A0))*(unwrap2[1]-unwrap1[1]);
       fy=  0.5* Xi *(1-(curarea[molecule[i]-1]/A0))*(unwrap1[0]-unwrap2[0]);

        


      
    

 
 
  
   
    
    



     
    
    
    

   /* fprintf(screen, "fx= %f \n", fx);
    fprintf(screen, "fy= %f \n", fy); */

    f[i][0] += fx;  
    f[i][1] += fy; 


      }

    }

    }

}


/* 
A =Σxi(yi+1 - yi-1)
*/

