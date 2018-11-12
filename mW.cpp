/* This program performs Grand Canonical Monte Carlo simulations of monatomic
 * water (mW). It also has a Lennard Jones potential functiomn, which can be 
 * used to test the program. The program uses a cell structure, where the length
 * of a cell is the cut-off radius of interaction, a*sigma. 
 * The parameters for the program are epsilon, the energy term in the potential
 * and mu, the chemical potential. Epsilon is a reduced unit here, therefore the
 * user should take the epsilon stated in the Molinero and Moore 2009 paper and
 * divide by kBT, where T is the required temperature of simulation. Some useful
 * ones are:
 * T=298K:  epsilon = 10.451
 * T=450K:	epsilon = 6.921
 * T=910K:  epislon = 3.422
 * T=0.95*Tc=864.5K	epsilon = 3.6026
 * T=0.85*Tc=773.5K 	epsilon = 4.0264
 * T = 915K		epsilon = 3.4037

 * T = 750K		epsilon = 4.1526
 * T=800K 		epsilon = 3.8930
 * T=850K 		epsilon = 3.6640
 * T=890K		epsilon = 3.4993
 * T=840K		epsilon = 3.7077
 * For LJ: epsilon = 3.2747 mu = -2.7... a = 2.5 for critical point
 * Useful chemical potentials to follow.
 * Coordinates of molecules are in terms of the length of the molecule, however
 * from the equations we can cancel sigma, therefore when performing calculations
 * we can essential set sigma = 1 and only multiply by a. For more information
 * see the supporting documentation (to follow). 
 *
 * This program was written by Mary Coe, University of Bristol. 
 * E-mail: m.k.coe@bristol.ac.uk*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <deque>
#include <string>
#include <sstream>
#include <iomanip>
extern "C" {
#include "ranvec.h"
}
#include <random>
#include <ctime>
#include <cstdlib>

using namespace std;

/* These parameters allow the user to change how the simulation runs fundamentally. This is the 
 * only point at which the user should have to edit the main bulk of the program. Note that currently,
 * SURFACE and SPHERICAL do nothing within the program. These are place holders for future development. */
 
//#define LJ				// Uncomment to run LJ potential to test GCMC.
//#define SURFACE			// Uncomment to run with surface along one side of simulation box. Default is planar surface. 
//#define SPHERICAL		// Uncomment to run with spherical surface. 
//#define TRANSLATION		// Uncomment to perform translational moves in addition to insertion/deletion.
					
#define MULTI
//#define MULTI_HAND


#define NRANDS 10000	 // These are the parameters needed for both mW and LJ simulations
#define mu -5.3668
#define PI 3.14159

#define LJ_epsilon 3.2747  // The LJ parameters are defined in case of using LJ potential
#define LJ_a 2.5
#define LJ_sigma 1.0

#define epsilon 3.6026 // These parameters are needed only for mW simulations
#define sigma 2.3925
#define A 7.049556277
#define B 0.6022245584
#define gamma 1.2
#define a 1.8
#define aa 3.24
#define lambda 23.15 
#define MW 18.015
#define NA 6.022
#define max_N 200



class water_molecule {
	
	/*Class structure for the 'water molecules' which are modelled as spheres*/
	
	/* Describes fraction of insertion if using staged insertion and deletion.
	 * Also keeps track of the energy of the particle for unbiased insertion/deletion
	 * to prevent constant recalculating. */
	#ifdef MULTI
	private: double energy; int stage;
	#endif
	
	public: int id; // index in global molecule list
			double x_pos,y_pos,z_pos; // position of molecule within cell
			int cell_id; //index of occupied cell in cell_list
			
			water_molecule() {};
			water_molecule(int ID, double x_pos_in, double y_pos_in, double z_pos_in, int cell_id, bool initial);
};

class cell {
	
	/*Class structure for a cell*/
	
	public: int id; //id gives place of cell in cell_list
			int occupancy[50]; //ids of occupants of the cell
			int num_occupants; //number of occupants in cell
			int i,j,k; //position of cell in real space w.r.t other cells (imagine 3D grid)
			int neighbours[27][4]; //Gives information of neighbouring cells within cell list
			
			cell() {};
			cell(int ID, int i_in, int j_in, int k_in, int num_occupants_in); //constructor for cell
};

class sim_box {
	
	/*Class structure for the simulation box. */
	
	private: double volume, running_energy;
			 double *random_numbers=new double[NRANDS];
			 int num_cells,num_cells_row, mca, mcr, num_mols, seed, current_steps, num_steps, num_steps_eqb, num_gr, save_frq, rnd_check;
			 string file_path, in_file_path;
			#ifdef MULTI
				double Collection[max_N][3], Totals_Collect[max_N], Transition[max_N][3], Weights[max_N];
				//double total; // Keeps track of the total unbiased probability in the transition matrix.
				int max_mols; 
			#endif

			#ifdef MULTI_HAND
			double Weights[max_N];
			#endif							
			
	
	public: 	double return_volume(void) {return volume;}
			
			/* Returns cell information. */
			int return_num_cells_row(void) {return num_cells_row;}
			int return_num_cells(void) {return num_cells;}
			
			/* Returns information about steps and saving frequencies*/
			int return_seed(void) {return seed;}
			int return_current_steps(void) {return current_steps;}
			int return_steps_eqb(void) {return num_steps_eqb;}
			int return_steps(void) {return num_steps;}
			int return_num_gr(void) {return num_gr;}
			int return_save_frq(void) {return save_frq;}
			string return_file_path(void) {return file_path;}
			string return_in_file_path(void) {return in_file_path;}
			
			/* Provides details on random numbers in system */
			int return_rnd_check(void) {return rnd_check;}
			void check_rnd(void) {
				if (NRANDS-rnd_check<2) {
					vector_random_generator(NRANDS, random_numbers); 
					rnd_check = 0; 
				}
			}
			double get_rnd(void) {
				double rnd;
				rnd = random_numbers[rnd_check];
				rnd_check = rnd_check+1; 
				check_rnd();
				return rnd;
			}
			
			/* Provides details of number of accepted and rejected moves. */
			void inc_mca(void) {mca=mca+1;} 
			void inc_mcr(void) {mcr=mcr+1;}
			int return_mca(void) {return mca;}
			int return_mcr(void) {return mcr;}
			
			/* Provides details of the running energy of the system. */
			void inc_running_energy(double energy) {running_energy=running_energy+energy;}
			void dec_running_energy(double energy) {running_energy=running_energy-energy;}
			double return_running_energy(void) {return running_energy;}
			void set_running_energy(double energy) {running_energy = energy;}
			
			/* Provides details of number of molecules in the system */
			void inc_mols(void) {num_mols+=+1;}
			void dec_mols(void) {num_mols-=1;}
			int return_mols(void) {return num_mols;}
			void set_mols(int mols) {num_mols=mols;}
			
		
			/*TMMC matrices*/
			#ifdef MULTI
			// Collection matrix stores all the probabilities of transitions
			void update_Collection(int N, int Na, double prob, bool acc) { 
				//Uses N as current num of mols. Na:  0 = N-1, 1 = N, 2 = N+1
				if(acc) {Collection[N][Na]+=prob;  Totals_Collect[N]+=prob;}
				else {Collection[N][Na]+=(1-prob); Totals_Collect[N]+=(1-prob);}
			}
			void update_Transition(void) {
				//Uses convention j=N; 0 = N-1, 0 = N, 1 = N+1
				for(int j = 0;j<max_mols; j++) {
					for(int i=0; i<3; i++) {
					 	Transition[j][i] = (Collection[j][i]+1.0)/(Totals_Collect[j]+1.0);
						if(isnan(Transition[j][i])) cout<<"j: "<<j<<" i: "<<i<<" C: "<<Collection[j][i]<<" total: "<<Totals_Collect[j]<<" max mols "<<max_mols<<endl;
					}
				}
			}
			void update_weights(void) {
				double total=0.0;
				cout<<"max mols is "<<max_mols<<endl;
				Weights[0] = 1.0;
				total = 1.0;
				for(int i=1;i<max_N;i++) {
					
					Weights[i] = Weights[i-1]*(Transition[i][2]/Transition[i+1][0]);
					//Weights[i] = -1.0*log(P_N);
					//cout<<i<<"->"<<i+1<<": "<<Transition[i][2]<<" "<<i+1<<"->"<<i<<": "<<Transition[i+1][0]<<" P: "<<P<<" P_N: "<<P_N<<" weight: "<<Weights[i]<<" total "<<total<<endl;
					//total+=Weights[i];
					//cout<<"TOTAL IS "<<total<<" WEIGHTS ARE "<<Weights[i]<<endl;
				}
				//for(int i=max_mols;i<max_N;i++) Weights[i] = Weights[max_mols-1];
				//cout<<"TOTAL IS "<<total<<" WEIGHTS ARE "<<Weights[i]<<endl;
				//for(int i=0;i<max_N;i++) Weights[i]/=total;	
			}
			//void update_hits(void) {Hits[num_mols]+=1;cout<<"MOLS IS "<<num_mols<<endl;}
			void write_weights(void) {
				int i, j;
				ofstream hist; string out_name;

				out_name = file_path + "weights.txt";
				hist.open(out_name, ios::out);
				hist<<"Collection Matrix"<<endl;
				for(i = 0;i<max_N; i++) {
					for(j=0; j<3;j++) {
					 	hist<<Collection[i][j]<<" ";
					}
					hist<<endl;
				}
				hist<<endl;
				hist<<"Transition Matrix"<<endl;
				for(i = 0;i<max_N; i++) {
					for(j=0; j<3;j++) {
					 	hist<<Transition[i][j]<<" ";
					}
					hist<<endl;
				}
				hist<<endl;
				hist<<"Weights"<<endl;
				for(i = 0;i<max_N; i++) hist<<Weights[i]<<" "<<endl;
				hist.close();
			}
			double return_weight(int N) {return Weights[N];}
			
			void inc_max_mols(void) {max_mols++;}
			int return_max_mols(void) {return max_mols;}
			#endif
			
			#ifdef MULTI_HAND

			double return_weight(int num_mols) {return Weights[num_mols];}
			#endif

	sim_box() {};
	sim_box(int num_cells_row_in, int mols, int mca_in, int mcr_in, int current_in, int steps_eqb, int steps, int gr, int frq, int seed_in, string file_path_in, string load_file_path, double weights_in[]);
};

/*These constructors set up the classes within the simulation. Sim_box is called only once, and cell is
 * called only in initialisation of the system. */
 
water_molecule::water_molecule(int ID, double x_pos_in, double y_pos_in, double z_pos_in, int cell_id_in, bool initial) {

	id = ID; x_pos = x_pos_in; y_pos=y_pos_in; z_pos = z_pos_in; cell_id = cell_id_in;
	
}

cell::cell(int ID, int i_in, int j_in, int k_in, int num_occupants_in) {

	id = ID; i = i_in; j=j_in;k=k_in; num_occupants = num_occupants_in;

}

sim_box::sim_box(int num_cells_row_in, int mols, int mca_in, int mcr_in, int current_in, int steps_eqb, int steps, int gr, int frq, int seed_in, string file_path_in, string load_file_path, double weights_in[]) {
	
		volume = pow((num_cells_row_in*a),3); num_cells = pow(num_cells_row_in,3); 
		num_cells_row = num_cells_row_in; mca = mca_in; mcr = mcr_in; num_mols = mols;
		seed = seed_in; current_steps = current_in; num_steps_eqb = steps_eqb;
		num_steps = steps; num_gr = gr; save_frq = frq;
		rnd_check=0;
		file_path = file_path_in;
		in_file_path = load_file_path;
		
		init_vector_random_generator(seed, NRANDS);
		vector_random_generator(NRANDS, random_numbers); 
		
		#ifdef MULTI
		for(int i=0;i<max_N;i++) {
			Weights[i]=1.0; Totals_Collect[i] = 0.0;
			for(int j=0;j<max_N;j++) {
				Collection[j][i]=0; Transition[j][i]=1; 
			}
		}
		max_mols = mols; 
		#endif
		
		#ifdef MULTI_HAND
		for(int i=0;i<max_N;i++) {
			Weights[i] = weights_in[i];
			cout<<"Weight: "<<Weights[i]<<" and in "<<weights_in[i]<<endl;
		}
		#endif
}

water_molecule Find_New_Pos(double dis[], cell cell_list[], water_molecule mol, sim_box sim_env) {
	
	int cell_id, cell_dis[3] = {0,0,0},i, num_cells_row = sim_env.return_num_cells_row(), boundary = num_cells_row-1;
	double half_box = num_cells_row/2.;
		
	if(dis[0]>0) {
		if(mol.x_pos + dis[0]<1) {mol.x_pos+=dis[0];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[0]-=(1-mol.x_pos); cell_dis[0]+=1;
			for(i=0;i<=half_box;i++) {
				if (dis[0]>=1) {cell_dis[0]+=1;dis[0]-=1;}
				else mol.x_pos = dis[0];
			}
		}
	}
	else if (dis[0]<0) {
		if(mol.x_pos + dis[0]>0) {mol.x_pos+=dis[0];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[0]+=mol.x_pos; cell_dis[0]-=1;
			for(i=0;i<=half_box;i++) {
				if (dis[0]<=-1) {cell_dis[0]-=1;dis[0]+=1; }
				else mol.x_pos = 1+dis[0];
			}
		}
	}
	
	if(dis[1]>0) {
		if(mol.y_pos + dis[1]<1) {mol.y_pos+=dis[1];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[1]-=(1-mol.y_pos); cell_dis[1]+=1;
			for(i=0;i<=half_box;i++) {
				if (dis[1]>=1) {cell_dis[1]+=1;dis[1]-=1;}
				else mol.y_pos = dis[1];
			}
		}
	}
	else if (dis[1]<0) {
		if(mol.y_pos + dis[1]>0) {mol.y_pos+=dis[1];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[1]+=mol.y_pos; cell_dis[1]-=1;
			for(i=0;i<=half_box;i++) {
				if (dis[1]<=-1) {cell_dis[1]-=1;dis[1]+=1; }
				else mol.y_pos = 1+dis[1];
			}
		}
	}	
	
	if(dis[2]>0) {
		if(mol.z_pos + dis[2]<1) {mol.z_pos+=dis[2];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[2]-=(1-mol.z_pos); cell_dis[2]+=1;
			for(i=0;i<=half_box;i++) {
				if (dis[2]>=1) {cell_dis[2]+=1;dis[2]-=1;}
				else mol.z_pos = dis[2];
			}
		}
	}
	else if (dis[2]<0) {
		if(mol.z_pos + dis[2]>0) {mol.z_pos+=dis[2];} //If the displacement does not take you out of the cell, then simply add
		else {
			dis[2]+=mol.z_pos; cell_dis[2]-=1;
			for(i=0;i<=half_box;i++) {
				if (dis[2]<=-1) {cell_dis[2]-=1;dis[2]+=1; }
				else mol.z_pos = 1+dis[2];
			}
		}
	}
	
	cell_dis[0] = cell_list[mol.cell_id].i+cell_dis[0]; 
	cell_dis[1] = cell_list[mol.cell_id].j+cell_dis[1];
	cell_dis[2] = cell_list[mol.cell_id].k+cell_dis[2];
	
	
	if(cell_dis[0]>boundary) cell_dis[0]-=num_cells_row;
	else if (cell_dis[0]<0) cell_dis[0]+=num_cells_row;
	if(cell_dis[1]>boundary) cell_dis[1]-=num_cells_row;
	else if (cell_dis[1]<0) cell_dis[1]+=num_cells_row;
	if(cell_dis[2]>boundary) cell_dis[2]-=num_cells_row;
	else if (cell_dis[2]<0) cell_dis[2]+=num_cells_row;
	
	
	for(int c=0;c<sim_env.return_num_cells();c++) {
		if(cell_dis[0]==cell_list[c].i and cell_dis[1]==cell_list[c].j and cell_dis[2]==cell_list[c].k) {
			cell_id = c; break;
		}
	}
	
	return mol;
}

double Find_distance(water_molecule mol_a, water_molecule mol_b, int i, int j, int k) {
	
	/*This function returns the squared distance between two molecules in neighbouring cells only. It is not
	 * Uses the formula: r = 1+index*(mol_2.r - mol_1.r) when index is not 0 and r = mol_2.r - mol_1.r when 
	 * index is 0. */
	
	double r2,ri,rj,rk;
	
	if(i==0) ri=mol_b.x_pos - mol_a.x_pos;
	else ri=1+i*(mol_b.x_pos-mol_a.x_pos);
	if(j==0) rj=mol_b.y_pos - mol_a.y_pos;
	else rj=1+j*(mol_b.y_pos-mol_a.y_pos);
	if(k==0) rk=mol_b.z_pos - mol_a.z_pos;
	else rk=1+k*(mol_b.z_pos-mol_a.z_pos);
	
	r2 = ri*ri+rj*rj+rk*rk;
	
	return r2;
	
}

double Find_distance(water_molecule mol_a, water_molecule mol_b, int i, int j, int k, int num_cells_row, bool gr) {
	
	/* This function returns the squared distance between two molecules which are not necessarily in neighbouring cells, 
	 * for example this would be the case for some three body interactions and also for g(r) calculations. The general
	 * form is to work out the number of cells between the two molecules, then work out the shortest distance between
	 * then, before finally adding back in the number of cells between them. Each component is found separately and 
	 * at the end the squared distance is found.*/
	 
	double r2, ri,rj,rk;
	int ci=0,cj=0,ck=0;

	if (i>0) {
		if (i>1) ci = i-1;
		ri = 1+(mol_b.x_pos - mol_a.x_pos);
		if (ci>=(int)floor(num_cells_row/2)-1) {
			if ((1-(mol_b.x_pos - mol_a.x_pos))<ri) ri = 1-(mol_b.x_pos - mol_a.x_pos);
		}
		
	}
	else if (i<0) {
		if(i<-1) ci = -1*(i+1);
		ri = 1-(mol_b.x_pos - mol_a.x_pos);
		
		if (ci>=(int)floor(num_cells_row/2)-1) {
			if ((1+(mol_b.x_pos - mol_a.x_pos))<ri) ri = 1+(mol_b.x_pos - mol_a.x_pos);
		}
	}
	else ri = mol_b.x_pos - mol_a.x_pos;
	
	if (j>0) {
		if (j>1) cj = j-1;
		rj = 1+(mol_b.y_pos - mol_a.y_pos);
		if (cj>=(int)floor(num_cells_row/2)-1) {
			if ((1-(mol_b.y_pos - mol_a.y_pos))<rj) rj = 1-(mol_b.y_pos - mol_a.y_pos);
		}
	}
	else if (j<0) {
		if(j<-1) cj = -1*(j+1);
		rj = 1-(mol_b.y_pos - mol_a.y_pos);
		if (cj>=(int)floor(num_cells_row/2)-1) {
			if ((1+(mol_b.y_pos - mol_a.y_pos))<rj) rj = 1+(mol_b.y_pos - mol_a.y_pos);
		}
	}
	else rj = mol_b.y_pos - mol_a.y_pos;
	if (k>0) {
		if (k>1) ck = k-1;
		rk = 1+(mol_b.z_pos - mol_a.z_pos);
		if (ck>=(int)floor(num_cells_row/2)-1) {
			if ((1-(mol_b.z_pos - mol_a.z_pos))<rk) rk = 1-(mol_b.z_pos - mol_a.z_pos);
		}
	}
	else if (k<0) {
		if(k<-1) ck = -1*(k+1);
		rk = 1-(mol_b.z_pos - mol_a.z_pos);
		if (ck>=(int)floor(num_cells_row/2)-1) {
			if ((1+(mol_b.z_pos - mol_a.z_pos))<rk) rk = 1+(mol_b.z_pos - mol_a.z_pos);
		}
	}
	else rk = mol_b.z_pos - mol_a.z_pos;
	
	r2 = ((ri+ci)*(ci+ri)) + ((cj+rj)*(cj+rj)) + ((ck+rk)*(ck+rk));
	
	return r2;
}

/*This inline function returns the cos(theta) of the angle between three molecules. For use
 * with three body potential. */
inline double Find_angle(double ar2, double br, double cr) { return ((cr*cr + br*br - ar2)/(2*br*cr)); }

double LJPotential (water_molecule mol, water_molecule gmol_list[], cell cell_list[]) {
	
	/*Returns the Lennard Jones contribution of one molecule to the system. Loops over all neighbouring
	 * cells, finding molecules within cut-off, then calculates the 6-12 potential.*/	
	
	double LJ=0, r2,r6,r12;
	int i,j,k,pcell=mol.cell_id;
	deque<water_molecule> occupancy;
	cell mcell;

	for(int n=0;n<27;n++) {
		mcell = cell_list[cell_list[pcell].neighbours[n][0]];
		for (int m=0;m<mcell.num_occupants;m++) {
			if(mol.id == mcell.occupancy[m]) continue;
			r2 = Find_distance(mol, gmol_list[mcell.occupancy[m]],cell_list[mol.cell_id].neighbours[n][1],cell_list[mol.cell_id].neighbours[n][2],cell_list[mol.cell_id].neighbours[n][3]);
			if (r2<=1) {
				r2=1/(r2*LJ_a*LJ_a);
				r6=r2*r2*r2; r12=r6*r6;
				LJ+=((LJ_sigma*r12)-(LJ_sigma*r6));
			}
		}
	} 
	
	return LJ_epsilon*LJ;
}

double SWPotential (water_molecule gmol_list[], cell cell_list[], int num_cells_row, water_molecule mol) {
	/* This function is used to find the complete energy contribution of one molecule to the system. There
	 * are three parts to the contribution:
	 * 1. Pairwise interaction
	 * 2. Three body interactions for the angles around the molecule in question.
	 * 3. Three body interactions from the angles not around the molecule, but which the molecule contributes
	 *    to.
	 * This function is used every move attempt, and NOT for the total energy of the system. Used to improve
	 * efficiency. */
	 
	 double rij2, r4, two_body=0, three_body=0, rik2, rjk2, costheta, neighbour_info[3][50];
	 int neighbour_counter=0, pcell = mol.cell_id, ci, cj, ck, mol_nei;
	 cell nei_cell, mol_cell = cell_list[mol.cell_id];
	 water_molecule nei_mol; 
	 
	 /* First we find the pairwise interaction, and in doing so create a list of neighbouring molecules to
	  * use for the three body interactions. */
	 for (int n=0;n<27;n++) {
		nei_cell = cell_list[cell_list[pcell].neighbours[n][0]];
		for(int m=0;m<nei_cell.num_occupants;m++){
			if(mol.id == nei_cell.occupancy[m]) continue; //If we are comparing the same molecule, continue.
			rij2 = Find_distance(mol,gmol_list[nei_cell.occupancy[m]],cell_list[pcell].neighbours[n][1],cell_list[pcell].neighbours[n][2],cell_list[pcell].neighbours[n][3]);
			if (rij2<=1.) {
				rij2*=aa;
				r4=rij2*rij2;
				rij2 = sqrt(rij2);
				two_body+=(B*(1./r4)-1.)*exp(1./(rij2-a)); //No need to divide by 2 here as adding or removing a molecule will either add the entire energy or remove it
				
				/* Add information to neighbours array*/
				neighbour_info[0][neighbour_counter] = nei_cell.occupancy[m]; //record the neighbouring molecule id
				neighbour_info[1][neighbour_counter] = n; //record the corresponding position in the cell neighbour list
				neighbour_info[2][neighbour_counter] = rij2; //record the distance between the two
				neighbour_counter++;
			}
		}
	}
	
	/* Now we loop back over this neighbour list working out the three body potential for all
	 * sets of triplets.*/
	for(int m=0;m<neighbour_counter;m++) {
		for(int n=m+1;n<neighbour_counter;n++) {
			ci = cell_list[pcell].neighbours[(int)neighbour_info[1][m]][1]-cell_list[pcell].neighbours[(int)neighbour_info[1][n]][1];
			cj = cell_list[pcell].neighbours[(int)neighbour_info[1][m]][2]-cell_list[pcell].neighbours[(int)neighbour_info[1][n]][2];
			ck = cell_list[pcell].neighbours[(int)neighbour_info[1][m]][3]-cell_list[pcell].neighbours[(int)neighbour_info[1][n]][3]; 
			rjk2 = Find_distance(gmol_list[(int)neighbour_info[0][n]],gmol_list[(int)neighbour_info[0][m]],ci,cj,ck, num_cells_row, true);
			rjk2*=aa;
			costheta = Find_angle(rjk2,neighbour_info[2][m],neighbour_info[2][n]);
			
			three_body+=(costheta + (1./3.))*(costheta + (1./3.))*exp((gamma/(neighbour_info[2][m]-a))+(gamma/(neighbour_info[2][n]-a)));
			
		}
	}
	
	/* Now we also need to take into account the energy contribution of the particle
	 * to making tetrahedral angles for its neighbours. */
	for(int wm =0;wm<neighbour_counter;wm++) {
		nei_mol = gmol_list[(int)neighbour_info[0][wm]]; //find the neighbouring molecule
		mol_cell = cell_list[nei_mol.cell_id]; //and its cell
		for(int c=0;c<27;c++) { if(mol_cell.neighbours[c][0] == mol.cell_id) {mol_nei = c; break;}} //Find the neighbouring molecules cell in the cell lists neighbour lists
		for(int n=0;n<27;n++) { //for all of the neighbouring molecules neighbouring cells...
			nei_cell = cell_list[mol_cell.neighbours[n][0]]; //locate the neighbours neighbouring cell...
			for(int m=0;m<nei_cell.num_occupants;m++) {//for all the neighbours neighbouring cell's occupants
				if(mol.id==nei_cell.occupancy[m]) continue; //the angle jij will be 0 and we dont want to add this to the total energy
				if(neighbour_info[0][wm] == nei_cell.occupancy[m]) continue;
				rjk2 = Find_distance(nei_mol, gmol_list[nei_cell.occupancy[m]], mol_cell.neighbours[n][1],mol_cell.neighbours[n][2],mol_cell.neighbours[n][3]); //Should this use gr version?!
				if(rjk2<=1.0) {			
					rjk2=sqrt(rjk2)*a;
					ci = mol_cell.neighbours[n][1]-mol_cell.neighbours[mol_nei][1]; 
					cj = mol_cell.neighbours[n][2]-mol_cell.neighbours[mol_nei][2];
					ck = mol_cell.neighbours[n][3]-mol_cell.neighbours[mol_nei][3]; 
					rik2 = Find_distance(mol,gmol_list[nei_cell.occupancy[m]],ci,cj,ck, num_cells_row, true);
					rik2*=aa;
					costheta = Find_angle(rik2,rjk2,neighbour_info[2][wm]);
					three_body+=(costheta + (1./3.))*(costheta + (1./3.))*exp((gamma/(rjk2-a))+(gamma/(neighbour_info[2][wm]-a)));
				}
			}
		}
	}
	

	return epsilon*(A*two_body + lambda*three_body);
}

double Find_System_Energy_LJ (water_molecule gmol_list[], cell cell_list[], int num_mols) {
	
	/* Returns the system energy when using LJ potential. */
	double LJ=0;
	
	for(int m=0;m<num_mols;m++) LJ+=LJPotential(gmol_list[m], gmol_list,cell_list);
	
	return LJ;
}

double Find_System_Energy_SW (water_molecule gmol_list[], cell cell_list[], int num_cells_row, int num_mols) {
	
	/* Works out the total system energy using the original Stillinger-Weber potential equation (see Molinero and
	 * Moore 2009 paper.*/
	
	double two_body=0, three_body=0, r2,r4, nnr2,ctheta;
	deque<water_molecule> neighbours;
	int  ci=0,cj=0,ck=0, ncell, mcell, neighbour_counter=0, pcell;
	cell rcell;
	water_molecule mol;
	deque<int> cell_neighbour_pos;
	deque<double>  ndistances;
	
	/* Loop over all molecules and find first the two body interaction then the three body interaction. */
	for(int wm=0;wm<num_mols;wm++) {
		mol = gmol_list[wm];
		pcell = mol.cell_id;
		for (int n=0;n<27;n++) {
			rcell = cell_list[cell_list[pcell].neighbours[n][0]];
			for(int m=0;m<rcell.num_occupants;m++){
				if(mol.id == rcell.occupancy[m]) continue;
				r2 = Find_distance(mol,gmol_list[rcell.occupancy[m]],cell_list[pcell].neighbours[n][1],cell_list[pcell].neighbours[n][2],cell_list[pcell].neighbours[n][3]);
				if (r2<=1.0) {
					neighbours.push_back(gmol_list[rcell.occupancy[m]]);
					cell_neighbour_pos.push_back(n);
					ndistances.push_back(sqrt(r2)*a);
					r2*=aa;
					r4=r2*r2;
					two_body+=(B*(1./r4)-1.)*exp(1./((sqrt(r2)-a)));
				
				}
			}
		}
	
		for(int m=0;m<neighbours.size();m++) {
			for(int n=m+1;n<neighbours.size();n++) {
				mcell = cell_neighbour_pos[m];
				ncell = cell_neighbour_pos[n];
			
				ci = cell_list[pcell].neighbours[mcell][1]-cell_list[pcell].neighbours[ncell][1];
				cj = cell_list[pcell].neighbours[mcell][2]-cell_list[pcell].neighbours[ncell][2];
				ck = cell_list[pcell].neighbours[mcell][3]-cell_list[pcell].neighbours[ncell][3]; 
			
				nnr2 = Find_distance(neighbours[n],neighbours[m],ci,cj,ck, num_cells_row, true);
				nnr2*=aa;
				ctheta = Find_angle(nnr2,ndistances[n],ndistances[m]);
			
				three_body+=(ctheta + (1./3.))*(ctheta + (1./3.))*exp((gamma/(ndistances[m]-a))+(gamma/(ndistances[n]-a)));
			}
		}
		while(!neighbours.empty()) {neighbours.pop_back();cell_neighbour_pos.pop_back();ndistances.pop_back();}
	}

	return epsilon*((0.5*A*two_body)+(lambda*three_body));
}

void insert_molecule(cell cell_list[], water_molecule gmol_list[], sim_box* sim_env_ptr, water_molecule temp, double energy) {
	/* Inserts molecule into the relevant molecule and cell lists. */
	
	gmol_list[sim_env_ptr->return_mols()] = temp; //Add the new particle
	cell_list[temp.cell_id].occupancy[cell_list[temp.cell_id].num_occupants] = sim_env_ptr->return_mols();
	cell_list[temp.cell_id].num_occupants+=1;
	sim_env_ptr->inc_mols();
	sim_env_ptr->inc_running_energy(energy);
	#ifdef MULTI
	if(sim_env_ptr->return_mols()>sim_env_ptr->return_max_mols()) sim_env_ptr->inc_max_mols();
	#endif

}

void decide_isertion (cell cell_list[],water_molecule gmol_list[], water_molecule temp, sim_box* sim_env_ptr){
	
	/* Calls the potential function to caculate the energy the new particle would contribute to the system,
	 * then works out the GC probability of acceptance. If this is greater than the random number we sent 
	 * to the function, perform the move and update the cell to reflect the new particle. The new particle 
	 * is added to the end of the gmol_list. Uses Metropolis method. */
		
	double deltaE, prob, prob_unbiased;

	#ifdef LJ
	  deltaE = LJPotential(temp, gmol_list, cell_list);
	#else
	  deltaE = SWPotential(gmol_list,cell_list,sim_env_ptr->return_num_cells_row(),temp);
	#endif
	 
	prob = (sim_env_ptr->return_volume()/((sim_env_ptr->return_mols())+1))*exp((mu-deltaE)); 
	
	#ifdef MULTI
	prob_unbiased = prob; 
	prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()+1);
	//prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()+1);
	//if(isnan(prob)) cout<<"Prob_unbiased "<<prob_unbiased<<" eta(N) "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols())<<" eta(N+1) "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols()+1)<<endl;

	#endif
	
	#ifdef MULTI_HAND
	prob_unbiased = prob;
	prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()+1);
	//cout<<"Prob before: "<<prob_unbiased<<" Prob after: "<<prob<<" Insertion with weights N: " <<sim_env_ptr->return_weight(sim_env_ptr->return_mols())<<" N+1: "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols()+1)<<endl;

	#endif
	

	if (isnan(prob)) {
		//cout<<"Volume: "<<sim_env_ptr->return_volume()<<" Num mols: "<<sim_env_ptr->return_mols()<<" deltaE: "<<deltaE<<endl;
	}

	if (prob>sim_env_ptr->get_rnd()) {
		//Perform the move
		#ifdef MULTI
		sim_env_ptr->update_Collection(sim_env_ptr->return_mols(),2,prob_unbiased, true);
		//sim_env_ptr->update_hits();
		#endif
		//cout<<"Insertion prob was "<<prob<<endl;
		insert_molecule(cell_list, gmol_list, sim_env_ptr, temp, deltaE);
		sim_env_ptr->inc_mca();
	}
	else {
		//Increment the rejected move count
		#ifdef MULTI
		sim_env_ptr->update_Collection(sim_env_ptr->return_mols(),1,prob_unbiased,false);
		//sim_env_ptr->update_hits();
		#endif
		sim_env_ptr->inc_mcr();
	}

}

void delete_molecule(cell cell_list[], water_molecule gmol_list[], sim_box* sim_env_ptr, water_molecule mol, double energy) {

	/* Deletes molecule and reorganises cell and molecule lists. */
	
	int cell_id, mol_id, to_shift, mol_space, error_cell, i;

	if(sim_env_ptr->return_mols()==1) {
		sim_env_ptr->dec_mols();
		cell_list[mol.cell_id].occupancy[0] = 0;
		cell_list[mol.cell_id].num_occupants-=1;
		sim_env_ptr->inc_mca();
		sim_env_ptr->dec_running_energy(energy); 
	}
	else {
		cell_id = gmol_list[sim_env_ptr->return_mols()-1].cell_id; //Find last molecule in list
			
		/*Update the cell of the last molecule to include its new id in the global
		* molecule list*/
		for(i=0;i<cell_list[cell_id].num_occupants;i++) {
			if(cell_list[cell_id].occupancy[i]==sim_env_ptr->return_mols()-1) {
			cell_list[cell_id].occupancy[i]=mol.id; 
			break;
			}
		}
		
		gmol_list[mol.id] = gmol_list[sim_env_ptr->return_mols()-1]; //Move the last particle to the empty position
		gmol_list[mol.id].id = mol.id;
		
		/*The particle at the end remains there, however as the number of molecules is
		* used to find particles to add and delete, the program will no longer
		* see the last molecule*/
		sim_env_ptr->dec_mols(); 
		for (i=0;i<cell_list[mol.cell_id].num_occupants;i++) {
		/*Finds the position of the molecule in the cell occupancy list and stores
		* this to move a molecule into it later*/ 
			if(cell_list[mol.cell_id].occupancy[i]==mol.id) {
				to_shift  = i; 
				break;
			}	
		}
		
		cell_list[mol.cell_id].occupancy[to_shift] = cell_list[mol.cell_id].occupancy[cell_list[mol.cell_id].num_occupants-1]; //Move molecule at end of list to newly empty space
		cell_list[mol.cell_id].occupancy[cell_list[mol.cell_id].num_occupants-1] = 0; //Remove molecule at end from occupancy
		cell_list[mol.cell_id].num_occupants-=1; //Decrement occupancy of cell
	
		sim_env_ptr->dec_running_energy(energy);
	}
}

void decide_deletion(cell cell_list[],water_molecule gmol_list[], water_molecule mol, sim_box* sim_env_ptr) {

	/* Calls the potential energy function to find the energy the particle provides to the system, 
	 * and calulates the probability of acceptance of removal. If this is greater than the random number we
	 * sent to the function, perform the move, and update the global lists to reflect the missing particle.
	 * The cell occupancy list must also be updated, the particle at the bottom of the list moved to the 
	 * newly empty position. The number of occupants of the cell must be decremented. */
	 
	double deltaE, prob, prob_unbiased;
	
	#ifdef LJ
	  deltaE = LJPotential(temp, gmol_list, cell_list);
	#else
	  deltaE = SWPotential(gmol_list,cell_list,sim_env_ptr->return_num_cells_row(),mol);
	#endif
	
	prob = (((sim_env_ptr->return_mols()))/ sim_env_ptr->return_volume())*exp(deltaE-mu); //Calculate probability

	#ifdef MULTI
	prob_unbiased = prob;
	//prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()-1);
	prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()-1);
if (isnan(prob)) {
		//cout<<"Prob_unbiased "<<prob_unbiased<<" eta(N) "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols())<<"eta(N-1) "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols()-1)<<endl;
	}
	#endif

	#ifdef MULTI_HAND
	prob_unbiased = prob;
	prob*=sim_env_ptr->return_weight(sim_env_ptr->return_mols())/sim_env_ptr->return_weight(sim_env_ptr->return_mols()-1);
	//cout<<"Prob before: "<<prob_unbiased<<" Prob after: "<<prob<<" Deletion with weights N: " <<sim_env_ptr->return_weight(sim_env_ptr->return_mols())<<" N-1: "<<sim_env_ptr->return_weight(sim_env_ptr->return_mols()-1)<<endl;
	
	#endif

	
	
	if (prob>sim_env_ptr->get_rnd())  { //Accept the move
		#ifdef MULTI
		sim_env_ptr->update_Collection(sim_env_ptr->return_mols(),0,prob_unbiased, true);
		//sim_env_ptr->update_hits();
		#endif
		//cout<<"Deletion prob was "<<prob<<endl;
		sim_env_ptr->inc_mca();
		delete_molecule(cell_list, gmol_list, sim_env_ptr, mol, deltaE);
		
	}
	
	else { //Reject the move
		#ifdef MULTI
		sim_env_ptr->update_Collection(sim_env_ptr->return_mols(),1,prob_unbiased, false);
		//sim_env_ptr->update_hits();
		#endif
		sim_env_ptr->inc_mcr();
	}
}

void decide_translation(water_molecule mol, water_molecule temp, water_molecule gmol_list[], cell cell_list[], sim_box* sim_env_ptr, double rnd, double dis[], int mol_cell) {
	
	/*Gives molecule a random displacement which is less than or equal to half the box size, then uses the GC
	 * acceptance rules to determine whether to move the particle or not. */
	 
	 double energy_before, energy_after, prob;
	 int num_cells_row = sim_env_ptr->return_num_cells_row();
	 
	#ifdef LJ
	  energy_before = LJPotential(mol, gmol_list, cell_list);
	#else
	  energy_before = SWPotential(gmol_list,cell_list, num_cells_row, mol); //Find energy before translation
	#endif
	 
	 mol = temp;
	 
	#ifdef LJ
	  energy_after = LJPotential(mol, gmol_list, cell_list);
	#else
	  energy_after = SWPotential(gmol_list,cell_list, num_cells_row, mol); //Find energy after translation
	#endif
	
	 prob = exp(energy_before-energy_after);
	 if(prob>sim_env_ptr->get_rnd()) {
		 gmol_list[mol.id] = mol;
		 sim_env_ptr->inc_running_energy(energy_after-energy_before);
		 sim_env_ptr->inc_mca();
		 cout<<"mca is "<<sim_env_ptr->return_mca()<<endl;
	}
	else sim_env_ptr->inc_mcr();
	 
}

void output_positions(water_molecule gmol_list[], int num_mols, string file_path) {

	/* Outputs the positions of each molecule*/
	ofstream positions; string out_name;

	out_name = file_path + "positions.dat";
	positions.open(out_name, ios::out);
	positions<<fixed<<setprecision(9);
	
	/* Go through the molecule list and output all data needed */
	for (int m=0;m<num_mols;m++) positions<<m<<" "<<gmol_list[m].x_pos<<" "<<gmol_list[m].y_pos<<" "<<gmol_list[m].z_pos<<" "<<gmol_list[m].cell_id<<endl;
	
	positions.close();
}

void output_positions_xyz(cell cell_list[], string file_path, water_molecule gmol_list[], int num_mols) {
	/*Positions in real space output as xyz file for VMD visualisation*/
	
	double x_pos, y_pos, z_pos; cell mcell;
	ofstream positions; string out_name;

	out_name = file_path + "positions.xyz";
	positions.open(out_name, ios::out);
	positions<<num_mols<<endl;
	positions<<"monatomic water"<<endl;
	
	/* Go through the molecule list and output all data needed */
	for (int m=0;m<num_mols;m++) {
		mcell = cell_list[gmol_list[m].cell_id];
		x_pos = (mcell.i+gmol_list[m].x_pos)*a*sigma;
		y_pos = (mcell.j+gmol_list[m].y_pos)*a*sigma;
		z_pos = (mcell.k+gmol_list[m].z_pos)*a*sigma;
		positions<<"O"<<" "<<setw(5)<<x_pos<<" "<<setw(5)<<y_pos<<" "<<setw(5)<<z_pos<<endl;
	}
	
	positions.close();
		
}

void output_data(int time_step, double save_dens, double save_mols, double save_energy, double acc_ratio, string file_path) {
	
	/*Outputs the densities of the system which were recorded throughout simulation*/
	ofstream data_out;
	data_out.open(file_path + "data.dat",ios::app);
	data_out<<time_step<<" "<<save_dens<<" "<<save_mols<<" "<<save_energy<<" "<<acc_ratio<<endl;
	data_out.close();
}

void output_ang_dist(water_molecule gmol_list[], sim_box sim_env, cell cell_list[], int iter) {

	/* Using same method as in Molinero and Moore 2009 paper, we aim to output the angular distribution of the 8 nearest 
	 * neighbours of each molecule. */
	 
	 int num_bins = 100, ncount = 0, max_r2_index = 0, ci, cj, ck, boundary = sim_env.return_num_cells_row()-1, num_cells_row = sim_env.return_num_cells_row();
	 int num_mols = sim_env.return_mols();
	 double bin_size = PI/num_bins, max_r2=1.0, r2;
	 double hist[num_bins][2], nntheta[num_mols][2][8], nnangles[num_mols][28], costheta;
	 cell ncell, mcell;
	 //Sets up the bins of the histogram
	 for (int i=0;i<num_bins;i++) {hist[i][0] = (i+1)*bin_size; hist[i][1]=0;}

	 /* Find the 8 nearest neighbouring particles and then calculate the angles between them and
	  * update the histogram*/
	 for(int m=0;m<num_mols;m++) { //Loop over all molecules
		 mcell = cell_list[gmol_list[m].cell_id]; //Cell of molecule we are looking at
		 for(int n=0;n<27;n++) { //Loop over all neighbouring cells
			 ncell = cell_list[mcell.neighbours[n][0]]; //Cell of neighbour
			 for(int mm=0;mm<ncell.num_occupants;mm++) { //Loop over all occupants of neighbouring cell
				if (m == ncell.occupancy[mm]) continue; //If the molecule in question is the same as the neighbour being investigated, move on
				r2 = Find_distance(gmol_list[m], gmol_list[ncell.occupancy[mm]], mcell.neighbours[n][1], mcell.neighbours[n][2], mcell.neighbours[n][3]);
				if(ncount<8) {nntheta[m][0][ncount]=ncell.occupancy[mm];nntheta[m][1][ncount] = r2; ncount++;if(r2>max_r2) {max_r2 = r2; max_r2_index=ncount-1;}}
				else {
					if(r2<max_r2) {
						nntheta[m][0][max_r2_index] = ncell.occupancy[mm];
						nntheta[m][1][max_r2_index] = r2;//need to initially set max_r2 and max_r2_index
						max_r2 = r2;
						for(int j=0;j<8;j++) { if(nntheta[m][1][j]>max_r2) {max_r2 = j; max_r2_index = j;} }
					}
				}
			}
		}
		
		if(ncount==8) {ncount=0; max_r2_index=0; max_r2=1; continue;} //If we have found the 8 nearest neighbours, we can move on
		else {
			/*Else we need to go the next nearest cells to find the neighbours. THIS DOESNT WORK HAS BUG*/
			for(int n=0;n<27;n++) {
				ncell=cell_list[mcell.neighbours[n][0]]; 
				for(int nn=0;nn<27;nn++) { //loop over neighbouring cells neighbouring cells
					if (mcell.neighbours[n][0]==ncell.neighbours[nn][0]) continue; //skip cells we have searched before
					else {
						for(int mm=0;mm<cell_list[ncell.neighbours[nn][0]].num_occupants;mm++) {
							ci = ncell.neighbours[nn][1]-mcell.neighbours[n][1]; 
							cj = ncell.neighbours[nn][2]-mcell.neighbours[n][2]; 
							ck = ncell.neighbours[nn][3]-mcell.neighbours[n][3]; 
							
							if(ci == boundary) ci = -1;
							if(ci == -1*boundary) ci = 1;
							if(cj == boundary) cj = -1;
							if(cj == -1*boundary) cj = 1;
							if(ck == boundary) ck = -1;
							if(ck == -1*boundary) ck = 1;
							
							r2 = Find_distance(gmol_list[m], gmol_list[cell_list[ncell.neighbours[nn][0]].occupancy[mm]], ci, cj, ck, num_cells_row, true);
							if(ncount<8) {nntheta[m][0][ncount]=ncell.occupancy[mm]; nntheta[m][1][ncount] = r2; ncount++; if(r2>max_r2) {max_r2 = r2; max_r2_index=ncount-1;}}
							else {
								if(r2<max_r2) {
									nntheta[m][0][max_r2_index] = ncell.occupancy[mm];
									nntheta[m][1][max_r2_index] = r2;//need to initially set max_r2 and max_r2_index
									max_r2 = r2;
									for(int j=0;j<8;j++) { if(nntheta[m][1][j]>max_r2) { max_r2 = j; max_r2_index = j;} }
								}
							}
						}
					}
				}
			}
		}
		ncount = 0; max_r2=1.0; max_r2_index=0;
	}
	/* Once we have 8 nearest neighbours, we can find the 28 angles they make*/
	for(int m=0;m<num_mols;m++) {
		ncount = 0;
		for (int i=0;i<8;i++) {
			for (int j=i+1;j<8;j++) {
				ci = cell_list[gmol_list[(int)nntheta[m][0][j]].cell_id].i-cell_list[gmol_list[(int)nntheta[m][0][i]].cell_id].i; 
				cj = cell_list[gmol_list[(int)nntheta[m][0][j]].cell_id].j-cell_list[gmol_list[(int)nntheta[m][0][i]].cell_id].j;
				ck = cell_list[gmol_list[(int)nntheta[m][0][j]].cell_id].k-cell_list[gmol_list[(int)nntheta[m][0][i]].cell_id].k;
				
				if(ci == boundary) ci = -1;
				if(ci == -1*boundary) ci = 1;
				if(cj == boundary) cj = -1;
				if(cj == -1*boundary) cj = 1;
				if(ck == boundary) ck = -1;
				if(ck == -1*boundary) ck = 1;
				
				//Overall I think the problem here is that we still haven't found 8 neighbours hence costheta is going above one
				r2 = Find_distance(gmol_list[(int)nntheta[m][0][i]],gmol_list[(int)nntheta[m][0][j]],ci,cj,ck,num_cells_row, true);
				costheta = Find_angle(r2*aa,sqrt(nntheta[m][1][i]*aa),sqrt(nntheta[m][1][j]*aa));
				nnangles[m][ncount] = acos(costheta); //Outputs angle in radians
				ncount++;
			}
		}
	}
	/* And now we bin the angles accordingly */
	for(int m=0;m<num_mols;m++) {
		for(int i=0;i<28;i++) {
			for(int j=0;j<num_bins;j++) {
				if(j==0 and nnangles[m][i]<hist[j][0]) {hist[j][1]+=1; break;}
				else if(nnangles[m][i]>hist[j][0] and nnangles[m][i]<=hist[j+1][0]) { hist[j][1]+=1; break; } 
			}
		}
	}
	
	/* Now normalise the outputs for P(theta)*/
	for(int i=0;i<num_bins;i++) hist[i][1]/=(num_mols*28);
	
	/* And finally output to a file*/
	ofstream ang_dist;
	ang_dist.open(sim_env.return_file_path() + "angular_distribution_" + to_string(iter) + ".txt", ios::out);

	for (int m=0;m<num_bins;m++) ang_dist<<hist[m][0]*(180.0/PI)<<" "<<hist[m][1]<<endl;
	
	ang_dist.close();

	
}

double output_gr(water_molecule gmol_list[], sim_box sim_env, cell cell_list[], int iter) {
	
	/*Goes through the global molecule list and finds the g(r) curves for each particle. Normalises these and puts them
	 * into bins. Outputs this into a file for gnuplot or python visulisation.  */
	 
	cell mcell,ncell;
	int list_size, num_nn=0, ci, cj, ck, num_cells_row = sim_env.return_num_cells_row(), boundary = num_cells_row-1;
	int num_mols = sim_env.return_mols();
	double r2, min_r=0, max_r=(num_cells_row/2)*a*sigma, bin_size, nn;
	deque<double> r;
	double bin_vol;

	/*Goes through molecule list and works out the distance between each molecule pair*/
	for (int m=0;m<num_mols;m++) {
		mcell = cell_list[gmol_list[m].cell_id];
		for (int n=0;n<num_mols;n++) {
			if (m == n) continue;
			ncell = cell_list[gmol_list[n].cell_id];
			
			ci = cell_list[gmol_list[n].cell_id].i-cell_list[gmol_list[m].cell_id].i; 
			cj = cell_list[gmol_list[n].cell_id].j-cell_list[gmol_list[m].cell_id].j;
			ck = cell_list[gmol_list[n].cell_id].k-cell_list[gmol_list[m].cell_id].k;
				
			if(ci == boundary) ci = -1;
			if(ci == -1*boundary) ci = 1;
			if(cj == boundary) cj = -1;
			if(cj == -1*boundary) cj = 1;
			if(ck == boundary) ck = -1;
			if(ck == -1*boundary) ck = 1;
			
			r2 = Find_distance(gmol_list[m], gmol_list[n], ci, cj, ck, num_cells_row, true);
			r2 = r2*a*a*sigma*sigma;
			if(r2<=12.25 and r2!=0) num_nn++;
			if (r2!=0) r.push_back(sqrt(r2));
		}
	}
	
	nn = num_nn/(double)num_mols; //Still something slightly wrong here
	cout<<"Num nearest neighbours "<<nn <<endl;
	
	int num_bins = 100;
	bin_size =  (max_r - min_r)/num_bins;
	list_size = r.size();
	double hist_bin[num_bins],hist_r[num_bins];;
	
	/*Sets bins*/
	for(int m=0;m<num_bins;m++) {
		hist_bin[m] = (m+1)*bin_size;
		hist_r[m]=0;
	}
	
	/*Sorts distances into bins*/
	for (int m=0;m<r.size();m++) {
		for(int n=0;n<num_bins;n++){
			if (r[m]>=hist_bin[n] and r[m]<hist_bin[n+1]) hist_r[n]++;
			else continue;
		}
	}
	
	double mid_bin;
	ofstream gr;
	gr.open(sim_env.return_file_path() + "gr" + to_string(iter) + ".txt", ios::app);
	double rho = (num_mols-1)/sim_env.return_volume();
	
	
	/*Normalises the bins and writes to a file*/
	for (int m=0;m<num_bins;m++) {
		mid_bin = hist_bin[m] - (bin_size/2);
		bin_vol = 4*PI*(mid_bin*mid_bin)*bin_size;
		hist_r[m]=(hist_r[m]/(num_mols));
		hist_r[m] = hist_r[m]/(rho*bin_vol);
		gr<<hist_bin[m]<<" "<<hist_r[m]<<endl;
	}
	
	gr.close();
	return nn;
	
}

void output_sim_data(sim_box sim_env, double end_energy, double num_nn){
		
		/*At the end of the simulation, output the details which would be needed to restart the simulation*/
		ofstream sim_data;
		sim_data.open(sim_env.return_file_path() +"OUTPUT",ios::out);
		#ifdef LJ
		sim_data<<"SIM_TYPE LJ"<<endl;
		#else
		sim_data<<"SIM_TYPE mW"<<endl;
		#endif
		sim_data<<"NUM_CELLS_ROW "<<sim_env.return_num_cells_row()<<endl;
		sim_data<<"NUM_MOLS "<<sim_env.return_mols()<<endl;
		sim_data<<"SEED "<<sim_env.return_seed()<<endl;
		sim_data<<"CURRENT_STEPS "<<sim_env.return_steps()+sim_env.return_current_steps()+sim_env.return_steps_eqb()<<endl;
		sim_data<<"MCA "<<sim_env.return_mca()<<endl;
		sim_data<<"MCR "<<sim_env.return_mcr()<<endl;
		sim_data<<"END_ENERGY "<<end_energy<<endl;
		sim_data<<"MU "<<mu<<endl;
		sim_data<<"NEAREST_NEIGHBOURS "<<num_nn<<endl;
		sim_data.close();
		
}

void Perform_MC_Eqb(cell cell_list[], water_molecule gmol_list[], sim_box* sim_env_ptr) {
	
	/*To speed up the program, we have a seperate equilibration function, which does not check if we need to output data.*/
	
	int rnd_check=0, mol_num=0, pcell,m, steps, num_cells_row = sim_env_ptr->return_num_cells_row(), num_cells = sim_env_ptr->return_num_cells();
	int save_frq = sim_env_ptr->return_save_frq();
	double x_pos, y_pos, z_pos, *dummy_ptr, dummy_energy, dis[3];
	water_molecule mol, temp, *mol_ptr;

	dummy_ptr = &dummy_energy;
	dummy_energy=0;

	for (steps=0;steps<sim_env_ptr->return_steps_eqb();steps++) {
		if(steps%save_frq==0) {
			cout<<steps<<" equilibration attempts complete"<<endl;
			#ifdef MULTI
			sim_env_ptr->update_Transition();
			//cout<<"Updating weights"<<endl;
			//sim_env_ptr->update_weights();
			#endif
		}

		#ifdef TRANSLATION
		if (sim_env_ptr->get_rnd() > 0.5) { //Attempt insertion/deletion
		#endif
			/*If random_numbers>0.5, attempt insertion. */
			if (sim_env_ptr->get_rnd() > (0.5)) {
				pcell = floor(sim_env_ptr->get_rnd()*(num_cells));
				x_pos = sim_env_ptr->get_rnd(); y_pos = sim_env_ptr->get_rnd(); z_pos = sim_env_ptr->get_rnd();
				temp = water_molecule(sim_env_ptr->return_mols(), x_pos,y_pos,z_pos,pcell, false);
				
				decide_isertion(cell_list,gmol_list,temp,sim_env_ptr);
				
			}
		
			/*else attempt deletion*/
			else {
				if(sim_env_ptr->return_mols() == 0) {
					sim_env_ptr->inc_mcr();
					#ifdef MULTI
					sim_env_ptr->update_Collection(0, 0, 1, true);
					#endif
					continue;
				}
			
				mol_num = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_mols()));
				temp = gmol_list[mol_num];
				decide_deletion(cell_list, gmol_list, temp, sim_env_ptr);
			}	
		#ifdef TRANSLATION
		}
		else { //Attempt translation
	
			/*Start by generating random displacements*/
			dis[0] = (sim_env_ptr->get_rnd()-0.5);
			dis[1] = (sim_env_ptr->get_rnd()-0.5);
			dis[2] = (sim_env_ptr->get_rnd()-0.5);
			
			/*Pick a random molecule to displace*/
			mol_num = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_mols()));
			mol = gmol_list[mol_num];
			
			/*Now work out what new position this corresponds to in the cell system*/
			temp = Find_New_Pos(dis,cell_list,mol,*sim_env_ptr);
		
			decide_translation(mol, temp, gmol_list, cell_list, sim_env, sim_env_ptr->get_rnd(), dis, pcell);
		}
		#endif
		
	}		
	cout<<"Equilibration Complete"<<endl;
}


void Perform_MC(cell cell_list[], water_molecule gmol_list[], sim_box* sim_env_ptr) {
	
	/*Performs the MC by attempting insertion and deletion. Also outputs the density, energy and g(r) at given intervals.*/
	
	int mol_num=0, pcell, save_frq = sim_env_ptr->return_save_frq(), output_gr_frq=floor(sim_env_ptr->return_steps()/sim_env_ptr->return_num_gr());
	int num_entries = floor(sim_env_ptr->return_steps()/save_frq),m, steps;
	int id, num_cells_row = sim_env_ptr->return_num_cells_row(), num_cells = sim_env_ptr->return_num_cells();
	double x_pos, y_pos, z_pos, avg_dens=0, end_energy, avg_energy=0, avg_nn=0, dens=0, dis[3], acc_ratio, volume = sim_env_ptr->return_volume();
	water_molecule mol, temp;
	
	for (steps=0;steps<sim_env_ptr->return_steps();steps++) {
		
		#ifdef TRANSLATION
		if(sim_env_ptr->get_rnd()>0.5) {
		#endif
			if (sim_env_ptr->get_rnd() > 0.5) {
				pcell = floor(sim_env_ptr->get_rnd()*(num_cells));
				x_pos = sim_env_ptr->get_rnd(); y_pos = sim_env_ptr->get_rnd(); z_pos = sim_env_ptr->get_rnd();
				temp = water_molecule(sim_env_ptr->return_mols(), x_pos,y_pos,z_pos,pcell, false);
				
				decide_isertion(cell_list, gmol_list, temp, sim_env_ptr);
			
			}
			else {
				
				if(sim_env_ptr->return_mols() == 0) {
					sim_env_ptr->inc_mcr();
					#ifdef MULTI
					sim_env_ptr->update_Collection(0, 0, 1, true);
					#endif
				}
				else {
				mol_num = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_mols()));
				temp = gmol_list[mol_num];
				id = temp.id;
				
				decide_deletion(cell_list, gmol_list, temp, sim_env_ptr);
				}
				
			
			}
		#ifdef TRANSLATION
		}
		
		else { //Attempt translation
		
			/*Start by generating random displacements*/
			dis[0] = (sim_env_ptr->get_rnd()-0.5);
			dis[1] = (sim_env_ptr->get_rnd()]-0.5);
			dis[2] = (sim_env_ptr->get_rnd()-0.5);
		
			/*Pick a random molecule to displace*/
			mol_num = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_mols()));
			mol = gmol_list[mol_num];
		
			/*Now work out what new position this corresponds to in the cell system*/
			temp = Find_New_Pos(dis,cell_list,mol,*sim_env_ptr);
		
			decide_translation(mol, temp, gmol_list, cell_list, sim_env_ptr, dis, pcell);
		}
		#endif	
		if (steps%save_frq == 0) {
			cout<<steps<<" attempts complete"<<" num_mols is "<<sim_env_ptr->return_mols()<<endl;
			#ifdef LJ
			dens = (double)sim_env_ptr->return_mols()/volume;
			#else
			dens=((double)sim_env_ptr->return_mols()/(volume*sigma*sigma*sigma))*(MW/NA)*10;
			#endif
			#ifdef MULTI
			sim_env_ptr->update_Transition();
			sim_env_ptr->update_weights();
			#endif
			avg_energy += sim_env_ptr->return_running_energy();
			avg_dens+=dens;
			acc_ratio = ((double) sim_env_ptr->return_mca()/((double)(sim_env_ptr->return_mca())+(double)(sim_env_ptr->return_mcr())))*100.0;
			output_data(steps, dens, sim_env_ptr->return_mols(), sim_env_ptr->return_running_energy(), acc_ratio, sim_env_ptr->return_file_path());
		}	
	
		if(steps%output_gr_frq==0) {
			avg_nn += output_gr(gmol_list, *sim_env_ptr, cell_list, steps);
			if (sim_env_ptr->return_mols()>75) output_ang_dist(gmol_list,*sim_env_ptr, cell_list,steps);
		}

	}

	/*At the end of the simulation we output some statistics to the screen for the user.*/
	cout<<"num entries is "<<num_entries<<endl;
	#ifdef LJ
	end_energy = Find_System_Energy_LJ(gmol_list,cell_list, sim_env_ptr->return_mols());
	#else
	end_energy = Find_System_Energy_SW(gmol_list, cell_list, num_cells_row, sim_env_ptr->return_mols());
	#endif
	cout<<"mca: "<<sim_env_ptr->return_mca()<<" mcr: "<<sim_env_ptr->return_mcr()<<endl;
	cout<<"steps is "<<steps<<" current num steps "<<sim_env_ptr->return_current_steps()<<endl;
	cout<<"Number of molecules was "<<sim_env_ptr->return_mols()<<endl;
	#ifdef MULTI
	cout<<"Max mols was "<<sim_env_ptr->return_max_mols()<<endl;
	#endif
	cout<<"Average density was "<<avg_dens/num_entries<<endl; 
	cout<<"Average energy was "<<avg_energy/num_entries<<endl;
	cout<<"Running energy was "<<sim_env_ptr->return_running_energy()<<endl;
	cout<<"End energy was "<<end_energy<<endl;
	cout<<"Number of nearest neighbours "<<avg_nn/(double)sim_env_ptr->return_num_gr()<<endl;
	
	//Here we output the final sim data so it can be reloaded in the future. Should we also put the energy calculations to compare later?
	output_sim_data(*sim_env_ptr, end_energy, avg_nn/(double)sim_env_ptr->return_num_gr()); 
	#ifdef MULTI
	sim_env_ptr->write_weights();
	#endif
	
}

/* This is quite outdated, because the load file actually contains the number of moelcules. UPDATE*/
void load_mols(string in_name, water_molecule gmol_list[], cell cell_list[]) {
	
	/*Load molecules from the positions file of a previous simulation.*/
	
	cout<<"Loading molecules"<<endl;
	
	ifstream positions(in_name); 
	int id, cell_id, m=0;
	float x_pos,y_pos,z_pos;

	while(positions>>id>>x_pos>>y_pos>>z_pos>>cell_id) {
		gmol_list[m] = water_molecule(id,x_pos,y_pos,z_pos,cell_id, true);
		cell_list[cell_id].occupancy[cell_list[cell_id].num_occupants] = m;
		cell_list[cell_id].num_occupants+=1;
		m++;
	}
	
	positions.close();
	
	cout<<"Loaded molecules successfully"<<endl;

}

void Setup_Sim_Box(cell cell_list[], sim_box sim_env) {
	
	/* We assume a square box. The simulation system must be set up so that we have the correct number
	 * of cells and every cell knows who's its neighbours are. */	
	 
	cout<<"This program assumes a square box."<<endl;
	cout<<"About to generate cells"<<endl;
	
	int num_neighbours=0, ci,cj,ck, correct_neighb=0, id=0;
	int num_cells_row = sim_env.return_num_cells_row(), num_cells_row_bound = num_cells_row-1;
	int num_cells_row_sqr = num_cells_row*num_cells_row; 
	deque<int> cell_distances;
	
	/* First we just generate the cell. Due to the way cell lists work, the length of the cell
	 * should be the cut off length of interaction, which is a*sigma here. Every cell has an
	 * id which is its position in the cell_list array. The number of occupants is initially
	 * set to 0, and its position in real space is stored (given by i,j,k). */
	 
	for (int i=0;i<num_cells_row;i++) {
		for (int j=0;j<num_cells_row;j++) {
			for (int k=0;k<num_cells_row;k++) {
				cell_list[id] = cell(id,i,j,k,0); 
				id++;
			}
		}
	}
	
	/* Now we need to find the neighbours for the cells. Memory is cheap so these are stored
	 * rather than found each time (significantly increasing the speed of the program). At the
	 * moment this was designed for a 4x4 box and whilst every effort to make it more general
	 * is being taken, this has yet to be tested for other systems. */
	
	for(int c =0; c<sim_env.return_num_cells(); c++) {
		num_neighbours=0;
		for(int cc=0; cc<sim_env.return_num_cells();cc++) {
			if (cc==c) {
				cell_list[c].neighbours[num_neighbours][0]=cc;
				cell_list[c].neighbours[num_neighbours][1]=0;
				cell_list[c].neighbours[num_neighbours][2]=0;
				cell_list[c].neighbours[num_neighbours][3]=0;
				num_neighbours++;
			}
			
			else {
				
				/* Here we assume periodic boundary conditions. */
				ci = ((int)floor(cell_list[cc].id/num_cells_row_sqr)) - ((int)floor(cell_list[c].id/num_cells_row_sqr));
				cj = ((int)floor(cell_list[cc].id/num_cells_row)%num_cells_row) - ((int)floor(cell_list[c].id/num_cells_row)%num_cells_row);
				ck = cell_list[cc].id%num_cells_row - cell_list[c].id%num_cells_row;
				
				/* If we are at the boundary, we need to manually change the neighbour positions to be correct. */
				if(ci==num_cells_row_bound) ci =-1; 
				if(ci==-1*num_cells_row_bound) ci=1;
				if(cj==num_cells_row_bound) cj=-1;
				if(cj==-1*num_cells_row_bound) cj=1;
				if(ck==num_cells_row_bound) ck=-1;
				if(ck==-1*num_cells_row_bound) ck=1;
				
				/* Finally, if the neighbour is found, add it to the neighbour list. */
				if (ci>=-1 and ci<=1 and cj>=-1 and cj<=1 and ck>=-1 and ck<=1) {
					cell_list[c].neighbours[num_neighbours][0]=cc;
					cell_list[c].neighbours[num_neighbours][1]=ci;
					cell_list[c].neighbours[num_neighbours][2]=cj;
					cell_list[c].neighbours[num_neighbours][3]=ck;
					num_neighbours++;
				
				}
			}
		}
	}
	
	cout<<"Cells generated successfully."<<endl;
}

/*Loading previous simulations will currently not work as the form of some classes and where things are kept has changed.
 * To be updated. */
void Place_Initial_Mols(water_molecule gmol_list[], sim_box* sim_env_ptr, bool load, cell cell_list[]){
	
	/* To start the simulation we either need to place molecules or load them from a previous configuration.*/
	
	int m=0, pcell, num_cells = sim_env_ptr->return_num_cells();
	double x_pos,y_pos,z_pos;
	string in_name;
	
	if(load) in_name = sim_env_ptr->return_in_file_path() + "positions.dat";
	else  in_name = sim_env_ptr->return_file_path() + "positions.dat";
	
	/* If we just need to load a previous configuration, the case is simple: */
	if (load) load_mols(in_name, gmol_list, cell_list);
	else {
		
		/*If we are placing molecules randomly, we just generate positions and fill the gmol_list. Currently
		 * we have nothing to prevent two particles overlapping, however the potential in the MC simulation
		 * normally sorts out any problems during equilibration. */

		cout<<"Placing initial molecules"<<endl;
		while (m<sim_env_ptr->return_mols()) {
		
			pcell = floor(sim_env_ptr->get_rnd()*num_cells);
			x_pos = sim_env_ptr->get_rnd(); y_pos = sim_env_ptr->get_rnd(); z_pos = sim_env_ptr->get_rnd();
			cell_list[pcell].occupancy[cell_list[pcell].num_occupants]=m; 
			cell_list[pcell].num_occupants+=1;
			gmol_list[m]=water_molecule(m, x_pos,y_pos,z_pos, pcell, true);
			m++; 
		}
	
	cout<<"Successfully placed molecules."<<endl;
	
	}
	
}

int main(int argc, char* argv[]) {
	
	time_t *seed_time = new time_t;
	
	int num_cells_row, seed, num_steps, num_steps_eqb=0, num_mols=0, mca=0, mcr=0, current_num_steps=0, save_frq, num_gr;
	string file_path, arg1, arg2, test_string, load_string, load_file_path;
	bool test, load;
	
	/*To start we load the simulation information from the INPUT file. This should be made into a function.*/
	cout<<"Loading input file"<<endl;
	file_path = argv[1];
	cout<<"FILE_PATH IS "<<file_path<<endl;
	cout<<"OPENING "<<file_path<<"INPUT"<<endl;

	ifstream input(file_path + "INPUT");
	
	while (input>>arg1>>arg2) {
		cout<<arg1<<" "<<arg2<<endl;
		if (!arg1.compare("LOAD")) {
			if(!arg2.compare("FALSE")) load=false;
			else { load=true; }
		}
		else if (!arg1.compare("FILE_PATH")) file_path = arg2;
		else if (!arg1.compare("IN_PATH")) load_file_path = arg2;
		else if (!arg1.compare("SAVE_FREQUENCY")) save_frq = atoi(arg2.c_str());
		else if (!arg1.compare("NUMBER_OF_GR")) num_gr = atoi(arg2.c_str());
		else if (!arg1.compare("NUM_STEPS_EQUILIBRATION")) num_steps_eqb = atoi(arg2.c_str());
		else if (!arg1.compare("NUM_STEPS")) {
			num_steps = atoi(arg2.c_str()); 
			if (load) break;
		}
		else if(!arg1.compare("NUM_CELLS_ROW")) {num_cells_row=atoi(arg2.c_str()); cout<<"num_cells_row is "<<num_cells_row<<endl;}
		else if (!arg1.compare("TEST")) {
			if(!arg2.compare("FALSE")) test=false;
			else test=true;
			cout<<"test is "<<test<<endl;
		}
		else if (!arg1.compare("NUM_MOLS")) num_mols = atoi(arg2.c_str());
		else {
			cout<<"Could not find argument, please check input file. Arguments are: "<<endl;
			cout<<"LOAD LOAD FILE_PATH IN_PATH SAVE_FRQUENCY NUMBER_OF_GR NUM_STEPS_EQUILIBRATION NUM_STEPS NUM_MOLS TEST NUM_CELLS_ROW"<<endl;
			exit(EXIT_FAILURE);

		}
	}
	input.close();
	if (load) {
		/* Currently incompatible with new layout. To be updated. */
		ifstream load_data(load_file_path+"OUTPUT");
		while (load_data>>arg1>>arg2) {
			cout<<arg1<<" "<<arg2<<endl;
			if(!arg1.compare("SIM_TYPE")) {
				#ifdef LJ
					if(arg2.compare("LJ")) {cout<<"Please define LJ in the program."<<endl; exit(0);}
				#endif
				continue;
			}
			else if(!arg1.compare("NUM_CELLS_ROW")) {num_cells_row=atoi(arg2.c_str()); cout<<"num_cells_row is "<<num_cells_row<<endl;}
			else if (!arg1.compare("NUM_MOLS")) num_mols = atoi(arg2.c_str());
			else if (!arg1.compare("SEED")) seed = atoi(arg2.c_str());
			else if (!arg1.compare("CURRENT_STEPS")) current_num_steps = atoi(arg2.c_str());
			else if (!arg1.compare("MCA")) mca = atoi(arg2.c_str());
			else if (!arg1.compare("MCR")) mcr = atoi(arg2.c_str());
			else if (!arg1.compare("END_ENERGY")) continue;
			else if (!arg1.compare("MU")) continue;
				/*if(atoi(arg2.c_str())!=mu){
					cout<<"Please enter the correct mu as found in the output file."<<endl;
					exit(0);
				}
			}*/
			else if (!arg1.compare("NEAREST_NEIGHBOURS")) continue;
			else {
				cout<<"Could not find argument, please check input file. Arguments are: "<<endl;
				cout<<"SIM_TYPE NUM_CELLS_ROW NUM_MOLS SEED CURRENT_NUM_STEPS MCA MCR NEAREST_NEIGHBOURS"<<endl;
				exit(EXIT_FAILURE);
			}
		}
		cout<<"seed is "<<seed<<endl;
		load_data.close();
	}
	
	cout<<"current num steps is "<<current_num_steps<<endl;
	cout<<" mca is "<<mca<<" and mcr "<<mcr<<endl;
	
	cout<<"file path is "<<file_path<<endl;
	
	/* The random number generator used is supplied by ranvec.h, and here we initialise it. */

	if (test) { seed = 1522245984;}
	else if (!load) {
		//seed = atof(argv[2]);
		time(seed_time);
		seed = (int)*seed_time;
		ofstream seeds;
		string out_seed_name;
		out_seed_name = file_path+ "seed.dat";
		seeds.open(out_seed_name,ios::app);
		seeds<<seed<<endl;
		seeds.close();
	}
	cout<<"num_mols is "<<num_mols<<endl;
	double dens, w, weights_in[max_N]; int m=0;
	ifstream weight_func(load_file_path + "weights"); 
		
	while(weight_func>>dens>>w) {
		if(abs(w)>0.000000001) weights_in[m] = w;
		else weights_in[m] = 0.000000001;
		m++;
		
	}
	
	weight_func.close();
	for(int i=m;i<max_N;i++) weights_in[i] = weights_in[m-1];
	sim_box sim_env(num_cells_row, num_mols, mca, mcr, current_num_steps, num_steps_eqb, num_steps, num_gr, save_frq, seed, file_path, load_file_path, weights_in);
	cout<<"Simulation environment created"<<endl;
	
	
	cell cell_list[sim_env.return_num_cells()];
	Setup_Sim_Box(cell_list, sim_env);

	int array_size=num_cells_row*1000;
	water_molecule gmol_list[array_size];
	sim_box* sim_env_ptr; sim_env_ptr = &sim_env;
	Place_Initial_Mols(gmol_list,sim_env_ptr,load,cell_list);
	cout<<"Starting Simulation"<<endl;
	
	Perform_MC_Eqb(cell_list, gmol_list, sim_env_ptr);
	double energy;
	#ifdef LJ
	energy = Find_System_Energy_LJ(gmol_list,cell_list,sim_env_ptr->return_mols());
	#else
	energy = Find_System_Energy_SW(gmol_list,cell_list,num_cells_row,sim_env_ptr->return_mols());
	#endif
	sim_env_ptr->set_running_energy(energy);
	
	cout<<"System energy is "<<sim_env.return_running_energy()<<endl;
	Perform_MC(cell_list, gmol_list, sim_env_ptr);
	output_positions(gmol_list, sim_env.return_mols(),file_path);
	output_positions_xyz(cell_list, file_path, gmol_list, sim_env.return_mols());
}

