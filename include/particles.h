#pragma once

#include <string.h>

#ifdef __SSE__
#include <pmmintrin.h>
#endif

#include "tinyvec3d.h"

#include "config.h"

#include "lib.h"
#include "constant.h"

#include "geometry.h"
#include "pdp3Time.h"
#include "chargeDensity.h"
#include "current.h"
#include "eField.h"
#include "hField.h"

using namespace std;

class ParticlesList;

class EField;

class HField;

class Particles
{
public:
  Particles(void);

  Particles(char const *p_name, double p_charge, double p_mass, int p_number,
            Geometry *geom);
  Particles(Particles& cp_particles);
  ~Particles();
public:
  // The specie name
  char *name;
  // The specie charge in electron charges
  double charge;
  // The specie  *mass
  double mass; // in electron masses
  double *mass_array;
  double *charge_array;
  double init_const_mass;
  // Number of particles
  unsigned int number;

  //! Array of particle position
  //! in format:
  //! \f$ [ [r_1, \phi_1, z_1 ], [ [r_2, \phi_2, z_2 ] ... ] \f$ and
  //! So, it is 2D array from particles
  //! which includes arrays of coordinate
  //! components per particle
  //!
  //! You can call: pos[particle_number][component_number]
  //! to get required coordinate component:
  //! \f$ 0 for r, 1 for \phi (zeros), 2 for z \f$ )
  double **pos;

  //! Array of particle old position
  //! (on previous step. Required for current weighting)
  //! in format:
  //! \f$ [ [r_1, \phi_1, z_1 ], [ [r_2, \phi_2, z_2 ] ... ] \f$ and
  //! So, it is 2D array from particles
  //! which includes arrays of coordinate
  //! components per particle
  //!
  //! You can call: pos_old[particle_number][component_number]
  //! to get required coordinate component:
  //! \f$ 0 for r, 1 for \phi (zeros), 2 for z \f$ )
  double **pos_old;

  //! Array of particle velocities
  //! in format:
  //! \f$ [ [v_{r_1}, v_{\phi_1}, v_{z_1} ], [ [v_{r_2}, v_{\phi_2}, v_{z_2} ] ... ] \f$
  //! So, it is 2D array from particles
  //! which includes arrays of velocity
  //! components per particle
  //!
  //! You can call: vel[particle_number][component_number]
  //! to get required velocity component:
  //! \f$ 0 for r, 1 for \phi, 2 for z \f$
  double **vel;

  //indicator if particle is still alive
  bool *is_alive;

  //current density
  //temporaty member of Particle class
  //created for the purpose of integration
  //with exsisting routine of Maxwell equations integration
  Geometry *geom1;
  ParticlesList *p_list;

private:

  //! service array to convert position from xy to rz pane
  double *sin_theta_r;

  //! service array to convert position from xy to rz pane
  double *cos_theta_r;

public:
/*   void charge_weighting(ChargeDensity *ro1); */
/*   void step_v(EField *e_fld, HField *h_fld, Time *t); */
/*   void half_step_pos(Time *t); */
/*   virtual void reflection(void); */
/*   void set_j_0(); */
/*   void set_v_0(); */
/*   void set_x_0(); */
/*   void velocity_distribution(double therm_vel); */
/*   void load_spatial_distribution(double n1, */
/*                                  double n2, */
/*                                  double left_plasma_boundary, */
/*                                  int type); */
/*   void load_spatial_distribution_with_variable_mass(double n1, */
/*                                                     double n2, */
/*                                                     double left_plasma_boundary, */
/*                                                     int type); */
/*   void load_velocity_distribution(double v_thermal); */
/*   void simple_j_weighting(Time *time1, */
/*                           Current *j1, */
/*                           double x1_new, */
/*                           double x3_new, */
/*                           double x1_old, */
/*                           double x3_old, */
/*                           int i_n, */
/*                           int k_n, */
/*                           int p_number); */
/*   void simple_constrho_j_weighting(Time *time1, */
/*                                    Current *j1, */
/*                                    double x1_new, */
/*                                    double x3_new, */
/*                                    double x1_old, */
/*                                    double x3_old, */
/*                                    int i_n, */
/*                                    int k_n, */
/*                                    int p_number); */
/*   void j_weighting(Time *time1, Current *j1); */
/*   void strict_motion_weighting(Time *time1, */
/*                                Current *j1, */
/*                                double x1_new, */
/*                                double x3_new, */
/*                                double x1_old, */
/*                                double x3_old, */
/*                                int p_number); */
/*   void azimuthal_j_weighting(Current *j1); */
/*   void set_simple_cell(int **cell_arr_jr, */
/*                        int **cell_arr_jz, */
/*                        int start_number, */
/*                        int i_new, */
/*                        int k_new); */
/*   void get_cell_numbers_jr_2(double x1_new, */
/*                              double x3_new, */
/*                              double x1_old, */
/*                              double x3_old, */
/*                              int **cell_arr_jr, */
/*                              int **cell_arr_jz, */
/*                              int *number); */
/*   void get_cell_numbers_jr_1(double x1_new, */
/*                              double x3_new, */
/*                              double x1_old, */
/*                              double x3_old, */
/*                              int *i_return, */
/*                              int *k_return, */
/*                              int *accur); */

/*   void back_position_to_rz(); */
/*   void back_velocity_to_rz(); */

/*   void dump_position_to_old(); */

/*   void step_v_single(EField *e_fld, HField *h_fld, Time *t, unsigned int i); */
/*   void dump_position_to_old_single(unsigned int i); */
/*   void half_step_pos_single(Time *t, unsigned int i); */
/*   void back_position_to_rz_single(unsigned int i); */
/*   void boris_pusher(EField *e_fld, HField *h_fld, Time *t); */

/*   void reflection_single(unsigned int i); */
/*   void move_half_reflect(Time *t); */

  void set_v_0();
  void set_x_0();
  void reflection();
  void half_step_pos(Time *t);
  void charge_weighting(ChargeDensity *ro1);
  void velocity_distribution(double tempr_ev);
  void load_spatial_distribution(double n1, double n2, double left_plasma_boundary, int type);
  void simple_j_weighting(Time *time1,
                          Current *j1,
                          double x1_new,
                          double x3_new,
                          double x1_old,
                          double x3_old,
                          int i_n,
                          int k_n,
                          int p_number);
  void simple_constrho_j_weighting(Time *time1,
                                   Current *j1,
                                   double x1_new,
                                   double x3_new,
                                   double x1_old,
                                   double x3_old,
                                   int i_n,
                                   int k_n,
                                   int p_number);
  void j_weighting(Time *time1, Current *j1);
  void azimuthal_j_weighting(Current *this_j);
  void strict_motion_weighting(Time *time1,
                               Current *j1,
                               double x1_new,
                               double x3_new,
                               double x1_old,
                               double x3_old,
                               int p_number);
  void back_position_to_rz();

  void back_velocity_to_rz();

  void dump_position_to_old();

  void boris_pusher(EField *e_fld, HField *h_fld, Time *t);
  virtual void reflection_single(unsigned int i);
  void move_half_reflect(Time *t);

  void full_j_weighting(Current *current, Time *t);

// protected:
  void half_step_pos_single(Time *t, unsigned int i);
  void back_position_to_rz_single(unsigned int i);
  void dump_position_to_old_single(unsigned int i);
  void back_velocity_to_rz_single(unsigned int i);
  void step_v_single(EField *e_fld, HField *h_fld,
                     Time *t, unsigned int i);
};
