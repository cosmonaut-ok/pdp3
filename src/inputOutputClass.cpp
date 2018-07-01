#include "inputOutputClass.h"

InputOutputClass::InputOutputClass(void)
{
}

InputOutputClass::InputOutputClass(char *c_pathres,char *c_pathdump)
{
  strcpy(path_result, c_pathres);
  strcpy(path_dump, c_pathdump);
}


InputOutputClass::~InputOutputClass(void)
{
}

void InputOutputClass::out_data(char const *comp_name,
                                  double **out_value,
                                  int step_number,
                                  int number,
                                  int r_step,
                                  int z_step)
{
  int inc_value = step_number/number;
  char st_name[100];
  strcpy(st_name, comp_name);
  char str_int [100];

  sprintf(str_int, "%d", inc_value);
  strcat(st_name, str_int);
  char path[100];
  strcpy(path, this->path_result);
  strcat(path, st_name);
  ofstream out_val(path,ios::app);
  // write  values  into file
  for (int i=0; i<r_step;i++)
    for(int k=0;k<z_step;k++)
      out_val<<out_value[i][k]<<" ";

  out_val.close();
}

void InputOutputClass::out_field_dump(char *comp_name,
                                        double **out_value,
                                        int r_step,
                                        int z_step)
{
  char st_name[100];

  strcpy(st_name, this->path_dump);

  strcat(st_name, comp_name);
  ofstream out_val(st_name);
  // write  values  into file
  for (int i=0; i<r_step;i++)
    for(int k=0;k<z_step;k++)
      out_val<<out_value[i][k]<<" ";

  out_val.close();
}

/////
void InputOutputClass::out_pos(char *comp_name,
                                   double *pos_r,
                                   double *pos_z,
                                   int step_number,
                                   int number,
                                   int particles_number)
{
  int inc_value = step_number/number;
  char st_name[50];
  strcpy(st_name, comp_name);
  char str_int [50];

  sprintf(str_int, "%d" ,inc_value);
  strcat(st_name, str_int);
  char path[100];

  strcpy(path, this->path_result);
  strcat(path,st_name);
  ofstream out_val(path, ios::app);

  // write  values  into file
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)
  {
    out_val<<pos_r[i]<<" ";
    out_val<<pos_z[i]<<" ";
  }

  out_val.close();
}

void InputOutputClass::out_pos_dump(char *comp_name ,
                                      double **pos,
                                      // double *pos_r,
                                      // double *pos_z,
                                      int particles_number)
{
  char st_name[100] ;

  strcat(strcpy(st_name,this->path_dump),"_poss_");

  strcat(st_name,comp_name);
  ofstream out_val(st_name);
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)
  {
    out_val<<pos[i][0]<<" ";
    out_val<<pos[i][2]<<" ";
  }
  out_val.close();
}

void InputOutputClass::out_velocity_dump(char *comp_name,
                                         double **vel,
                                         // double *v1,
                                         // double *v2,
                                         // double *v3,
                                         int particles_number)
{
  char st_name[100] ;

  strcat(strcpy(st_name,this->path_dump),"_velocities_");

  strcat(st_name, comp_name);
  ofstream out_val(st_name);
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)
  {
    for (int j=0; j<3; j++)
      out_val<<vel[i][j]<<" ";
    // out_val<<v1[i]<<" ";
    // out_val<<v2[i]<<" ";
    // out_val<<v3[i]<<" ";
  }

  out_val.close();

}
