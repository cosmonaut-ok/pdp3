#include "input_output_class.h"
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

input_output_class::input_output_class(void)
{
}

input_output_class::input_output_class(char* c_pathres,char*c_pathdump)
{
  strcpy(path_result,c_pathres);
  strcpy(path_dump,c_pathdump);
}


input_output_class::~input_output_class(void)
{
}

void input_output_class::out_data(char const *comp_name,
                                  double** out_value,
                                  int step_number,
                                  int number,
                                  int r_step,
                                  int z_step)
{
  //static variable for defining number of function calling
  //static int count =0;
  //static variable. increase text file name
  int inc_value = 0;
  //  if (step_number/number)
  //{
  //  i=i+1;
  ////char operations
  //char st_name[50];
  //strcpy(st_name,comp_name);
  //char str_int [50];
  //itoa(i,str_int,10);
  //char path[50] = "e:/Plazma/pdp3_result/";
  //
  //strcat(path,st_name);
  //ofstream out_val(path,ios::app);
  //}

  inc_value = step_number/number;
  char st_name[100];
  strcpy(st_name,comp_name);
  char str_int [100];
// (inc_value,str_int,10);
  sprintf(str_int,"%d",inc_value);
  strcat(st_name, str_int);
  char path[100] = "E:/Science_plasma/pdp3_result_____/____"; // TODO: WTF
  //char path[50] = "D:/pdp3 files/_results/";
  strcpy(path,this->path_result);
  strcat(path,st_name);
  ofstream out_val(path,ios::app);
  // write  values  into file
  for (int i=0; i<r_step;i++)
  {
    for(int k=0;k<z_step;k++)
    {
      out_val<<out_value[i][k]<<" ";
    }
  }
  out_val.close();
}
void input_output_class::out_field_dump(char* comp_name,
                                        double** out_value,
                                        int r_step,
                                        int z_step)
{
  char st_name[100];

  strcpy(st_name,this->path_dump);

  //char path[50] = "D:/pdp3 files/_results/";
  strcat(st_name,comp_name);
  ofstream out_val(st_name);
  // write  values  into file
  for (int i=0; i<r_step;i++)
  {
    for(int k=0;k<z_step;k++)
    {
      out_val<<out_value[i][k]<<" ";
    }
  }
  out_val.close();
}

/////
void input_output_class::out_coord(char* comp_name,
                                   double* coord_r,
                                   double* coord_z,
                                   int step_number,
                                   int number,
                                   int particles_number)
{
  //static variable for defining number of function calling
  //static int count =0;
  //static variable. increase text file name
  int inc_value = 0;
  /*if (step_number/number)
    {
    i=i+1;
    char st_name[50];
    strcpy(st_name,comp_name);
    char str_int [50];
    itoa(i,str_int,10);
    char path[50] = "e:/Plazma/pdp3_result/";

    strcat(path,st_name);
    ofstream out_val(path,ios::app);
    }*/
  inc_value = step_number/number;
  char st_name[50];
  strcpy(st_name,comp_name);
  char str_int [50];
  //itoa(inc_value,str_int,10);
  sprintf(str_int,"%d",inc_value);
  strcat(st_name, str_int);
  char path[100];
  //char path[50] = "D:/pdp3 files/_results/";
  strcpy(path,this->path_result);
  strcat(path,st_name);
  ofstream out_val(path,ios::app);
  // write  values  into file
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)

  {
    out_val<<coord_r[i]<<" ";
    out_val<<coord_z[i]<<" ";
  }

  out_val.close();
}

void input_output_class::out_coord_dump(char* comp_name ,
                                        double* coord_r,
                                        double* coord_z,
                                        int particles_number)
{
  char st_name[100] ;

  strcat(strcpy(st_name,this->path_dump),"_coords_");

  //char path[50] = "D:/pdp3 files/_results/";
  strcat(st_name,comp_name);
  ofstream out_val(st_name);
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)
  {
    out_val<<coord_r[i]<<" ";
    out_val<<coord_z[i]<<" ";
  }
  out_val.close();
}

void input_output_class:: out_velocity_dump(char* comp_name,
                                            double* v1,
                                            double* v2,
                                            double* v3,
                                            int particles_number)
{
  char st_name[100] ;

  strcat(strcpy(st_name,this->path_dump),"_velocities_");

  //char path[50] = "D:/pdp3 files/_results/";
  strcat(st_name,comp_name);
  ofstream out_val(st_name);
  out_val.setf(std::ios_base::scientific);
  out_val.precision(14);
  for (int i=0; i<particles_number;i++)
  {
    out_val<<v1[i]<<" ";
    out_val<<v2[i]<<" ";
    out_val<<v3[i]<<" ";
  }

  out_val.close();

}
