#include "ioText.h"

IOText::IOText(void)
{
}

IOText::IOText(char *c_pathres,char *c_pathdump, bool c_compress)
{
  strcpy(path_result, c_pathres);
  strcpy(path_dump, c_pathdump);
  if (c_compress)
    cerr << "WARNING! compression is not supported by this backend" << endl;
}


IOText::~IOText(void)
{
}

char* IOText::get_data_path(char const *comp_name, int number=-1)
{
  char st_name[100];
  strcpy(st_name, comp_name);

  if (number != -1)
  {
    char str_int [100];
    sprintf(str_int, "%d", number);
    strcat(st_name, str_int);
  }

  char* path = new char[100];
  strcpy(path, path_result);
  strcat(path, st_name);

  return path;
}

char* IOText::get_dump_path(char const *comp_name)
{
  char st_name[100];
  strcpy(st_name, comp_name);
  char* path = new char[100];
  strcpy(path, path_dump);
  strcat(path, st_name);

  return path;
}

void IOText::dump_data(char const *path,
                                 double **out_value,
                                 int r_step,
                                 int z_step,
                                 bool is_rewrite=false)
{
  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_val(path, omode);

  // write  values  into file
  for (int i=0; i<r_step;i++)
    for(int k=0;k<z_step;k++)
      out_val<<out_value[i][k]<<" ";

  out_val.close();
}

void IOText::dump_components(char const *path,
                                       double **components,
                                       int size,
                                       int dimensions,
                                       bool is_rewrite=false)
{
  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_val(path, omode);
  // write  values  into file
  for (int i=0; i < size; i++)
  {
    for(int k=0; k < dimensions; k++)
      out_val << components[i][k] << ",";
    out_val << " ";
  }

  out_val.close();
}

void IOText::out_data(char const *comp_name,
                                  double **out_value,
                                  int step_number,
                                  int number,
                                  int r_step,
                                  int z_step)
{
  char* path = get_data_path(comp_name, floor(step_number/number));
  dump_data(path, out_value, r_step, z_step);

  delete [] path;
}

void IOText::out_field_dump(char const *comp_name,
                                        double **out_value,
                                        int r_step,
                                        int z_step)
{
  char* path = get_dump_path(comp_name);
  dump_data(path, out_value, r_step, z_step, true);

  delete [] path;
}

void IOText::out_triple(char *comp_name,
                                  double **triple,
                                  int step_number,
                                  int number,
                                  int particles_number)
{
  char* path = get_data_path(comp_name, floor(step_number/number));
  dump_components(path, triple, particles_number, 3);

  delete [] path;
}

void IOText::out_pos_dump(char *comp_name,
                                      double **pos,
                                      int particles_number)
{
  char* path = get_dump_path(comp_name);
  strcat(path, "_position");
  dump_components(path, pos, particles_number, 3, true);

  delete [] path;
}

void IOText::out_velocity_dump(char *comp_name,
                                         double **vel,
                                         int particles_number)
{
  char* path = get_dump_path(comp_name);
  strcat(path, "_velocity");
  dump_components(path, vel, particles_number, 3, true);

  delete [] path;
}

void IOText::out_current_time_dump(double current_time)
{
  char* path = get_dump_path("current_time");
  ofstream out_val(path, ios::trunc);
  out_val<<current_time;
  out_val.close();

  delete [] path;
}
