#include "lib.h"
#include "writerPlain.h"

using namespace std;

WriterPlain::WriterPlain(void)
{
}

WriterPlain::WriterPlain(char *c_path, char *c_component, int c_type,
                         int c_start_r, int c_start_z, int c_end_r, int c_end_z,
                         bool c_compress, int c_compress_level, int c_schedule)
{
  strcpy(path, c_path);
  strcpy(component, c_component);

  type = c_type;
  start_r = c_start_r;
  start_z = c_start_z;
  end_r = c_end_r;
  end_z = c_end_z;

  schedule = c_schedule;

  // not implemented
  compress = c_compress;
  compress_level = c_compress_level;
  if (compress)
    cerr << "WARNING! compression is not supported by plaintext backend" << endl;

  switch (type)
    {
    case 0:
      sprintf(path_result, "%s/%s/frame_%d:%d_%d:%d", c_path, c_component,
              c_start_r, c_start_z, c_end_r, c_end_z);
      break;
    case 1:
      sprintf(path_result, "%s/%s/col_%d", c_path, c_component, c_start_z);
      break;
    case 2:
      sprintf(path_result, "%s/%s/row_%d", c_path, c_component, c_start_r);
      break;
    case 3:
      sprintf(path_result, "%s/%s/dot_%d_%d", c_path, c_component,
              c_start_r, c_start_z);
      break;
    }

  lib::makeDirectory(path_result);
}

WriterPlain::~WriterPlain(void)
{
}

char* WriterPlain::get_data_path(char *name)
{
  char* path = new char[100];
  strcpy(path, path_result);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".dat");

  return path;
}

void WriterPlain::write_frame(char *name, double **out_value, bool is_rewrite = false)
{
  char *path = get_data_path(name);

  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_file(path, omode);

  // write  values  into file
  for (int i = start_r; i < end_r; i++)
    for(int k = start_z; k < end_z; k++)
      out_file<<out_value[i][k]<<" ";

  out_file.close();
}

void WriterPlain::write_col(char *name, double **out_value, bool is_rewrite = false)
{
  char *path = get_data_path(name);

  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_file(path, omode);

  // write  values  into file
  for(int k = start_r; k < end_r; k++)
    out_file<<out_value[k][start_z]<<" ";

  out_file.close();
}

void WriterPlain::write_row(char *name, double **out_value, bool is_rewrite = false)
{
  char *path = get_data_path(name);

  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_file(path, omode);

  // write  values  into file
  for(int k = start_z; k < end_z; k++)
    out_file<<out_value[start_r][k]<<" ";

  out_file.close();
}

void WriterPlain::write_dot(char *name, double **out_value, bool is_rewrite = false)
{
  char *path = get_data_path(name);

  std::ofstream::openmode omode;
  if (is_rewrite)
    omode = ios::trunc;
  else
    omode = ios::app;

  ofstream out_file(path, omode);

  // write  values  into file
  out_file<<out_value[start_r][start_z]<<" ";

  out_file.close();
}

void WriterPlain::write(char *name, double **out_value)
{
  switch (type)
    {
    case 0:
      write_frame(name, out_value, false);
      break;
    case 1:
      write_col(name, out_value, false);
      break;
    case 2:
      write_row(name, out_value, false);
      break;
    case 3:
      write_dot(name, out_value, false);
      break;
    }
}

// char* WriterPlain::get_data_path(char const *comp_name, int number=-1)
// {
//   char st_name[100];
//   strcpy(st_name, comp_name);

//   if (number != -1)
//   {
//     char str_int [100];
//     sprintf(str_int, "%d", number);
//     strcat(st_name, str_int);
//   }

//   char* path = new char[100];
//   strcpy(path, path_result);
//   strcat(path, st_name);

//   return path;
// }

// char* WriterPlain::get_dump_path(char const *comp_name)
// {
//   char st_name[100];
//   strcpy(st_name, comp_name);
//   char* path = new char[100];
//   strcpy(path, path_dump);
//   strcat(path, st_name);

//   return path;
// }

// void WriterPlain::dump_data(char const *path,
//                                  double **out_value,
//                                  int r_step,
//                                  int z_step,
//                                  bool is_rewrite=false)
// {
//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_val(path, omode);

//   // write  values  into file
//   for (int i=0; i<r_step;i++)
//     for(int k=0;k<z_step;k++)
//       out_val<<out_value[i][k]<<" ";

//   out_val.close();
// }

// void WriterPlain::dump_components(char const *path,
//                                        double **components,
//                                        int size,
//                                        int dimensions,
//                                        bool is_rewrite=false)
// {
//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_val(path, omode);
//   // write  values  into file
//   for (int i=0; i < size; i++)
//   {
//     for(int k=0; k < dimensions; k++)
//       out_val << components[i][k] << ",";
//     out_val << " ";
//   }

//   out_val.close();
// }

// void WriterPlain::out_data(char const *comp_name,
//                                   double **out_value,
//                                   int step_number,
//                                   int number,
//                                   int r_step,
//                                   int z_step)
// {
//   char* path = get_data_path(comp_name, floor(step_number/number));
//   dump_data(path, out_value, r_step, z_step);

//   delete [] path;
// }

// void WriterPlain::out_field_dump(char const *comp_name,
//                                         double **out_value,
//                                         int r_step,
//                                         int z_step)
// {
//   char* path = get_dump_path(comp_name);
//   dump_data(path, out_value, r_step, z_step, true);

//   delete [] path;
// }

// void WriterPlain::out_triple(char *comp_name,
//                                   double **triple,
//                                   int step_number,
//                                   int number,
//                                   int particles_number)
// {
//   char* path = get_data_path(comp_name, floor(step_number/number));
//   dump_components(path, triple, particles_number, 3);

//   delete [] path;
// }

// void WriterPlain::out_pos_dump(char *comp_name,
//                                       double **pos,
//                                       int particles_number)
// {
//   char* path = get_dump_path(comp_name);
//   strcat(path, "_position");
//   dump_components(path, pos, particles_number, 3, true);

//   delete [] path;
// }

// void WriterPlain::out_velocity_dump(char *comp_name,
//                                          double **vel,
//                                          int particles_number)
// {
//   char* path = get_dump_path(comp_name);
//   strcat(path, "_velocity");
//   dump_components(path, vel, particles_number, 3, true);

//   delete [] path;
// }

// void WriterPlain::out_current_time_dump(double current_time)
// {
//   char* path = get_dump_path("current_time");
//   ofstream out_val(path, ios::trunc);
//   out_val<<current_time;
//   out_val.close();

//   delete [] path;
// }
