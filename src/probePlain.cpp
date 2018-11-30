#include "lib.h"
#include "probePlain.h"

using namespace std;

ProbePlain::ProbePlain(void)
{
}

ProbePlain::ProbePlain(char *c_path, char *c_component, int c_type,
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

ProbePlain::~ProbePlain(void)
{
}

char* ProbePlain::get_data_path(char *name)
{
  char* path = new char[100];
  strcpy(path, path_result);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".dat");

  return path;
}

void ProbePlain::write_frame(char *name, double **out_value, bool is_rewrite = false)
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

void ProbePlain::write_col(char *name, double **out_value, bool is_rewrite = false)
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

void ProbePlain::write_row(char *name, double **out_value, bool is_rewrite = false)
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

void ProbePlain::write_dot(char *name, double **out_value, bool is_rewrite = false)
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

void ProbePlain::write(char *name, double **out_value)
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
