#include "probeHDF5.h"

using namespace std;

ProbeHDF5::ProbeHDF5(void)
{
}

ProbeHDF5::ProbeHDF5(char *c_path, char *c_data_root, char *c_component, int c_type,
                         int c_start_r, int c_start_z, int c_end_r, int c_end_z,
                         bool c_compress, int c_compress_level, int c_schedule)
{
  hid_t file_id, gcpl;
  herr_t status;

  strcpy(path_result, c_path);
  strcpy(data_root, c_data_root);
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

  // make hdf5 file path
  strcpy(hdf5_file, path_result);
  strcat(hdf5_file, "/data.h5"); // yes, baby, it's hardcode

  #ifdef DEBUG
  cerr << "Creating directory " << path_result << endl;
#endif
  lib::makeDirectory(path_result);

  //////////////////////////////////////////////////////////////////////////////
  /* Save old error handler */
  herr_t (*old_func)(void*);
  void *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Probe. Likely to fail, but that’s okay */
  // status = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);
  // file_id = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);


  file_id = H5Fcreate (hdf5_file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  // file_id = H5Fcreate (hdf5_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //////////////////////////////////////////////////////////////////////////////

  // Create group creation property list and set it to allow creation
  // of intermediate groups.
  gcpl = H5Pcreate (H5P_LINK_CREATE);
  status = H5Pset_create_intermediate_group (gcpl, 1);

  // close
  status = H5Fclose(file_id);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  // switch (type)
  //   {
  //   case 0:
  //     sprintf(path_result, "%s/%s/frame_%d:%d_%d:%d", c_path, c_component,
  //             c_start_r, c_start_z, c_end_r, c_end_z);
  //     break;
  //   case 1:
  //     sprintf(path_result, "%s/%s/col_%d", c_path, c_component, c_start_z);
  //     break;
  //   case 2:
  //     sprintf(path_result, "%s/%s/row_%d", c_path, c_component, c_start_r);
  //     break;
  //   case 3:
  //     sprintf(path_result, "%s/%s/dot_%d_%d", c_path, c_component,
  //             c_start_r, c_start_z);
  //     break;
  //   case 4:
  //     sprintf(path_result, "%s/%s/mpframe_%d:%d_%d:%d", c_path, c_component,
  //             c_start_r, c_start_z, c_end_r, c_end_z);
  //     break;
  //   }
}

ProbeHDF5::~ProbeHDF5(void)
{
}

// char* ProbeHDF5::get_data_path(char *name)
// {
//   char* path = new char[100];
//   strcpy(path, path_result);
//   strcat(path, "/");
//   strcat(path, name);
//   strcat(path, ".dat");

//   return path;
// }

// void ProbeHDF5::write_frame(char *name, double **out_value, bool is_rewrite = false)
// {
//   char *path = get_data_path(name);

//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_file(path, omode);

//   // write  values  into file
//   for (int i = start_r; i < end_r; i++)
//     for(int k = start_z; k < end_z; k++)
//       out_file<<out_value[i][k]<<" ";
//   out_file.close();
// }

// void ProbeHDF5::write_col(char *name, double **out_value, bool is_rewrite = false)
// {
//   char *path = get_data_path(name);

//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_file(path, omode);

//   // write  values  into file
//   for(int k = start_r; k < end_r; k++)
//     out_file<<out_value[k][start_z]<<" ";
//   out_file.close();
// }

// void ProbeHDF5::write_row(char *name, double **out_value, bool is_rewrite = false)
// {
//   char *path = get_data_path(name);

//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_file(path, omode);

//   // write  values  into file
//   for(int k = start_z; k < end_z; k++)
//     out_file<<out_value[start_r][k]<<" ";

//   out_file.close();
// }

// void ProbeHDF5::write_dot(char *name, double **out_value, bool is_rewrite = false)
// {
//   char *path = get_data_path(name);

//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_file(path, omode);

//   // write  values  into file
//   out_file<<out_value[start_r][start_z]<<" ";

//   out_file.close();
// }

void ProbeHDF5::write(char *name, double **out_value)
{
  cout << "write" << endl;
//   switch (type)
//     {
//     case 0:
//       write_frame(name, out_value, false);
//       break;
//     case 1:
//       write_col(name, out_value, false);
//       break;
//     case 2:
//       write_row(name, out_value, false);
//       break;
//     case 3:
//       write_dot(name, out_value, false);
//       break;
//     }
}

void ProbeHDF5::mpwrite(char *name, Particles *p_specie)
{
  cout << "mpwrite" << endl;
//   bool is_rewrite = false;
//   char name_r[100];
//   char name_phi[100];
//   char name_z[100];
//   char name_pos_r[100];
//   char name_pos_z[100];

//   sprintf(name_r, "%s_vel_r", name);
//   sprintf(name_phi, "%s_vel_phi", name);
//   sprintf(name_z, "%s_vel_z", name);

//   sprintf(name_pos_r, "%s_pos_r", name);
//   sprintf(name_pos_z, "%s_pos_z", name);

//   char *path_r = get_data_path(name_r);
//   char *path_phi = get_data_path(name_phi);
//   char *path_z = get_data_path(name_z);
//   char *path_pos_r = get_data_path(name_pos_r);
//   char *path_pos_z = get_data_path(name_pos_z);

//   std::ofstream::openmode omode;
//   if (is_rewrite)
//     omode = ios::trunc;
//   else
//     omode = ios::app;

//   ofstream out_file_r(path_r, omode);
//   ofstream out_file_phi(path_phi, omode);
//   ofstream out_file_z(path_z, omode);
//   ofstream out_file_pos_r(path_pos_r, omode);
//   ofstream out_file_pos_z(path_pos_z, omode);

//   for (unsigned int i = 0; i < p_specie->number; i++)
//   {
//     double dr = p_specie->geom1->dr;
//     double dz = p_specie->geom1->dz;

//     unsigned int r_i = (int)ceil((p_specie->pos[i][0])/dr)-1;
//     unsigned int z_k = (int)ceil((p_specie->pos[i][2])/dz)-1;

//     if (r_i >= start_r && r_i <= end_r && z_k >= start_z && z_k <= end_z)
//     {
//       out_file_r<<p_specie->vel[i][0]<<" ";
//       out_file_phi<<p_specie->vel[i][1]<<" ";
//       out_file_z<<p_specie->vel[i][2]<<" ";
//       out_file_pos_r<<p_specie->pos[i][0]<<" ";
//       out_file_pos_z<<p_specie->pos[i][2]<<" ";
//     }
//   }
}
