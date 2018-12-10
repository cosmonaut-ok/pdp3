// follow this: https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_extend.c
#include "probeHDF5.h"

using namespace std;

ProbeHDF5::ProbeHDF5(void)
{
}

ProbeHDF5::ProbeHDF5(char *c_path, char *c_data_root, char *c_component, int c_type,
                     int c_start_r, int c_start_z, int c_end_r, int c_end_z,
                     bool c_compress, int c_compress_level, int c_schedule,
                     unsigned int c_r_size, unsigned int c_z_size)
{

  strcpy(path_result, c_path);
  strcpy(data_root, c_data_root);
  strcpy(component, c_component);

  type = c_type;
  start_r = c_start_r;
  start_z = c_start_z;
  end_r = c_end_r;
  end_z = c_end_z;
  schedule = c_schedule;
  compress = c_compress;
  compress_level = c_compress_level;

  // make hdf5 file path
  strcpy(hdf5_file, path_result);
  strcat(hdf5_file, "/data.h5"); // yes, baby, it's hardcode

  unsigned int r_rank;

  // set probe path and data set name
  switch (type)
  {
  case 0: // frame
    sprintf(probe_path, "%s/%s/frame", data_root, component);
    sprintf(probe_data_set_name, "%d:%d_%d:%d", start_r, start_z, end_r, end_z);
    // set dataset rank
    r_rank = 3;
    break;
  case 1: // column
    sprintf(probe_path, "%s/%s/col", data_root, component);
    sprintf(probe_data_set_name, "%d", start_z);
    // set dataset rank
    r_rank = 2;
    break;
  case 2: // row
    sprintf(probe_path, "%s/%s/row", data_root, component);
    sprintf(probe_data_set_name, "%d", start_r);
    // set dataset rank
    r_rank = 2;
    break;
  case 3: // dot
    sprintf(probe_path, "%s/%s/dot", data_root, component);
    sprintf(probe_data_set_name, "%d_%d", start_r, start_z);
    // set dataset rank
    r_rank = 1;
    break;
  case 4: // macroparticles frame
    sprintf(probe_path, "%s/%s/mpframe", data_root, component);
    sprintf(probe_data_set_name, "%d:%d_%d:%d", start_r, start_z, end_r, end_z);
    // set dataset rank
    r_rank = 3;
    break;
  }

  char dataset_path[200]; // will contain full dataset path
  sprintf(dataset_path, "%s/%s", probe_path, probe_data_set_name);

////////////////////////////////////////////////////////////////////////////////
  hid_t file;
  hid_t dataspace, dataset;
  hid_t filespace, memspace;
  hid_t prop;
  hid_t group;

  hsize_t dims[r_rank]; // = {3, 3, 3}; /* dataset dimensions at creation time */
  hsize_t maxdims[r_rank]; // = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED};
  herr_t status;
  hsize_t chunk_dims[r_rank]; //  = {2, 2, 5};

  switch (type)
  {
  case 0: // frame
    dims[0] = end_r - start_r;
    dims[1] = end_z - start_z;
    dims[2] = 1;
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = H5S_UNLIMITED;
    maxdims[2] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    chunk_dims[1] = 2;
    chunk_dims[2] = 2;
    break;
  case 1: // col
    dims[0] = c_z_size;
    dims[1] = 1;
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    chunk_dims[1] = 2;
    break;
  case 2: // row
    dims[0] = c_r_size;
    dims[1] = 1;
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    chunk_dims[1] = 2;
    break;
  case 3: // dot
    dims[0] = 1;
    maxdims[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    break;
  case 4: // mpframe
    dims[0] = end_r - start_r;
    dims[1] = end_z - start_z;
    dims[2] = 1;
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = H5S_UNLIMITED;
    maxdims[2] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    chunk_dims[1] = 1;
    chunk_dims[2] = 1;
    break;
  }

  /* Create the data space with unlimited dimensions. */
  dataspace = H5Screate_simple (r_rank, dims, maxdims);

  //// disable error printing
  // Save old error handler
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto1 (&old_func, &old_client_data);

  status = H5Eset_auto1(NULL, NULL);
  /* Create a new file. If file exists its contents will be overwritten. */
  file = H5Fcreate (hdf5_file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  //// enable error printing
  status = H5Eset_auto1(old_func, old_client_data);

  if (file < 0)
    file = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);
  else
  {
    // Create group creation property list and set it to allow creation
    // of intermediate groups.
    hid_t gcpl = H5Pcreate (H5P_LINK_CREATE);
    status = H5Pset_create_intermediate_group (gcpl, 1);
  }

  // catch HDF file creating/opening error
  if (file < 0)
  {
    cerr << "ERROR: can not create/open HDF5 data file ``" << hdf5_file << "''. Exiting." << endl;
    exit(1);
  }

  group = create_or_open_h5_group(file, probe_path);

  /* Modify dataset creation properties, i.e. enable chunking  */
  prop = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_chunk (prop, r_rank, chunk_dims);

  /* Create a new dataset within the file using chunk
     creation properties.  */
  dataset = H5Dcreate2 (file, dataset_path, H5T_NATIVE_INT, dataspace,
                        H5P_DEFAULT, prop, H5P_DEFAULT);


  status = H5Pclose (prop);
  status = H5Dclose (dataset);
  status = H5Gclose (group);
  // status = H5Sclose (filespace);
  // status = H5Sclose (memspace);
  status = H5Fclose (file);
}

ProbeHDF5::~ProbeHDF5(void)
{
}

hid_t ProbeHDF5::create_or_open_h5_group(hid_t file_id, char const *group_name)
{
  herr_t status;
  hid_t group_id, gcpl;

  // Save old error handler
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto1 (&old_func, &old_client_data);

  gcpl = H5Pcreate (H5P_LINK_CREATE);
  status = H5Pset_create_intermediate_group (gcpl, 1);
  status = H5Eset_auto1(NULL, NULL);
  status = H5Gget_objinfo (file_id, group_name, 0, NULL);

  if (status != 0)
  { // create group if not exists
#ifdef PERF_DEBUG
    cout << "create group " << group_name << endl;
#endif
    group_id = H5Gcreate (file_id, group_name, gcpl, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
#ifdef PERF_DEBUG
    cout << "open group " << group_name << endl;
#endif
    group_id = H5Gopen (file_id, group_name, H5P_DEFAULT);
  }

  status = H5Eset_auto1(old_func, old_client_data);

  return(group_id);
}

hid_t ProbeHDF5::create_or_open_h5_dataset(
  hid_t file_id, char const *path,
  int rank, hsize_t *dims)
{
  // hid_t datatype, dataspace, dataset;
  // hid_t filespace, memspace;
  // hid_t prop;
  // herr_t status;

  // // Create the data space with unlimited dimensions
  // // int rank = rank;
  // // hsize_t h_dims[rank] = dims;
  // // h_dims = dims; /* dataset dimensions at creation time */
  // hsize_t maxdims[rank] = {H5S_UNLIMITED, H5S_UNLIMITED};
  // dataspace = H5Screate_simple (rank, dims, maxdims);

  // datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  // status = H5Tset_order(datatype, H5T_ORDER_LE);

  // // Modify dataset creation properties, i.e. enable chunking
  // hsize_t chunk_dims[2] = {2, 5}; // TODO: WTF
  // prop = H5Pcreate (H5P_DATASET_CREATE);
  // status = H5Pset_chunk (prop, rank, chunk_dims);

  // //
  // // create or open dataset
  // //
  // // Save old error handler
  // herr_t (*old_func)(void*);
  // void *old_client_data;
  // H5Eget_auto1 (&old_func, &old_client_data);

  // status = H5Eset_auto1(NULL, NULL);

  // dataset = H5Dopen2 (file_id, path, H5P_DEFAULT);

  // if (dataset == -1)
  //   dataset = H5Dcreate2 (file_id, path, datatype, dataspace,
  //                         H5P_DEFAULT, prop, H5P_DEFAULT);

  // // cout << "create_or_open_h5_dataset " << dataset << endl;
  // status = H5Eset_auto1(old_func, old_client_data);

  // return(dataset);
}

////////////////////////////////////////////////////////////////////////////////



  // int data[3][3] = { {1, 1, 1},    /* data to write */
  //                    {1, 1, 1},
  //                    {1, 1, 1} };

  // /* Variables used in extending and writing to the extended portion of dataset */
  // hsize_t size[2];
  // hsize_t offset[2];
  // hsize_t dimsext[2] = {7, 3};         /* extend dimensions */
  // int dataext[7][3] = { {2, 3, 4},
  //                       {2, 3, 4},
  //                       {2, 3, 4},
  //                       {2, 3, 4},
  //                       {2, 3, 4},
  //                       {2, 3, 4},
  //                       {2, 3, 4} };

  // /* Variables used in reading data back */
  // hsize_t chunk_dimsr[2];
  // hsize_t dimsr[2];
  // hsize_t i, j;
  // int rdata[10][3];
  // herr_t status_n;
  // int rank, rank_chunk;

////
  // /* Create a new dataset within the file using chunk
  //    creation properties.  */
  // dataset = H5Dcreate2 (file, DATASETNAME, H5T_NATIVE_INT, dataspace,
  //                      H5P_DEFAULT, prop, H5P_DEFAULT);

  // /* Write data to dataset */
  // status = H5Dwrite (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
  //                    H5P_DEFAULT, data);

  // /* Extend the dataset. Dataset becomes 10 x 3  */
  // size[0] = dims[0]+ dimsext[0];
  // size[1] = dims[1];
  // status = H5Dset_extent (dataset, size);

  // /* Select a hyperslab in extended portion of dataset  */
  // filespace = H5Dget_space (dataset);
  // offset[0] = 3;
  // offset[1] = 0;
  // status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
  //                               dimsext, NULL);

  // /* Define memory space */
  // memspace = H5Screate_simple (RANK, dimsext, NULL);

  // /* Write the data to the extended portion of dataset  */
  // status = H5Dwrite (dataset, H5T_NATIVE_INT, memspace, filespace,
  //                    H5P_DEFAULT, dataext);

  // /* Close resources */
  // status = H5Dclose (dataset);
  // status = H5Pclose (prop);
  // status = H5Sclose (dataspace);
  // status = H5Sclose (memspace);
  // status = H5Sclose (filespace);
  // status = H5Fclose (file);

  // /********************************************
  //  * Re-open the file and read the data back. *
  //  ********************************************/

  // file = H5Fopen (FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  // dataset = H5Dopen2 (file, DATASETNAME, H5P_DEFAULT);

  // filespace = H5Dget_space (dataset);
  // rank = H5Sget_simple_extent_ndims (filespace);
  // status_n = H5Sget_simple_extent_dims (filespace, dimsr, NULL);

  // prop = H5Dget_create_plist (dataset);

  // if (H5D_CHUNKED == H5Pget_layout (prop))
  //    rank_chunk = H5Pget_chunk (prop, rank, chunk_dimsr);

  // memspace = H5Screate_simple (rank, dimsr, NULL);
  // status = H5Dread (dataset, H5T_NATIVE_INT, memspace, filespace,
  //                   H5P_DEFAULT, rdata);

  // printf("\n");
  // printf("Dataset: \n");
  // for (j = 0; j < dimsr[0]; j++)
  // {
  //    for (i = 0; i < dimsr[1]; i++)
  //        printf("%d ", rdata[j][i]);
  //    printf("\n");
  // }

  // status = H5Pclose (prop);
  // status = H5Dclose (dataset);
  // status = H5Sclose (filespace);
  // status = H5Sclose (memspace);
  // status = H5Fclose (file);

// }

void ProbeHDF5::write_frame(char *name, double **out_value, bool is_rewrite = false)
{
  cout << "write" << endl;
  herr_t status;
  hid_t file_id, group_id, ds_id;

  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", probe_path, probe_data_set_name);

  file_id = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);

  int rank = 2;
  hsize_t dims[rank] = {254, 2046};

  group_id = create_or_open_h5_group(file_id, probe_path);
  ds_id = create_or_open_h5_dataset(file_id, dataset_name, rank, dims);

  double data[3][3] = { {1, 1, 1},    /* data to write */
                        {1, 1, 1},
                        {1, 1, 1} };



  cout << "AAAAA :" << ds_id << endl;
  status = H5Dwrite(ds_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, out_value);

  cout << "BBBBB :" << status << endl;
  status = H5Dclose(ds_id);
  cout << "CCCCC :" << status << endl;
  status = H5Gclose(group_id);
  cout << "DDDDD :" << status << endl;
  status = H5Fclose(file_id);
  cout << "EEEEE :" << status << endl;

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
}

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
  cout << "/me dummy" << endl;
  // switch (type)
  // {
  // case 0:
  //   write_frame(name, out_value, false);
  //   break;
//     case 1:
//       write_col(name, out_value, false);
//       break;
//     case 2:
//       write_row(name, out_value, false);
//       break;
//     case 3:
//       write_dot(name, out_value, false);
//       break;
// }
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
