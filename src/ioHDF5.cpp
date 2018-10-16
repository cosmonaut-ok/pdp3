#include "ioHDF5.h"

IOHDF5::IOHDF5(void)
{
}

IOHDF5::IOHDF5(char *c_pathres, char *c_pathdump, bool c_compress)
{
  hid_t file_id, gcpl;

  compress = c_compress;

  strcpy(path_result, c_pathres);
  strcpy(path_dump, c_pathdump);

  strcpy(path_data, c_pathres);

  strcpy(hdf5_file, path_data);
  strcat(hdf5_file, "data.h5");

  strcpy(root_group_name,  "/pdp3");

  sprintf(result_group_name, "%s/%s", root_group_name, "result");
  sprintf(state_group_name, "%s/%s", root_group_name, "dump");

  // Create a new file using the default properties.
  file_id = H5Fcreate (hdf5_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Create group creation property list and set it to allow creation
  // of intermediate groups.
  gcpl = H5Pcreate (H5P_LINK_CREATE);
  status = H5Pset_create_intermediate_group (gcpl, 1);

  // close
  status = H5Fclose(file_id);
}

IOHDF5::~IOHDF5(void)
{
  // status = H5Fclose(file_id);
}

hid_t IOHDF5::create_or_open_h5_group(hid_t file_id, char const *group_name)
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

void IOHDF5::dump_h5_dataset(char const *group_name,
                             char const *name,
                             double **out_value,
                             int r_step,
                             int z_step)
{
  hid_t file_id, group_id, dataset, datatype, dataspace, plist_id;   /* declare identifiers */
  hsize_t dimsf[2];              /* dataset dimensions */

  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", group_name, name);

  double dv[r_step][z_step];
  for (int i=0; i< r_step; i++)
    for (int j=0; j< z_step; j++)
      dv[i][j] = out_value[i][j];

  file_id = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);

  group_id = create_or_open_h5_group(file_id, group_name);

  dimsf[0] = r_step;
  dimsf[1] = z_step;

  dataspace = H5Screate_simple(2, dimsf, NULL);

  plist_id  = H5P_DEFAULT;
  if (compress)
  {
    plist_id  = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk (plist_id, 2, dimsf);
    status = H5Pset_deflate (plist_id, compress_level);
  }

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  status = H5Tset_order(datatype, H5T_ORDER_LE);

  // Save old error handler
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto1 (&old_func, &old_client_data);

  status = H5Eset_auto1(NULL, NULL);

  dataset = H5Dopen(file_id, name, H5P_DEFAULT);

  if (dataset == -1)
    dataset = H5Dcreate(file_id, name, datatype, dataspace,
                        H5P_DEFAULT, plist_id, H5P_DEFAULT);

  status = H5Eset_auto1(old_func, old_client_data);

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, dv);

  status = H5Sclose(dataspace);
  status = H5Tclose(datatype);
  status = H5Dclose(dataset);
  status = H5Gclose(group_id);
  status = H5Fclose(file_id);
}

void IOHDF5::dump_h5_nd_components(char const *group_name,
                                   char const *name,
                                   double **components,
                                   int size,
                                   int dimensions)
//! dump array from N-dimensional vectors
{
  hid_t file_id, group_id, dataset, datatype, dataspace;
  hsize_t dimsf[size];

  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", group_name, name);

  file_id = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);

  group_id = create_or_open_h5_group(file_id, group_name);

  // fill dimsf with arrays
  for (int i=0; i < dimensions; i++)
    dimsf[i] = dimensions;

  dataspace = H5Screate_simple(2, dimsf, NULL);
  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  status = H5Tset_order(datatype, H5T_ORDER_LE);

  // Save old error handler
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto1 (&old_func, &old_client_data);

  status = H5Eset_auto1(NULL, NULL);

  dataset = H5Dopen(file_id, name, H5P_DEFAULT);

  if (dataset == -1)
    dataset = H5Dcreate(file_id, name, datatype, dataspace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Eset_auto1(old_func, old_client_data);

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, &components);

  status = H5Sclose(dataspace);
  status = H5Tclose(datatype);
  status = H5Dclose(dataset);
  status = H5Gclose(group_id);
}

void IOHDF5::out_data(char const *comp_name,
                                  double **out_value,
                                  int step_number,
                                  int number,
                                  int r_step,
                                  int z_step)
{
  char subgroup_name[200];
  char dataset_name[200];
  sprintf(subgroup_name, "%s/%s", result_group_name, comp_name);
  sprintf(dataset_name, "%s/%d", subgroup_name, step_number);
  (void)(number); // just avoid unused argument warning

  dump_h5_dataset(subgroup_name, dataset_name, out_value, r_step, z_step);
}

void IOHDF5::out_field_dump(char const *comp_name,
                                        double **out_value,
                                        int r_step,
                                        int z_step)
{
  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", state_group_name, comp_name);

  dump_h5_dataset(state_group_name, dataset_name, out_value, r_step, z_step);
}

void IOHDF5::out_triple(char *comp_name,
                        double **triple,
                        int step_number,
                        int number,
                        int particles_number)
{
  out_data(comp_name, triple, step_number, number, particles_number, 3);
}

void IOHDF5::out_pos_dump(char *comp_name,
                                      double **pos,
                                      int particles_number)
{
  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", state_group_name, comp_name);
  strcat(dataset_name, "_position");

  dump_h5_nd_components(state_group_name, dataset_name, pos, particles_number, 3);
}

void IOHDF5::out_velocity_dump(char *comp_name,
                                         double **vel,
                                         int particles_number)
{
  char dataset_name[200];
  sprintf(dataset_name, "%s/%s", state_group_name, comp_name);
  strcat(dataset_name, "_velocity");

  dump_h5_nd_components(state_group_name, dataset_name, vel, particles_number, 3);
}

void IOHDF5::out_current_time_dump(double current_time)
{
  hid_t file_id, attribute_id, dataspace_id, group_id;
  hsize_t dims;
  double attr_data[1];

  file_id = H5Fopen (hdf5_file, H5F_ACC_RDWR, H5P_DEFAULT);

  group_id = create_or_open_h5_group(file_id, state_group_name);

  // Initialize the attribute data
  attr_data[0] = current_time;

  // Create the data space for the attribute
  dims = 1;
  dataspace_id = H5Screate_simple(1, &dims, NULL);

  // Save old error handler
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto1 (&old_func, &old_client_data);

  status = H5Eset_auto1(NULL, NULL);

  attribute_id = H5Aopen (group_id, "timestamp", H5P_DEFAULT);

  /* Create a dataset attribute. */
  if (attribute_id == -1)
    attribute_id = H5Acreate (group_id, "timestamp", H5T_NATIVE_DOUBLE, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT);

  status = H5Eset_auto1(old_func, old_client_data);

  /* Write the attribute data. */
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, attr_data);

  status = H5Aclose(attribute_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
}
