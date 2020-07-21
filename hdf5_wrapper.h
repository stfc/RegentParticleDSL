#include <hdf5.h>


hid_t WRAP_H5F_ACC_TRUNC = 2;
hid_t WRAP_H5T_STD_I32LE = 0;
hid_t WRAP_H5T_STD_I64LE = 0;
hid_t WRAP_H5T_IEEE_F64LE = 0;
hid_t WRAP_H5P_DEFAULT = 0;
hid_t WRAP_H5S_SIMPLE = 0;
hid_t WRAP_H5S_ALL = 0;
hid_t WRAP_H5T_NATIVE_DOUBLE = 0;
hid_t WRAP_H5T_NATIVE_INT = 0;
hid_t WRAP_H5F_ACC_RDONLY = 0;

void init_wrapper(){
 WRAP_H5F_ACC_TRUNC = H5F_ACC_TRUNC;
 WRAP_H5T_STD_I32LE = H5T_STD_I32LE;
 WRAP_H5T_STD_I64LE = H5T_STD_I64LE;
 WRAP_H5T_IEEE_F64LE = H5T_IEEE_F64LE;
 WRAP_H5P_DEFAULT = H5P_DEFAULT;
 WRAP_H5S_SIMPLE = H5S_SIMPLE;
 WRAP_H5S_ALL = H5S_ALL;
 WRAP_H5T_NATIVE_DOUBLE = H5T_NATIVE_DOUBLE;
 WRAP_H5T_NATIVE_INT = H5T_NATIVE_INT;
 WRAP_H5F_ACC_RDONLY = H5F_ACC_RDONLY;
}
