#include <iostream>
#include <hdf5.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Missing file name" << std::endl;
        return 1;
    }

    hid_t file_id = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Fclose(file_id);

    std::cout << "Status: " << status << std::endl;

    return 0;
}
