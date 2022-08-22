#include <iostream>

#include <cvodes/cvodes.h>

int main()
{
    printf("Sundials installed. CV_SUCCESS: %d, CV_TSTOP_RETURN: %d, CV_ROOT_RETURN: %d.\n", CV_SUCCESS, CV_TSTOP_RETURN, CV_ROOT_RETURN);
}