#include <iostream>
#include <SAMRAI/SAMRAI_config.h>

int main() {
    printf("SAMRAI config working. Version: %d.%d.%d\n", (int) SAMRAI_VERSION_MAJOR, (int) SAMRAI_VERSION_MINOR, (int) SAMRAI_VERSION_PATCHLEVEL);
}
