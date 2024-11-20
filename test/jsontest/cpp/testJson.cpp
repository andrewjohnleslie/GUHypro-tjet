#include "core/systemModule.h"
#include "core/propModule.h"
// #include <sstream>

int main(int argc, char *argv[]) {
    // stringstream ss;
    std::string filename = "filename";
    std::shared_ptr<hypro::propModule> sys1 = hypro::propModule::loadJson(filename);
    sys1->calculate();

}