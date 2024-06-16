#include <PythonFacade.h>
#include <Bathymetry.h>

void Facade::initDomain(
        const Storage<3>& geometry,
        const Topology& topology) {
    this->_bathymetry.reset(new Bathymetry(mesh));
}

void Facade::initFields(const Storage<3>& initialValues) {
    
}

    void initSpaceDisc(
            const SpaceDisc::Type discType,
            const double coriolisForce,
            const double frictionForce);

    void initTimeDisc(
            const TimeDisc::Type discType);

    Storage<3> run(const double time);
}



