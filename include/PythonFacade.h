#include <Includes.h>
#include <SpaceDisc.h>
#include <TimeDisc.h>
#include <memory>

class Facade {

    std::unique_ptr<Bathymetry> _bathymetry;
    std::unique_ptr<SpaceDisc> _spaceDisc;
    std::unique_ptr<TimeDisc> _timeDisc;

public:

    void initDomain(
            const Storage<3>& geometry,
            const Topology& topology);

    void initFields(const Storage<3>& initialValues);

    void initSpaceDisc(
            const SpaceDisc::Type discType,
            const double coriolisForce,
            const double frictionForce);

    void initTimeDisc(
            const TimeDisc::Type discType);

    Storage<3> run(const double time);
};
