#include "DimensionManager.h"

enum class scales { height, length, velocity, source, time };

constexpr double g = 9.8;

DimensionManager::DimensionManager(Parser const& Parser)
: Parser_(Parser)
, h0_(Parser.get<double>("Common", "vertical_measure"))
, l0_(Parser.get<double>("Common", "horizontal_measure"))
, c0_(sqrt(g*h0_))
{}