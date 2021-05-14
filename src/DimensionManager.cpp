#include <DimensionManager.h>

constexpr double g = 1.0;

DimensionManager::DimensionManager(const Parser& Parser)
: m_h0(Parser.Get("Common", "vertical_measure"))
, m_l0(Parser.Get("Common", "horizontal_measure"))
, m_c0(sqrt(g*m_h0))
{}
