//==============================================================================
//!
//! \file TestSIMPoisson.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "Poisson.h"
#include "SIMPoisson.h"

#include "SIM2D.h"
#include "Vec3.h"

#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinRel;


TEST_CASE("TestSIMPoisson.Parse")
{
  SIMPoisson<SIM2D> sim;
  REQUIRE(sim.read("Square.xinp"));

  const Poisson& poisson = static_cast<const Poisson&>(*sim.getProblem());

  REQUIRE_THAT(poisson.getHeat(Vec3()), WithinRel(M_PI*M_PI*2.0));
  REQUIRE_THAT(poisson.getMaterial(), WithinRel(1.0));
  REQUIRE(poisson.getNoGalerkin() == 1);
}
