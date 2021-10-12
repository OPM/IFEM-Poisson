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
#include <string>

#include "gtest/gtest.h"


TEST(TestSIMPoisson, Parse)
{
  SIMPoisson<SIM2D> sim;
  EXPECT_TRUE(sim.read("Square.xinp"));

  const Poisson& poisson = static_cast<const Poisson&>(*sim.getProblem());

  ASSERT_FLOAT_EQ(poisson.getHeat(Vec3()), M_PI*M_PI*2.0);
  ASSERT_FLOAT_EQ(poisson.getMaterial(), 1.0);
  EXPECT_EQ(poisson.getNoGalerkin(), 1U);
}
