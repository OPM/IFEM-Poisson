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
#include "PoissonSolutions.h"
#include "SIMPoisson.h"

#include "SIM2D.h"
#include "Vec3.h"
#include "tinyxml2.h"

#include <cmath>

#include "Catch2Support.h"

namespace {

Vec3 point(double x, double y, double z = 0.0)
{
  Vec3 X;
  X.x = x;
  X.y = y;
  X.z = z;
  return X;
}


class ConstantRealFunc : public RealFunc
{
public:
  explicit ConstantRealFunc(double c) : value(c) {}

protected:
  double evaluate(const Vec3&) const override { return value; }

private:
  double value;
};


class ConstantVecFunc : public VecFunc
{
public:
  explicit ConstantVecFunc(const Vec3& c) : value(c) {}

protected:
  Vec3 evaluate(const Vec3&) const override { return value; }

private:
  Vec3 value;
};


double derivativeX(const VecFunc& flux, const Vec3& X, size_t component,
                   double h = 1.0e-6)
{
  Vec3 xp(X), xm(X);
  xp.x += h;
  xm.x -= h;
  return (flux(xp)[component] - flux(xm)[component])/(2.0*h);
}


double derivativeY(const VecFunc& flux, const Vec3& X, size_t component,
                   double h = 1.0e-6)
{
  Vec3 yp(X), ym(X);
  yp.y += h;
  ym.y -= h;
  return (flux(yp)[component] - flux(ym)[component])/(2.0*h);
}


double derivativeZ(const VecFunc& flux, const Vec3& X, size_t component,
                   double h = 1.0e-6)
{
  Vec3 zp(X), zm(X);
  zp.z += h;
  zm.z -= h;
  return (flux(zp)[component] - flux(zm)[component])/(2.0*h);
}


double divergence2D(const VecFunc& flux, const Vec3& X)
{
  return derivativeX(flux,X,0) + derivativeY(flux,X,1);
}


double divergence3D(const VecFunc& flux, const Vec3& X)
{
  return derivativeX(flux,X,0) + derivativeY(flux,X,1) + derivativeZ(flux,X,2);
}


double derivative(const RealFunc& field, const Vec3& X, char axis,
                  double h = 1.0e-6)
{
  Vec3 xp(X), xm(X);
  switch (axis)
  {
    case 'x': xp.x += h; xm.x -= h; break;
    case 'y': xp.y += h; xm.y -= h; break;
    default:  xp.z += h; xm.z -= h; break;
  }

  return (field(xp) - field(xm))/(2.0*h);
}

}


TEST_CASE("TestSIMPoisson.Parse")
{
  SIMPoisson<SIM2D> sim;
  REQUIRE(sim.read("Square.xinp"));

  const Poisson& poisson = static_cast<const Poisson&>(*sim.getProblem());

  REQUIRE_THAT(poisson.getHeat(Vec3()), WithinRel(M_PI*M_PI*2.0));
  REQUIRE_THAT(poisson.getMaterial(), WithinRel(1.0));
  REQUIRE(poisson.getNoGalerkin() == 1);
}


TEST_CASE("Poisson public API")
{
  Poisson poisson(2);

  REQUIRE(poisson.getHeat(Vec3()) == 0.0);
  REQUIRE(poisson.getFlux(Vec3(),Vec3()) == 0.0);
  REQUIRE_THAT(poisson.getMaterial(), WithinRel(1.0));
  REQUIRE_FALSE(poisson.constrainIntgSol());
  REQUIRE(poisson.getNoGLMs() == 0);

  poisson.setMaterial(2.5);
  REQUIRE_THAT(poisson.getMaterial(), WithinRel(2.5));

  ConstantRealFunc conductivity(4.25);
  poisson.setMaterial(&conductivity);
  REQUIRE_THAT(poisson.getMaterial(point(0.2,0.4)), WithinRel(4.25));

  ConstantVecFunc traction(point(2.0,-1.0));
  poisson.setTraction(&traction);
  REQUIRE_THAT(poisson.getFlux(point(0.1,0.2),point(0.0,1.0)), WithinRel(-1.0));

  ConstantRealFunc flux(6.0);
  poisson.setTraction(&flux);
  REQUIRE_THAT(poisson.getFlux(point(0.1,0.2),point(0.0,1.0)), WithinRel(6.0));

  ConstantRealFunc heat(7.5);
  poisson.setSource(&heat);
  REQUIRE_THAT(poisson.getHeat(point(0.3,0.6)), WithinRel(7.5));

  poisson.constrainIntgSol(true);
  REQUIRE(poisson.constrainIntgSol());
  REQUIRE(poisson.getNoGLMs() == 1);

  REQUIRE(poisson.getNoFields(1) == 1);
  REQUIRE(poisson.getNoFields(2) == 4);
  REQUIRE(poisson.getField1Name(0,nullptr) == "u");
  REQUIRE(poisson.getField1Name(0,"prefix") == "prefix u");
  REQUIRE(poisson.getField2Name(0,nullptr) == "q_x");
  REQUIRE(poisson.getField2Name(1,"flux") == "flux q_y");
  REQUIRE(poisson.getField2Name(2,nullptr) == "(q_ex)_x");
  REQUIRE(poisson.getField2Name(3,nullptr) == "(q_ex)_y");
  REQUIRE(poisson.getField2Name(4,nullptr).empty());
  REQUIRE_FALSE(poisson.suppressOutput(1,ASM::SECONDARY));
  REQUIRE(poisson.suppressOutput(2,ASM::SECONDARY));
}


TEST_CASE("Poisson extraction and galerkin helpers")
{
  Poisson poisson(2);
  ConstantRealFunc extract1(1.0);
  ConstantRealFunc extract2(2.0);
  ConstantVecFunc vectorField(point(1.0,2.0));

  poisson.addExtrFunction(&extract1);
  poisson.addExtrFunction(&extract2);
  REQUIRE(poisson.numExtrFunction() == 2);

  poisson.setMode(SIM::STATIC);
  REQUIRE(poisson.getExtractionField(1) != nullptr);
  REQUIRE(poisson.getExtractionField(3) == nullptr);

  poisson.addExtrFunction(&vectorField);
  REQUIRE(poisson.numExtrFunction() == 0);
  poisson.setMode(SIM::STATIC);
  REQUIRE(poisson.getExtractionField(1) == nullptr);

  tinyxml2::XMLDocument doc;
  REQUIRE(doc.Parse("<root><galerkin>x;y;0</galerkin><residual/></root>") == tinyxml2::XML_SUCCESS);
  const tinyxml2::XMLElement* root = doc.FirstChildElement("root");
  REQUIRE(root != nullptr);
  REQUIRE(poisson.parse(root->FirstChildElement("galerkin")));
  REQUIRE(poisson.getNoGalerkin() == 1);
  REQUIRE(poisson.parse(root));
  poisson.clearGalerkinProjections();
  REQUIRE(poisson.getNoGalerkin() == 0);

  tinyxml2::XMLDocument invalid;
  REQUIRE(invalid.Parse("<foo/>") == tinyxml2::XML_SUCCESS);
  REQUIRE_FALSE(poisson.parse(invalid.FirstChildElement("foo")));
}


TEST_CASE("Poisson analytic solutions match expected values")
{
  const Vec3 X2 = point(0.25,1.0);
  const Vec3 Xsin = point(0.125,0.125);
  const Vec3 Xline = point(0.5,0.0);

  const Square2D squareFlux;
  const Square2DHeat squareHeat;
  const Vec3 q2 = squareFlux(X2);
  REQUIRE_THAT(q2.x, WithinRel(-M_PI/std::sqrt(2.0)));
  REQUIRE_THAT(q2.y, WithinRel(-1.0/std::sqrt(2.0)));
  REQUIRE_THAT(squareHeat(X2), WithinRel(-M_PI*M_PI/std::sqrt(2.0)));

  const SquareSinus sinusFlux;
  const SquareSinusSource sinusHeat;
  const Vec3 qsin = sinusFlux(Xsin);
  REQUIRE_THAT(qsin.x, WithinRel(M_PI));
  REQUIRE_THAT(qsin.y, WithinRel(M_PI));
  REQUIRE_THAT(sinusHeat(Xsin), WithinRel(-4.0*M_PI*M_PI));

  const PoissonCube cubeFlux;
  const PoissonCubeSource cubeHeat;
  const Vec3 qcube = cubeFlux(point(0.0,0.5,0.5));
  REQUIRE_THAT(qcube.x, WithinRel(-M_PI));
  REQUIRE(std::abs(qcube.y) < 1.0e-12);
  REQUIRE(std::abs(qcube.z) < 1.0e-12);
  REQUIRE(std::abs(cubeHeat(point(0.0,0.5,0.5))) < 1.0e-12);

  const PoissonLine lineFlux(2.0);
  const PoissonLineSource lineHeat(2.0);
  REQUIRE_THAT(lineFlux(Xline).x, WithinRel(-0.5*M_PI/std::sqrt(2.0)));
  REQUIRE_THAT(lineHeat(Xline), WithinRel(0.25*M_PI*M_PI/std::sqrt(2.0)));

  const LshapePoisson lshape;
  const Vec3 qL = lshape(Vec3());
  REQUIRE(std::isfinite(qL.x));
  REQUIRE(std::isfinite(qL.y));
}


TEST_CASE("Poisson sources stay consistent with analytic fields")
{
  const Square2D squareFlux;
  const Square2DHeat squareHeat;
  const Vec3 X2 = point(0.31,0.42);
  REQUIRE(std::abs(divergence2D(squareFlux,X2) - squareHeat(X2)) < 1.0e-4);

  const SquareSinus sinusFlux;
  const SquareSinusSource sinusHeat;
  const Vec3 Xsin = point(0.19,0.37);
  REQUIRE(std::abs(divergence2D(sinusFlux,Xsin) - sinusHeat(Xsin)) < 1.0e-3);

  const PoissonCube cubeFlux;
  const PoissonCubeSource cubeHeat;
  const Vec3 X3 = point(0.21,0.37,0.41);
  REQUIRE(std::abs(divergence3D(cubeFlux,X3) - cubeHeat(X3)) < 1.0e-4);

  const PoissonLine lineFlux(2.0);
  const PoissonLineSource lineHeat(2.0);
  const Vec3 XL = point(0.37,0.0);
  REQUIRE(std::abs(derivativeX(lineFlux,XL,0) - lineHeat(XL)) < 1.0e-5);

  const PoissonInteriorLayer interiorFlux(12.0);
  const PoissonInteriorLayerSol interiorSol(12.0);
  const Vec3 Xi = point(0.8,0.15);
  const Vec3 qi = interiorFlux(Xi);
  REQUIRE(std::abs(qi.x + derivative(interiorSol,Xi,'x')) < 1.0e-4);
  REQUIRE(std::abs(qi.y + derivative(interiorSol,Xi,'y')) < 1.0e-4);

  const PoissonWaterfall waterfallFlux(0.01);
  const PoissonWaterfallSol waterfallSol(0.01);
  const Vec3 Xw = point(0.2,0.15,0.4);
  const Vec3 qw = waterfallFlux(Xw);
  REQUIRE(std::abs(qw.x + derivative(waterfallSol,Xw,'x')) < 1.0e-4);
  REQUIRE(std::abs(qw.y + derivative(waterfallSol,Xw,'y')) < 1.0e-4);
  REQUIRE(std::abs(qw.z + derivative(waterfallSol,Xw,'z')) < 1.0e-4);
}
