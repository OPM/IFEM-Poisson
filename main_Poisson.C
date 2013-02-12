// $Id$
//==============================================================================
//!
//! \file main_Poisson.C
//!
//! \date 20 May 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Main program for the isogeometric solver for the Poisson equation.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMPoisson.h"
#include "AdaptiveSIM.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <cstdio>


/*!
  \brief Main program for the NURBS-based isogeometric Poisson equation solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -1D : Use one-parametric simulation driver
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
  \arg -DGL2 : Estimate error using discrete global L2 projection
  \arg -CGL2 : Estimate error using continuous global L2 projection
  \arg -SCR : Estimate error using Superconvergent recovery at Greville points
  \arg -VDSA: Estimate error using Variational Diminishing Spline Approximations
  \arg -LSQ : Estimate error using through Least Square projections
  \arg -QUASI : Estimate error using Quasi-interpolation projections
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  std::vector<int> ignoredPatches;
  size_t adaptor = 0;
  int  i, iop = 0;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  char ndim = 3;
  char* infile = 0;

  int myPid = IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strncmp(argv[i],"-adap",5))
    {
      iop = 10;
      if (strlen(argv[i]) > 5)
        adaptor = atoi(argv[i]+5);
    }
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-free] [-lag|-spec|-LR] [-1D|-2D] [-adap[<i>]] [-nGauss <n>]"
	      <<"\n       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
	      <<"       [-DGL2] [-CGL2] [-SCR] [-VDLSA] [-LSQ] [-QUASI]\n"
	      <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup]"
	      <<" [-checkRHS] [-check]\n";
    return 0;
  }

  if (myPid == 0)
  {
    const SIMoptions& opts = IFEM::getOptions();
    std::cout <<"\n >>> IFEM Poisson equation solver <<<"
	      <<"\n ====================================\n"
	      <<"\n Executing command:\n";
    for (i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< opts.solver
	      <<"\nNumber of Gauss points: "<< opts.nGauss[0];
    if (opts.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (opts.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< opts.nViz[0];
      if (ndim > 1) std::cout <<" "<< opts.nViz[1];
      if (ndim > 2) std::cout <<" "<< opts.nViz[2];
    }

    if (opts.eig > 0)
      std::cout <<"\nEigenproblem solver: "<< opts.eig
		<<"\nNumber of eigenvalues: "<< opts.nev
		<<"\nNumber of Arnoldi vectors: "<< opts.ncv
		<<"\nShift value: "<< opts.shift;
    if (opts.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (opts.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (opts.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    if (SIMbase::ignoreDirichlet)
      std::cout <<"\nSpecified boundary conditions are ignored";
    if (fixDup)
      std::cout <<"\nCo-located nodes will be merged";
    if (checkRHS)
      std::cout <<"\nCheck that each patch has a right-hand coordinate system";
    if (!ignoredPatches.empty())
    {
      std::cout <<"\nIgnored patches:";
      for (size_t i = 0; i < ignoredPatches.size(); i++)
	std::cout <<" "<< ignoredPatches[i];
    }
    std::cout << std::endl;
  }
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the simulation model
  SIMbase* model;
  if (ndim == 1)
    model = new SIMPoisson1D();
  else if (ndim == 2)
    model = new SIMPoisson2D();
  else
    model = new SIMPoisson3D(checkRHS);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = 0;
  if (iop == 10)
  {
    theSim = aSim = new AdaptiveSIM(model);
    IFEM::getOptions().discretization = ASM::LRSpline;
  }

  // Read in model definitions
  if (!theSim->read(infile))
    return 1;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (model->opt.eig != 4 && model->opt.eig != 6)
    SIMbase::ignoreDirichlet = false;

  // Load vector visualization is not available when using additional viz-points
  for (i = 0; i < 3; i++)
    if (i >= ndim)
      model->opt.nViz[i] = 1;
    else if (model->opt.nViz[i] > 2)
      vizRHS = false;

  if (myPid == 0)
  {
    std::cout <<"\n\nEquation solver: "<< model->opt.solver
	      <<"\nNumber of Gauss points: "<< model->opt.nGauss[0]
	      <<" "<< model->opt.nGauss[1];
    if (model->opt.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (model->opt.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< model->opt.nViz[0];
      if (ndim > 1) std::cout <<" "<< model->opt.nViz[1];
      if (ndim > 2) std::cout <<" "<< model->opt.nViz[2];
    }
    if (model->opt.eig > 0)
      std::cout <<"\nEigenproblem solver: "<< model->opt.eig
		<<"\nNumber of eigenvalues: "<< model->opt.nev
		<<"\nNumber of Arnoldi vectors: "<< model->opt.ncv
		<<"\nShift value: "<< model->opt.shift;
    if (model->opt.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (model->opt.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (model->opt.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    std::cout << std::endl;
  }

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess(ignoredPatches,fixDup))
    return 1;

  if (model->opt.discretization < ASM::Spline && !model->opt.hdf5.empty())
  {
    std::cout <<"\n ** HDF5 output is available for spline discretization only."
	      <<" Deactivating...\n"<< std::endl;
    model->opt.hdf5.clear();
  }

  SIMoptions::ProjectionMap& pOpt = model->opt.project;
  SIMoptions::ProjectionMap::const_iterator pit;

  // Set default projection method (tensor splines only)
  bool staticSol = iop + model->opt.eig == 0 || iop == 10;
  if (model->opt.discretization < ASM::Spline || !staticSol)
    pOpt.clear(); // No projection if Lagrange/Spectral or no static solution
  else if (model->opt.discretization == ASM::Spline && pOpt.empty())
    pOpt[SIMoptions::GLOBAL] = "Greville point projection";

  const char* prefix[pOpt.size()];
  if (model->opt.format >= 0 || model->opt.dumpHDF5(infile))
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++)
      prefix[i] = pit->second.c_str();

  model->setQuadratureRule(model->opt.nGauss[0],true);

  Matrix eNorm, ssol;
  Vector sol, load;
  Vectors projs(pOpt.size()), gNorm;
  std::vector<Mode> modes;
  int iStep = 1, nBlock = 0;
  bool iterate = true;

  if (aSim)
    aSim->setupProjections();

  DataExporter* exporter = NULL;
  if (model->opt.dumpHDF5(infile) && staticSol)
  {
    if (myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< model->opt.hdf5
                <<".hdf5"<< std::endl;

    // Include secondary results only if no projection has been requested.
    // The secondary results will be projected anyway, but without the
    // nodal averaging across patch boundaries in case of multiple patches.
    int results = DataExporter::PRIMARY | DataExporter::NORMS;
    if (pOpt.empty()) results |= DataExporter::SECONDARY;

    exporter = new DataExporter(true);
    exporter->registerField("u","heat",DataExporter::SIM,results);
    exporter->setFieldValue("u",model, aSim ? &aSim->getSolution() : &sol);
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++) {
      exporter->registerField(prefix[i], "projected", DataExporter::SIM,
                              DataExporter::SECONDARY, prefix[i]);
      exporter->setFieldValue(prefix[i], model,
                              aSim ? &aSim->getProjection(i) : &projs[i]);
    }
    exporter->registerWriter(new HDF5Writer(model->opt.hdf5));
    exporter->registerWriter(new XMLWriter(model->opt.hdf5));
    exporter->setNormPrefixes(prefix);
  }

  switch (iop+model->opt.eig) {
  case 0:
    model->setMode(SIM::STATIC);
    model->initSystem(model->opt.solver,1,1);
    model->setAssociatedRHS(0,0);
    if (!model->assembleSystem())
      return 2;
    else if (vizRHS)
      model->extractLoadVec(load);

    // Solve the linear system of equations
    if (!model->solveSystem(sol,1))
      return 3;

    // Project the FE stresses onto the splines basis
    model->setMode(SIM::RECOVERY);
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++)
      if (!model->project(ssol,sol,pit->first))
	return 4;
      else
	projs[i] = ssol;

    if (myPid == 0 && !pOpt.empty())
      std::cout << std::endl;

    // Evaluate solution norms
    model->setQuadratureRule(model->opt.nGauss[1]);
    if (!model->solutionNorms(Vectors(1,sol),projs,eNorm,gNorm))
      return 4;

    if (myPid == 0)
    {
      model->printNorms(gNorm,std::cout);
      size_t j = 2;
      for (pit = pOpt.begin(); pit != pOpt.end() && j < gNorm.size(); pit++,j++)
      {
	std::cout <<"\n\n>>> Error estimates based on "<< pit->second <<" <<<";
	std::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm[j](1);
	std::cout <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm[j](2);
	std::cout <<"\n- relative error (% of |u^r|) : "
		  << gNorm[j](2)/gNorm[j](1)*100.0;
	if (model->haveAnaSol() && j <= gNorm.size())
	{
	  std::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm[j](3)
		    <<"\n- relative error (% of |u|)   : "
		    << gNorm[j](3)/gNorm[0](3)*100.0;
	  std::cout <<"\nEffectivity index             : "
		    << gNorm[j](2)/gNorm[0](4);
        }
      }
    }
    break;

  case 10:
    // Adaptive simulation
    if (!aSim->initAdaptor(adaptor,2))
      break;

    if (exporter)
      exporter->setNormPrefixes(aSim->getNormPrefixes());

    while (iterate) {
      char iterationTag[256];
      sprintf(iterationTag, "Adaptive step #%03d", iStep);
      utl::profiler->start(iterationTag);
      if (!aSim->solveStep(infile,iStep))
        return 5;
      else if (!aSim->writeGlv(infile,iStep,nBlock,2))
        return 6;
      else if (exporter)
        exporter->dumpTimeLevel(NULL,true);
      utl::profiler->stop(iterationTag);
      sprintf(iterationTag, "Refinement step #%03d", iStep);
      utl::profiler->start(iterationTag);
      iterate = aSim->adaptMesh(++iStep);
      utl::profiler->stop(iterationTag);
    }

  case 100:
    break; // Model check

  default:
    // Free vibration: Assemble coefficient matrix [K]
    model->setMode(SIM::VIBRATION);
    model->initSystem(model->opt.solver,1,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes))
      return 6;
  }

  utl::profiler->start("Postprocessing");

  if (iop != 10 && model->opt.format >= 0)
  {
    // Write VTF-file with model geometry
    if (!model->writeGlvG(nBlock,infile))
      return 7;

    // Write boundary tractions, if any
    if (!model->writeGlvT(iStep,nBlock))
      return 8;

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      return 8;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",iStep,nBlock))
      return 9;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(sol,iStep,nBlock))
      return 10;

    // Write projected solution fields to VTF-file
    size_t i = 0;
    int iBlk = 100;
    for (pit = pOpt.begin(); pit != pOpt.end(); pit++, i++, iBlk += 10)
      if (!model->writeGlvP(projs[i],iStep,nBlock,iBlk,pit->second.c_str()))
        return 11;
      else
	prefix[i] = pit->second.c_str();

    // Write eigenmodes
    bool isFreq = model->opt.eig==3 || model->opt.eig==4 || model->opt.eig==6;
    for (i = 0; i < modes.size(); i++)
      if (!model->writeGlvM(modes[i],isFreq,nBlock))
	return 11;

    // Write element norms
    if (!model->writeGlvN(eNorm,iStep,nBlock,prefix))
      return 12;

    model->writeGlvStep(1,0.0,1);
  }
  model->closeGlv();
  if (exporter && iop != 10)
    exporter->dumpTimeLevel();

  utl::profiler->stop("Postprocessing");
  delete theSim;
  delete exporter;
  return 0;
}