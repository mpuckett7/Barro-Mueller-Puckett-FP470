/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef FD_BOUNDARIES_2D_HH
#define FD_BOUNDARIES_2D_HH

#include "boundaryPostProcessors2D.h"

#include "utilities/finiteDifference2D.h"
#include "core/util.h"

#include "dynamics/dynamics.h"
#include "dynamics/lbm.h"

namespace olb {

///////////  StraightFdBoundaryProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction,int orientation>
template<CONCEPT(Cell) CELL>
void StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::apply(CELL& cell)
{
  using namespace olb::util::tensorIndices2D;

  T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
  T rho, u[DESCRIPTOR::d];

  auto& dynamics = cell.getDynamics();

  cell.computeRhoU(rho,u);

  interpolateGradients<0>(cell, dx_u);
  interpolateGradients<1>(cell, dy_u);

  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];
  T omega = dynamics.getOmegaOrFallback(std::numeric_limits<T>::signaling_NaN());
  T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
  T pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = dynamics.computeEquilibrium(iPop,rho,u)
               + equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
  }
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
template<int deriveDirection, typename CELL>
void StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(CELL& cell, T velDeriv[DESCRIPTOR::d]) const
{
  fd::DirectedGradients2D<T,DESCRIPTOR,direction,orientation,direction==deriveDirection>
    ::interpolateVector(velDeriv, cell);
}

////////  StraightConvectionBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
StraightConvectionBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), uAv(uAv_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  this->getName() = "StraightConvectionBoundaryProcessor2D";
  saveCell = new T** [(size_t)(x1_-x0_+1)];
  for (int iX=0; iX<=x1_-x0_; ++iX) {
    saveCell[iX] = new T* [(size_t)(y1_-y0_+1)];
    for (int iY=0; iY<=y1_-y0_; ++iY) {
      saveCell[iX][iY] = new T [(size_t)(DESCRIPTOR::q)];
      for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        // default set to -1 in order to avoid wrong results at first call
        saveCell[iX][iY][iPop] = T(-1);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
~StraightConvectionBoundaryProcessor2D()
{
  for (int iX=0; iX<=x1-x0; ++iX) {
    for (int iY=0; iY<=y1-y0; ++iY) {
      delete [] saveCell[iX][iY];
    }
    delete [] saveCell[iX];
  }
  delete [] saveCell;
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY);
        for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
          if (descriptors::c<DESCRIPTOR>(iPop,direction)==-orientation) {
            // using default -1 avoids wrong first call
            if (!util::nearZero(1 + saveCell[iX-newX0][iY-newY0][iPop]) ) {
              cell[iPop] = saveCell[iX-newX0][iY-newY0][iPop];
            }
          }
        }

        T rho0, u0[2];
        T rho1, u1[2];
        T rho2, u2[2];
        if (direction==0) {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX-orientation,iY).computeRhoU(rho1,u1);
          blockLattice.get(iX-orientation*2,iY).computeRhoU(rho2,u2);
        }
        else {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX,iY-orientation).computeRhoU(rho1,u1);
          blockLattice.get(iX,iY-orientation*2).computeRhoU(rho2,u2);
        }

        // rho0 = T(1); rho1 = T(1); rho2 = T(1);

        T uDelta[2];
        T uAverage = rho0*u0[direction];
        if (uAv!=nullptr) {
          uAverage = *uAv;
        }
        uDelta[0]=-uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
        uDelta[1]=-uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);

        for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
          if (descriptors::c<DESCRIPTOR>(iPop,direction) == -orientation) {
            saveCell[iX-newX0][iY-newY0][iPop] = cell[iPop] + descriptors::invCs2<T,DESCRIPTOR>()*descriptors::t<T,DESCRIPTOR>(iPop)*(uDelta[0]*descriptors::c<DESCRIPTOR>(iPop,0)+uDelta[1]*descriptors::c<DESCRIPTOR>(iPop,1));
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  StraightConvectionBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction,int orientation>
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR, direction,orientation>::
StraightConvectionBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), uAv(uAv_)
{ }

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessor2D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, uAv);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, uAv);
}


////////  SlipBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessor2D<T,DESCRIPTOR>::
SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX, int discreteNormalY)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  this->getName() = "SlipBoundaryProcessor2D";
  OLB_PRECONDITION(x0==x1 || y0==y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY < 0) {
      // std::cout << "-----" <<s td::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX*discreteNormalY==0) {
        mult = 2;
      }
      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalX);
      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 ) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          if (reflectionPop[iPop]!=0) {
            //do reflection
            blockLattice.get(iX,iY)[iPop] = blockLattice.get(iX,iY)[reflectionPop[iPop]];
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  SlipBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::
SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new SlipBoundaryProcessor2D<T,DESCRIPTOR>
         ( this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

////////  PartialSlipBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessor2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX, int discreteNormalY)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), tuner(tuner_)
{
  this->getName() = "PartialSlipBoundaryProcessor2D";
  OLB_PRECONDITION(x0==x1 || y0==y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY < 0) {
      //std::cout << "-----" <<std::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX*discreteNormalY==0) {
        mult = 2;
      }
      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalX);
      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 ) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          if (reflectionPop[iPop]!=0) {
            //do reflection
            blockLattice.get(iX,iY)[iPop] = tuner*blockLattice.get(iX,iY)[reflectionPop[iPop]];
          }
        }
        for (int iPop = 1; iPop < DESCRIPTOR::q/2 ; ++iPop) {
          T provv = blockLattice.get(iX,iY)[descriptors::opposite<DESCRIPTOR>(iPop)];
          blockLattice.get(iX,iY)[descriptors::opposite<DESCRIPTOR>(iPop)] +=
            (1.-tuner)*blockLattice.get(iX,iY)[iPop];
          blockLattice.get(iX,iY)[iPop] += (1.-tuner)*provv;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  PartialSlipBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessorGenerator2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), tuner(tuner_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

/////////// OuterVelocityCornerProcessor2D /////////////////////////////////////

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
template <CONCEPT(Cell) CELL>
void OuterVelocityCornerProcessor2D<T,DESCRIPTOR,xNormal,yNormal>::apply(CELL& cell)
{
  using namespace olb::util::tensorIndices2D;

  T rho10 = cell.neighbor({-1*xNormal, -0*yNormal}).computeRho();
  T rho01 = cell.neighbor({-0*xNormal, -1*yNormal}).computeRho();

  T rho20 = cell.neighbor({-2*xNormal, -0*yNormal}).computeRho();
  T rho02 = cell.neighbor({-0*xNormal, -2*yNormal}).computeRho();

  T rho = (T)2/(T)3*(rho01+rho10) - (T)1/(T)6*(rho02+rho20);

  T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
  fd::DirectedGradients2D<T, DESCRIPTOR, 0, xNormal, true>::interpolateVector(dx_u, cell);
  fd::DirectedGradients2D<T, DESCRIPTOR, 1, yNormal, true>::interpolateVector(dy_u, cell);
  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];

  auto& dynamics = cell.getDynamics();
  T omega = dynamics.getOmegaOrFallback(std::numeric_limits<T>::signaling_NaN());

  T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
  T pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  T u[DESCRIPTOR::d];
  cell.computeU(u);

  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = dynamics.computeEquilibrium(iPop,rho,u)
               + equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
  }
}


////////  FreeEnergyWallProcessor2D ////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
template <CONCEPT(Cell) CELL, typename PARAMETERS>
void FreeEnergyWallProcessor2D<T,DESCRIPTOR,NORMAL_X,NORMAL_Y>::apply(CELL& cell, PARAMETERS& vars) any_platform {

  T addend = cell.template getField<descriptors::ADDEND>();

  T rhoAvg = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();
  T rhoTmp = 0.;

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    rhoTmp += cell[iPop];
  }

  T rho = rhoAvg + addend;
  rho -= rhoTmp;

  cell[0] = rho - 1.;
}


////////  FreeEnergyChemPotBoundaryProcessor2D ////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
template <CONCEPT(Cell) CELL>
void FreeEnergyChemPotBoundaryProcessor2DA<T,DESCRIPTOR,NORMAL_X,NORMAL_Y>::apply(CELL& cell) any_platform {

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  T rho0 = cell.computeRho();
  T rho1 = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.template getField<descriptors::CHEM_POTENTIAL>() + (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>());
}

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
template <CONCEPT(Cell) CELL>
void FreeEnergyChemPotBoundaryProcessor2DB<T,DESCRIPTOR,NORMAL_X,NORMAL_Y>::apply(CELL& cell) any_platform {

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());
}


////////  FreeEnergyConvectiveProcessor2D ////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
template <CONCEPT(Cell) CELL>
void FreeEnergyConvectiveProcessor2D<T,DESCRIPTOR,NORMAL_X,NORMAL_Y>::apply(CELL& cell) any_platform {

  T rho, rho0, rho1, u[2];

  rho0 = cell.computeRho();

  cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRhoU(rho1, u);

  T uPerp = 0;

  Vector<T,2> normalVec({NORMAL_X,NORMAL_Y});

  if (normalVec[0] == 0) {
    uPerp = normalVec[1] * u[1];
  } else if (normalVec[1] == 0) {
          uPerp = normalVec[0] * u[0];
  }

  rho = (rho0 + uPerp * rho1) / (1. + uPerp);
  cell.defineRho(rho);
}

}  // namespace olb

#endif
