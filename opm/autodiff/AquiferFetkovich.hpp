/*
Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
Copyright 2017 Statoil ASA.

This file is part of the Open Porous Media project (OPM).

OPM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OPM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_AQUIFETP_HEADER_INCLUDED
#define OPM_AQUIFETP_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/Aquifetp.hpp>
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>
#include <opm/autodiff/BlackoilAquiferModel.hpp>
#include <opm/common/utility/numeric/linearInterpolation.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>

#include <vector>
#include <algorithm>

namespace Opm
{

  template<typename TypeTag>
  class AquiferFetkovich
  {

  public:

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    static const int numEq = BlackoilIndices::numEq;
    typedef double Scalar;

    typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;
    typedef Opm::BlackOilFluidState<Eval, FluidSystem, enableTemperature, enableEnergy, BlackoilIndices::gasEnabled, BlackoilIndices::numPhases> FluidState;

    static const auto waterCompIdx = FluidSystem::waterCompIdx;
    static const auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    AquiferFetkovich( const Aquifetp::AQUFETP_data& aqufetp_data, const Aquancon::AquanconOutput& connection, Simulator& ebosSimulator)
    : ebos_simulator_ (ebosSimulator),
    aqufetp_data_ (aqufetp_data),
    connection_ (connection)
    {}

      void initialSolutionApplied()
      {
        initQuantities(connection_);
      }

      void beginTimeStep()
      {
        ElementContext elemCtx(ebos_simulator_);
        auto elemIt = ebos_simulator_.gridView().template begin<0>();
        const auto& elemEndIt = ebos_simulator_.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
          const auto& elem = *elemIt;

          elemCtx.updatePrimaryStencil(elem);

          int cellIdx = elemCtx.globalSpaceIndex(0, 0);
          int idx = cellToConnectionIdx_[cellIdx];
          if (idx < 0)
          continue;

          elemCtx.updateIntensiveQuantities(0);
          const auto& iq = elemCtx.intensiveQuantities(0, 0);
          pressure_previous_[idx] = Opm::getValue(iq.fluidState().pressure(waterPhaseIdx));
        }
      }

      template <class Context>
      void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx)
      {
        unsigned cellIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        int idx = cellToConnectionIdx_[cellIdx];
        if (idx < 0)
        return;

        // We are dereferencing the value of IntensiveQuantities because cachedIntensiveQuantities return a const pointer to
        // IntensiveQuantities of that particular cell_id
        const IntensiveQuantities intQuants = context.intensiveQuantities(spaceIdx, timeIdx);
        // This is the pressure at td + dt
        updateCellPressure(pressure_current_,idx,intQuants);
        updateCellDensity(idx,intQuants);
        calculateInflowRate(idx, context.simulator());

        rates[BlackoilIndices::conti0EqIdx + FluidSystem::waterCompIdx] +=
        Qai_[idx]/context.dofVolume(spaceIdx, timeIdx);
      }

      void endTimeStep()
      {
        for (const auto& Qai: Qai_) {
          W_flux_ += Qai*ebos_simulator_.timeStepSize();
        }
      }
    private:
      Simulator& ebos_simulator_;

      // Grid variables
      std::vector<size_t> cell_idx_;
      std::vector<Scalar> faceArea_connected_;

      // Quantities at each grid id
      std::vector<Scalar> cell_depth_;
      std::vector<Scalar> pressure_previous_;
      std::vector<Eval> pressure_current_;
      std::vector<Eval> Qai_;
      std::vector<Eval> rhow_;
      std::vector<Scalar> alphai_;
      std::vector<int> cellToConnectionIdx_;

      // Variables constants
      const Aquifetp::AQUFETP_data aqufetp_data_;
      const Aquancon::AquanconOutput connection_;

      Scalar mu_w_    , //water viscosity
      Tc_,
      pa0_     , // initial aquifer pressure
      aquifer_pressure_;

      Eval W_flux_;

      Scalar gravity_() const
      { return ebos_simulator_.problem().gravity()[2]; }

      inline void initQuantities(const Aquancon::AquanconOutput& connection)
      {
        // We reset the cumulative flux at the start of any simulation, so, W_flux = 0
        W_flux_ = 0.;

        // We next get our connections to the aquifer and initialize these quantities using the initialize_connections function
        initializeConnections(connection);

        calculateAquiferCondition();

        pressure_previous_.resize(cell_idx_.size(), 0.);
        pressure_current_.resize(cell_idx_.size(), 0.);
        Qai_.resize(cell_idx_.size(), 0.0);
      }

      inline void updateCellPressure(std::vector<Eval>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
      {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx);
      }

      inline void updateCellPressure(std::vector<Scalar>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
      {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx).value();
      }

      inline void updateCellDensity(const int idx, const IntensiveQuantities& intQuants)
      {
        const auto& fs = intQuants.fluidState();
        rhow_.at(idx) = fs.density(waterPhaseIdx);
      }

      inline Scalar dpai(int idx)
      {
        Scalar dp = aquifer_pressure_ - pressure_previous_.at(idx)+ rhow_.at(idx).value()*gravity_()*(cell_depth_.at(idx) - aqufetp_data_.d0) ;
        std::cout << "Qu =" << idx << " IndexNumber " << std::endl;
        std::cout << "Qu =" <<  pressure_previous_.at(idx) << " Pressure Previous at Index Number" << std::endl;
        //std::cout << "Qu =" <<  pressure_current_.at(idx).value() << " Pressure Current at Index Number" << std::endl;
        //std::cout << "Qu =" <<  rhow_.at(idx).value() << " rho at  Index Number" << std::endl;
        //std::cout << "Qu =" <<  cell_depth_.at(idx) << " Cell depth at  Index Number" << std::endl;
        //std::cout << "Qu =" <<  aqufetp_data_.d0 << " Datum Depth at  Index Number" << std::endl;
        //std::cout << "Qu =" <<  rhow_.at(idx).value()*gravity_*(cell_depth_.at(idx) - aqufetp_data_.d0) << " Constant part at  Index Number" << std::endl;
        return dp;
      }
      // This function implements Eq 5.12 of the EclipseTechnicalDescription
      inline Scalar aquiferPressure()
      {
        Scalar Flux = W_flux_.value();
        Scalar pa_ = pa0_ - Flux / ( aqufetp_data_.C_t * aqufetp_data_.V0 );
        return pa_;
      }
      // This function implements Eq 5.14 of the EclipseTechnicalDescription
      inline void calculateInflowRate(int idx, const Simulator& simulator)
      {
        Tc_ = ( aqufetp_data_.C_t * aqufetp_data_.V0 ) / aqufetp_data_.J ;
        Scalar td_Tc_ = simulator.timeStepSize() / Tc_ ;
        Scalar exp_ = (1 - exp(-td_Tc_)) / td_Tc_;
        Qai_.at(idx) = alphai_.at(idx) * aqufetp_data_.J * dpai(idx) * exp_;
        //std::cout << "Qu =" << idx << " IndexNumber " << std::endl;
        std::cout << "Qu =" << dpai(idx) << " Pressure Difference at Index Number" << std::endl;
        std::cout << "Qu =" << exp_ << " Exponential Constant " << std::endl;
        std::cout << "Qu =" << Qai_.at(idx) << " Inflow at Index Number" << std::endl;
      }

      // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
      inline void initializeConnections(const Aquancon::AquanconOutput& connection)
      {
        const auto& eclState = ebos_simulator_.vanguard().eclState();
        const auto& ugrid = ebos_simulator_.vanguard().grid();
        const auto& grid = eclState.getInputGrid();

        cell_idx_ = connection.global_index;
        auto globalCellIdx = ugrid.globalCell();

        assert( cell_idx_ == connection.global_index);
        assert( (cell_idx_.size() == connection.influx_coeff.size()) );
        assert( (connection.influx_coeff.size() == connection.influx_multiplier.size()) );
        assert( (connection.influx_multiplier.size() == connection.reservoir_face_dir.size()) );

        // We hack the cell depth values for now. We can actually get it from elementcontext pos
        cell_depth_.resize(cell_idx_.size(), aqufetp_data_.d0);
        alphai_.resize(cell_idx_.size(), 1.0);
        faceArea_connected_.resize(cell_idx_.size(),0.0);
        Scalar faceArea;

        auto cell2Faces = Opm::UgGridHelpers::cell2Faces(ugrid);
        auto faceCells  = Opm::AutoDiffGrid::faceCells(ugrid);

        // Translate the C face tag into the enum used by opm-parser's TransMult class
        Opm::FaceDir::DirEnum faceDirection;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas = 0.;
        for (size_t idx = 0; idx < cell_idx_.size(); ++idx)
        {
          auto cellFacesRange = cell2Faces[cell_idx_.at(idx)];

          for(auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter)
          {
            // The index of the face in the compressed grid
            const int faceIdx = *cellFaceIter;

            // the logically-Cartesian direction of the face
            const int faceTag = Opm::UgGridHelpers::faceTag(ugrid, cellFaceIter);

            switch(faceTag)
            {
              case 0: faceDirection = Opm::FaceDir::XMinus;
              break;
              case 1: faceDirection = Opm::FaceDir::XPlus;
              break;
              case 2: faceDirection = Opm::FaceDir::YMinus;
              break;
              case 3: faceDirection = Opm::FaceDir::YPlus;
              break;
              case 4: faceDirection = Opm::FaceDir::ZMinus;
              break;
              case 5: faceDirection = Opm::FaceDir::ZPlus;
              break;
              default: OPM_THROW(Opm::NumericalIssue,"Initialization of Aquifer Carter Tracy problem. Make sure faceTag is correctly defined");
            }

            if (faceDirection == connection.reservoir_face_dir.at(idx))
            {
              // Check now if the face is outside of the reservoir, or if it adjoins an inactive cell
              // Do not make the connection if the product of the two cellIdx > 0. This is because the
              // face is within the reservoir/not connected to boundary. (We still have yet to check for inactive cell adjoining)
              faceArea = (faceCells(faceIdx,0)*faceCells(faceIdx,1) > 0)? 0. : Opm::UgGridHelpers::faceArea(ugrid, faceIdx);
              faceArea_connected_.at(idx) = faceArea;
              denom_face_areas += ( connection.influx_multiplier.at(idx) * faceArea_connected_.at(idx) );
            }
          }
          auto cellCenter = grid.getCellCenter(cell_idx_.at(idx));
          cell_depth_.at(idx) = cellCenter[2];
        }

        for (size_t idx = 0; idx < cell_idx_.size(); ++idx)
        {
          alphai_.at(idx) = ( connection.influx_multiplier.at(idx) * faceArea_connected_.at(idx) )/denom_face_areas;
        }
      }

      inline void calculateAquiferCondition()
      {
        int pvttableIdx = aqufetp_data_.pvttableID - 1;
        rhow_.resize(cell_idx_.size(),0.);
        if (!aqufetp_data_.p0)
        {
          pa0_ = calculateReservoirEquilibrium();
        }
        else
        {
          pa0_ = *(aqufetp_data_.p0);
        }

        aquifer_pressure_ = pa0_ ;
        // use the thermodynamic state of the first active cell as a
        // reference. there might be better ways to do this...
        ElementContext elemCtx(ebos_simulator_);
        auto elemIt = ebos_simulator_.gridView().template begin</*codim=*/0>();
        elemCtx.updatePrimaryStencil(*elemIt);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);

        // Initialize a FluidState object first
        FluidState fs_aquifer;
        // We use the temperature of the first cell connected to the aquifer
        // Here we copy the fluidstate of the first cell, so we do not accidentally mess up the reservoir fs
        fs_aquifer.assign( iq0.fluidState() );
        Eval temperature_aq, pa0_mean;
        temperature_aq = fs_aquifer.temperature(0);
        pa0_mean = pa0_;

        Eval mu_w_aquifer = FluidSystem::waterPvt().viscosity(pvttableIdx, temperature_aq, pa0_mean);

        mu_w_ = mu_w_aquifer.value();
      }

      inline Scalar calculateReservoirEquilibrium()
      {
        // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
        std::vector<Scalar> pw_aquifer;
        Scalar water_pressure_reservoir;

        ElementContext elemCtx(ebos_simulator_);
        const auto& gridView = ebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
          const auto& elem = *elemIt;
          elemCtx.updatePrimaryStencil(elem);

          size_t cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
          int idx = cellToConnectionIdx_[cellIdx];
          if (idx < 0)
          continue;

          elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
          const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
          const auto& fs = iq0.fluidState();

          water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
          rhow_[idx] = fs.density(waterPhaseIdx);
          pw_aquifer.push_back( (water_pressure_reservoir - rhow_[idx].value()*gravity_()*(cell_depth_[idx] - aqufetp_data_.d0))*alphai_[idx] );
        }

        // We take the average of the calculated equilibrium pressures.
        Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.)/pw_aquifer.size();
        return aquifer_pres_avg;
      }
    }; //Class AquiferFetkovich
  } // namespace Opm

  #endif
