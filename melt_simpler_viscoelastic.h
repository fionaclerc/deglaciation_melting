/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_melt_simpler_viscoelastic_h
#define _aspect_material_model_melt_simpler_viscoelastic_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that takes the implentation of melt production
     * from the melt simple script (dry peridotite solidus of Katz et al. (2003))
     * expansion of melt fraction partial derivatives from Dannberg et al. (2016),
     * and employs the viscoelastic rheology by tracking old stresses as
     * compositional fields.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltSimplerViscoelastic : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:
        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;

        double reference_darcy_coefficient () const override;

        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * @}
         */

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        double reference_rho_s;
        double reference_rho_f;
        double reference_T;
        double eta_0;
        double mu_0;
        double xi_0;
        double K_phi;
        double q_phi;
        double m_phi;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;
        double max_viscosity_scaling;
        double min_viscosity_scaling;
        double max_bulk_viscosity_scaling;
        double min_bulk_viscosity_scaling;
        double thermal_expansivity;
        bool constant_density;
        double reference_specific_heat;
        double thermal_conductivity;
        double reference_permeability;
        double alpha_phi;
        double extraction_depth;
        double compressibility;
        double melt_compressibility;
        double melt_bulk_modulus_derivative;
        double depletion_density_change;
        double depletion_solidus_change;
        double depletion_viscosity_change;
        bool model_is_compressible;
        bool fractional_melting;
        double freezing_rate;
        double melting_time_scale;
        double eruptive_time_scale;
        bool melt_extraction_dikes;
        double reference_tensile_strength;
        double thermal_tensile_strength_exponent;
        double max_tensile_strength_scaling;
        double min_tensile_strength_scaling;
        double eruption_rate;
        double diking_timescale;
        double use_effective_stress;
        bool update_porosity_after_eruption;
        double limiting_time_scale;
        double limiting_time_scale_time_step_factor;
        bool eruption_and_melting;

        bool weaken_ridge;
        double weaken_length;
        double weaken_gradient;

        double cohesion;
        double phi;
        double reference_strain_rate;
        double minimum_strain_rate;

        double activation_energy;
        double activation_volume;
        double viscosity_scale;
        double min_viscosity;
        double distance_melting_suppressed;
        /**
         *
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        // entropy change upon melting
        double peridotite_melting_entropy_change;

        /**
         * Percentage of material that is molten for a given @p temperature and
         * @p pressure (assuming equilibrium conditions). Melting model after Katz,
         * 2003, for dry peridotite.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

        /**
         * Compute the change in entropy due to melting for a given @p temperature
         * and @p pressure, and under the assumption that a fraction
         * @p maximum_melt_fraction of the material has already been molten
         * previously. The entropy change is computed with respect to temperature
         * or pressure, depending on @p dependence.
         * This is needed to calculate the latent heat of melt.
         */
        virtual
        double
        entropy_change (const double temperature,
                        const double pressure,
                        const double maximum_melt_fraction,
                        const NonlinearDependence::Dependence dependence) const;

        Rheology::Elasticity<dim> elastic_rheology;

        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;
        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;


    };

  }
}

#endif
