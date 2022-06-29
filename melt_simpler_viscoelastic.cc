/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#include <aspect/simulator_access.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/material_model/melt_simpler_viscoelastic.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <numeric>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>
#include <aspect/melt.h>
#include <aspect/geometry_model/box.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MeltSimplerViscoelastic<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }

    template <int dim>
    double
    MeltSimplerViscoelastic<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    MeltSimplerViscoelastic<dim>::
    is_compressible () const
    {
      return model_is_compressible;
    }

    template <int dim>
    double
    MeltSimplerViscoelastic<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;
    }


    template <int dim>
    double
    MeltSimplerViscoelastic<dim>::
    entropy_change (const double temperature,
                    const double pressure,
                    const double maximum_melt_fraction,
                    const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;

      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

      if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
          // melt fraction when clinopyroxene is still present
          double melt_fraction_derivative_temperature
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              / (T_lherz_liquidus - T_solidus);

          double melt_fraction_derivative_pressure
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                 + dT_lherz_liquidus_dp * (T_solidus - temperature))
              / pow(T_lherz_liquidus - T_solidus,2);

          // melt fraction after melting of all clinopyroxene
          const double R_cpx = r1 + r2 * std::max(0.0, pressure);
          const double F_max = M_cpx / R_cpx;

          if (melt_fraction(temperature, pressure) > F_max)
            {
              const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
              const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
              const double dT_max_dp = dT_solidus_dp
                                       + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                       + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

              melt_fraction_derivative_temperature
                = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  / (T_liquidus - T_max);

              melt_fraction_derivative_pressure
                = dF_max_dp
                  - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                  + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
            }

          double melt_fraction_derivative = 0;
          if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("not implemented"));

          if (melt_fraction(temperature, pressure) >= maximum_melt_fraction)
            entropy_gradient = melt_fraction_derivative * peridotite_melting_entropy_change;
        }
      return entropy_gradient;
    }


    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.temperature.size(); ++q)
        melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                this->get_adiabatic_conditions().pressure(in.position[q]));
    }


    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::initialize ()
    {
      if (this->include_melt_transport())
        {

          AssertThrow(this->get_parameters().use_operator_splitting,
                      ExcMessage("The material model ``Melt simpler viscoelastic'' can only be used with operator splitting!"));
          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt simpler viscoelastic only works if there is a "
                                 "compositional field called peridotite."));
          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Melt simpler viscoelastic with melt transport only "
                                 "works if there is a compositional field called porosity."));

        }
    }


    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();
      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);
      // feed this the old compaction pressures
      AdditionalMaterialOutputsStokesRHS<dim> *force = out.template get_additional_output<AdditionalMaterialOutputsStokesRHS<dim> >();


      const GeometryModel::Box<dim> &box_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->get_geometry_model());
      ComponentMask composition_mask(this->n_compositional_fields(), true);
      // assign compositional fields associated with viscoelastic stress a value of 0
      // assume these fields are listed first
      for (unsigned int i=0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        composition_mask.set(i,false);
      std::vector<double> average_elastic_shear_moduli (in.temperature.size());
      std::vector<double> elastic_shear_moduli(elastic_rheology.get_elastic_shear_moduli());
      std::vector<double> viscoelastic_viscosities (in.temperature.size());

      // declare these for compaction viscosity
      const double dte(elastic_rheology.elastic_timestep());
      const double dt = this->get_timestep();
      const double n_dislocation = 3.5;

      // Parameters to calculate melt prod
      std::vector<double> melt_production_rates(in.position.size());
      std::vector<double> old_pressures(in.position.size());
      std::vector<double> old_temperatures(in.position.size());
      std::vector<Tensor<1,dim> > temperature_gradients(in.position.size());
      std::vector<Tensor<1,dim> > pressure_gradients(in.position.size());
      std::vector<double> melt_fraction_temperature_derivative(in.position.size());
      std::vector<double> melt_fraction_pressure_derivative(in.position.size());
      std::vector<double> pressure_material_derivative(in.position.size());
      std::vector<double> temperature_material_derivative(in.position.size());

      if (in.current_cell.state() == IteratorState::valid)
        {
           std::vector<Point<dim> > quadrature_positions(in.position.size());
           for (unsigned int i=0; i < in.position.size(); ++i)
             quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);
           FEValues<dim> fe_values (this->get_mapping(),
                                    this->get_fe(),
                                    Quadrature<dim>(quadrature_positions),
                                    update_values | update_gradients);

           fe_values.reinit (in.current_cell);
 
           const FEValuesExtractors::Scalar extractor_temperature = this->introspection().variable("temperature").extractor_scalar();
           fe_values[extractor_temperature].get_function_gradients (this->get_solution(),
                                                                    temperature_gradients);
           fe_values[extractor_temperature].get_function_values (this->get_old_solution(),
                                                                 old_temperatures);
           const FEValuesExtractors::Scalar extractor_pressure = this->introspection().variable("pressure").extractor_scalar();
           fe_values[extractor_pressure].get_function_values (this->get_old_solution(),
                                                                old_pressures);
           fe_values[extractor_pressure].get_function_gradients (this->get_solution(),
                                                                    pressure_gradients);
         }

      for (unsigned int i=0; i<in.position.size(); ++i)
        {



          // calculate density first, we need it for the reaction term
          // first, calculate temperature dependence of density
          double temperature_dependence = 1.0;
          if (!constant_density)
            {
              if (this->include_adiabatic_heating ())
                {
                // temperature dependence is 1 - alpha * (T - T(adiabatic))
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
                }
            else
              temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
            }
          // calculate composition dependence of density

          const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                   :
                                   0.0;

          out.densities[i] = (reference_rho_s + delta_rho)
                             * temperature_dependence * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));

          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, composition_mask);

          equation_of_state.evaluate(in, i, eos_outputs);



          double visc_temperature_dependence = 1.0;
          double bulk_visc_temperature_dependence = 1.0;
          double tensile_strength_temperature_dependence = 1.0;
          double weaken_factor  = 1.0;
          const double visc_depletion_dependence = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   std::exp(depletion_viscosity_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")])
                                   :
                                   1.0;
          // calculate viscosity temperature dependence using explicit Arhenius relation.
          // note not modified in case ad. heating is turned off.
          visc_temperature_dependence = std::exp(std::max(activation_energy  + this->get_adiabatic_conditions().pressure(in.position[i]) * activation_volume, 0.0)/(n_dislocation * in.temperature[i] * constants::gas_constant));
          visc_temperature_dependence *= std::pow(viscosity_scale , -1/n_dislocation)/ std::pow(reference_strain_rate, (n_dislocation-1)/n_dislocation);// here v_scale is A d^-p fH2O^r
          visc_temperature_dependence = 1/(1/visc_temperature_dependence + 1/(max_viscosity_scaling * eta_0))/eta_0; // divide since will be multiplied back in.

          if (this->include_adiabatic_heating ())
            {
              const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
              bulk_visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),max_viscosity_scaling),min_viscosity_scaling);
              tensile_strength_temperature_dependence = std::max(std::min(std::exp(-thermal_tensile_strength_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),max_tensile_strength_scaling),min_tensile_strength_scaling);
              if (weaken_ridge)
                {
                  const double weaken_x = max_viscosity_scaling * (1.0+std::tanh((in.position[i][0] - weaken_length)/weaken_gradient))/2.0;
                  weaken_factor = 1.0 * weaken_x / std::min( std::max( std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])), weaken_x) , max_viscosity_scaling);
             //     visc_temperature_dependence *= weaken_factor;
                }
            }
          else
            {
              const double delta_temp = in.temperature[i]-reference_T;
              const double T_dependence = (thermal_viscosity_exponent == 0.0
                                           ?
                                           0.0
                                           :
                                           thermal_viscosity_exponent*delta_temp/reference_T);
              bulk_visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),max_bulk_viscosity_scaling),min_bulk_viscosity_scaling);
              tensile_strength_temperature_dependence = std::max(std::min(std::exp(-T_dependence),max_tensile_strength_scaling),min_tensile_strength_scaling);
              if (weaken_ridge)
                {
                  const double weaken_x = max_viscosity_scaling * (1.0+std::tanh((in.position[i][0] - weaken_length)/weaken_gradient))/2.0;
                  weaken_factor = 1.0 * weaken_x / std::min( std::max( std::exp(-T_dependence), weaken_x) , max_viscosity_scaling);
               //   visc_temperature_dependence *= weaken_factor;
                }
            }

          const double tensile_strength = reference_tensile_strength * max_tensile_strength_scaling / tensile_strength_temperature_dependence;

          double eq_melt_fraction = 0.0;
          double visc_porosity_dependence = 1.0;

          // Average viscosity and shear modulus
          average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions, elastic_shear_moduli, viscosity_averaging);
          double melting_on = 1.0;

          // for the subduction case - no melting in subducting mantle beneath oceanic plate.
          if (this->introspection().compositional_name_exists("subducting_mantle") && in.composition[i][this->introspection().compositional_index_for_name("subducting_mantle")] > 0.0)
            melting_on = 0.0;

          // for the plume case - no melting at the edges
          if (in.position[i][0] < distance_melting_suppressed || in.position[i][0] > box_geometry_model.get_extents()[0] - distance_melting_suppressed)
            melting_on = 0.0;

          if (dim==3 && ( in.position[i][1] > box_geometry_model.get_extents()[1] - distance_melting_suppressed))
            melting_on = 0.0;
            {

              out.entropy_derivative_pressure[i]    = melting_on * entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0, NonlinearDependence::pressure);
              out.entropy_derivative_temperature[i] = melting_on * entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0, NonlinearDependence::temperature);

             
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double old_porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              if (fractional_melting && melting_on > 0.0)
              {
                // solidus is lowered by previous melting events (fractional melting)
                const double solidus_change = (- old_porosity) * depletion_solidus_change;
                eq_melt_fraction = melt_fraction(in.temperature[i] - solidus_change, this->get_adiabatic_conditions().pressure(in.position[i]));
              }

              // no melting/freezing is used in the model --> set all reactions to zero

              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  out.reaction_terms[i][c] = 0.0;

                  if (reaction_rate_out != nullptr)
                    if (this->introspection().compositional_name_exists("porosity") && c == porosity_idx && this->get_timestep_number() > 0)
                      reaction_rate_out->reaction_rates[i][c] = (-old_porosity + eq_melt_fraction)/dt; 
                    else
                      reaction_rate_out->reaction_rates[i][c] = 0.0;

                }
              visc_porosity_dependence = exp(- alpha_phi * eq_melt_fraction);// theoretically should be the same as above, no negative values.

            }





          viscoelastic_viscosities[i] = elastic_rheology.calculate_viscoelastic_viscosity(eta_0 * visc_depletion_dependence * visc_temperature_dependence
                                                                                 * visc_porosity_dependence,
                                                                                 average_elastic_shear_moduli[i]);

          if (weaken_ridge)
            viscoelastic_viscosities[i] *= weaken_factor;

          // Fill the material properties that are part of the elastic additional outputs
          if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
            { 
              elastic_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
            }

          double viscosities_vep = viscoelastic_viscosities[i];

          if (in.strain_rate.size())
            {
              SymmetricTensor<2,dim> stress_old;
              for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = composition[j];

              const SymmetricTensor<2,dim> edot = 2. * (deviator(in.strain_rate[i])) + stress_old / (average_elastic_shear_moduli[i] * dte);
              const double edot_ii = ( (this->get_timestep_number() == 0 && in.strain_rate[i].norm() <= std::numeric_limits<double>::min() )
                                           ?
                                           reference_strain_rate
                                           :
                                           std::max(std::sqrt(std::fabs(second_invariant(edot))), minimum_strain_rate) );

              const double stress_ve = viscoelastic_viscosities[i] * edot_ii;

              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double old_porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              const double effective_pressure = in.pressure[i];
              const double transition_pressure = (cohesion * std::cos(phi) - tensile_strength) / (1.0 - std::sin(phi));
              double stress_yield;
               
              stress_yield = ( (dim==3)
                               ?
                               ( 6.0 * cohesion * std::cos(phi) + 2.0 * effective_pressure * std::sin(phi) )
                               / ( std::sqrt(3.0) * (3.0 + std::sin(phi) ) )
                               :
                               cohesion * std::cos(phi) + effective_pressure * std::sin(phi) );


              if (stress_ve >= stress_yield)
                viscosities_vep = stress_yield / edot_ii;

            }


          out.viscosities[i] = std::max(viscosities_vep,min_viscosity);
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = compressibility;

          AssertThrow (peridotite_melting_entropy_change < 0.0, ExcMessage ("Double check peridotite melting entropy change value as may be dividing by zero."));
          if (this->get_timestep_number() > 0)
            {
              melt_fraction_pressure_derivative[i] = out.entropy_derivative_pressure[i] / peridotite_melting_entropy_change;
              melt_fraction_temperature_derivative[i] = out.entropy_derivative_temperature[i] / peridotite_melting_entropy_change;
              pressure_material_derivative[i] =  (in.pressure[i] - old_pressures[i]) / dt  + (in.velocity[i] * pressure_gradients[i]);
              temperature_material_derivative[i] = (in.temperature[i] - old_temperatures[i]) / dt + (in.velocity[i] * temperature_gradients[i]);
              melt_production_rates[i] = melt_fraction_pressure_derivative[i] * pressure_material_derivative[i]  
                                   + melt_fraction_temperature_derivative[i] * temperature_material_derivative[i];
            } else {
              melt_fraction_pressure_derivative[i] = 0.0;
              melt_fraction_temperature_derivative[i] = 0.0;
              melt_production_rates[i] = 0.0;
              pressure_material_derivative[i] = 0.0;
              temperature_material_derivative[i] = 0.0;
           }
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      elastic_rheology.fill_elastic_force_outputs(in, average_elastic_shear_moduli, out);
      elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              melt_out->melt_production_rates[i] = melt_production_rates[i];
              melt_out->melt_fraction_temperature_derivative[i] = melt_fraction_temperature_derivative[i];
              melt_out->melt_fraction_pressure_derivative[i] = melt_fraction_pressure_derivative[i];
              melt_out->pressure_material_derivative[i] = pressure_material_derivative[i];
              melt_out->temperature_material_derivative[i] = temperature_material_derivative[i];
              // first, calculate temperature dependence of density
              double temperature_dependence = 1.0;
              if (!constant_density)
                {
                  if (this->include_adiabatic_heating ())
                    {
                      // temperature dependence is 1 - alpha * (T - T(adiabatic))
                      temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                                * thermal_expansivity;
                    }
                  else
                    temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
                }
            }
        }
    }


    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt simpler viscoelastic");
        {
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);
          Rheology::Elasticity<dim>::declare_parameters (prm);

          prm.declare_entry ("Reference solid density", "3000",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "2500",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $\\si{K}$.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference shear modulus","1e10",
                             Patterns::Double (0),
                             "The value of the constant  shear modulus $\\mu_0.$ Units: $Pa$");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa \\, s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Maximum viscosity scaling","1.e3",
                             Patterns::Double (0),
                             "For exponential temperature-dependent shear viscosity. Dimensionless");
          prm.declare_entry ("Minimum viscosity scaling","1.e-1",
                             Patterns::Double (0),
                             "For exponential temperature-dependent shear viscosity. Dimensionless");

          prm.declare_entry ("Maximum bulk viscosity scaling","1.e3",
                             Patterns::Double (0),
                             "For exponential temperature-dependent bulk viscosity. Dimensionless");
          prm.declare_entry ("Minimum bulk viscosity scaling","1.e-1",
                             Patterns::Double (0),
                             "For exponential temperature-dependent bulk viscosity. Dimensionless");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Constant density", "false",
                             Patterns::Bool (),
                             "Turn off density being controlled by thermal expansion"
                             "coefficient. Boolean, default false.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Solid compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the solid matrix. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Melt compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the melt. "
                             "Units: $1/Pa$.");

          prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                             Patterns::Double (0),
                             "The value of the pressure derivative of the melt bulk "
                             "modulus. "
                             "Units: None.");

          prm.declare_entry ("Use full compressibility", "false",
                             Patterns::Bool (),
                             "If the compressibility should be used everywhere in the code "
                             "(if true), changing the volume of material when the density changes, "
                             "or only in the momentum conservation and advection equations "
                             "(if false).");
          prm.declare_entry ("Use fractional melting", "false",
                             Patterns::Bool (),
                             "If fractional melting should be used (if true), including a solidus "
                             "change based on depletion (in this case, the amount of melt that has "
                             "migrated away from its origin), and freezing of melt when it has moved "
                             "to a region with temperatures lower than the solidus; or if batch "
                             "melting should be used (if false), assuming that the melt fraction only "
                             "depends on temperature and pressure, and how much melt has already been "
                             "generated at a given point, but not considering movement of melt in "
                             "the melting parameterization."
                             "\n\n"
                             "Note that melt does not freeze unless the 'Freezing rate' parameter is set "
                             "to a value larger than 0.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of "
                             "depleted material. Depletion is indicated by the compositional "
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Depletion solidus change", "200.0",
                             Patterns::Double (0),
                             "The solidus temperature change for a depletion of 100\\%. For positive "
                             "values, the solidus gets increased for a positive peridotite field "
                             "(depletion) and lowered for a negative peridotite field (enrichment). "
                             "Scaling with depletion is linear. Only active when fractional melting "
                             "is used. "
                             "Units: $\\si{K}$.");
          prm.declare_entry ("Depletion viscosity change","1.0",
                             Patterns::Double (0),
                             "Linear scaling of viscosity, which scales peridotite field.");
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Peridotite melting entropy change", "-300",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of peridotite. "
                             "Units: $J/(kg K)$.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition "),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Reference tensile strength", "1.e20",
                             Patterns::Double (),
                             "Tensile strength value used to calculate melt extraction by diking. Units Pa.");
          prm.declare_entry ("Thermal tensile strength exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the tensile strength. Dimensionless exponent. ");
          prm.declare_entry ("Maximum tensile strength scaling","1.e3",
                             Patterns::Double (0),
                             "For exponential temperature-dependent tensile strength. Dimensionless");
          prm.declare_entry ("Minimum tensile strength scaling","1.e-1",
                             Patterns::Double (0),
                             "For exponential temperature-dependent tensile strength. Dimensionless");
          prm.declare_entry ("Weaken ridge","false",
                             Patterns::Bool (),
                             "Allow region of low viscosity (e.g. at ridge axis)"
                             "in order to get corner flow");
          prm.declare_entry ("Weaken length","0.0",
                             Patterns::Double (),
                             "Length that goes into tanh to determine distance"
                             "from ridge axis where weakening stops");
          prm.declare_entry ("Weaken gradient","1.0",
                             Patterns::Double (),
                             "Gradient that goes into tanh to determine distance"
                             "from ridge axis where weakening stops");
          prm.declare_entry ("Cohesion","1.e50",
                             Patterns::Double (),
                             "Cohesion that determines failure in Mohr-Coulomb crit");
          prm.declare_entry ("Angle of internal friction","30.0",
                             Patterns::Double (),
                             "Phi that determines failure in Mohr-Coulomb crit");
          prm.declare_entry ("Reference strain rate","1.e-15",
                             Patterns::Double (),
                             "Strain rate when model is initialized");
          prm.declare_entry ("Minimum strain rate","1.e-17",
                             Patterns::Double (),
                             "Strain rate can't be lower than this");
          prm.declare_entry ("Activation energy","400.e3",
                             Patterns::Double (),
                             "For temperature/pressure dependence of viscosity. Units: J/mol");
          prm.declare_entry ("Activation volume","1.e-5",
                             Patterns::Double (),
                             "For temperature/pressure dependence of viscosity. Units: m^3");
          prm.declare_entry ("Viscosity scale","1.e-18",
                             Patterns::Double (),
                             "Value to scale the viscosity");
          prm.declare_entry ("Minimum viscosity","1.e15",
                             Patterns::Double (),
                             "Viscosity may not be lower than this value. Units: Pa s");
          prm.declare_entry ("Distance melting suppressed","0.",
                             Patterns::Double (),
                             "Distance from boundaries that melting is suppressed (2-D box). Units: km");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::parse_parameters (ParameterHandler &prm)
    {

      AssertThrow(this->get_parameters().enable_elasticity == true,
                  ExcMessage ("Material model Viscoelastic only works if 'Enable elasticity' is set to true"));

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt simpler viscoelastic");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          reference_T                = prm.get_double ("Reference temperature");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
          max_viscosity_scaling      = prm.get_double ("Maximum viscosity scaling");
          max_bulk_viscosity_scaling       = prm.get_double ("Maximum bulk viscosity scaling");
          min_viscosity_scaling      = prm.get_double ("Minimum viscosity scaling");
          min_bulk_viscosity_scaling       = prm.get_double ("Minimum bulk viscosity scaling");
          thermal_conductivity       = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
          constant_density           = prm.get_bool ("Constant density");
          alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
          compressibility            = prm.get_double ("Solid compressibility");
          melt_compressibility       = prm.get_double ("Melt compressibility");
          model_is_compressible      = prm.get_bool ("Use full compressibility");
          fractional_melting         = prm.get_bool ("Use fractional melting");
          melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
          depletion_density_change   = prm.get_double ("Depletion density change");
          depletion_solidus_change   = prm.get_double ("Depletion solidus change");
          depletion_viscosity_change = prm.get_double ("Depletion viscosity change");
          reference_tensile_strength           = prm.get_double ("Reference tensile strength");
          thermal_tensile_strength_exponent = prm.get_double ("Thermal tensile strength exponent");
          max_tensile_strength_scaling      = prm.get_double ("Maximum tensile strength scaling");
          min_tensile_strength_scaling      = prm.get_double ("Minimum tensile strength scaling");
          weaken_ridge               = prm.get_bool ("Weaken ridge");
          weaken_length              = prm.get_double ("Weaken length");
          weaken_gradient            = prm.get_double ("Weaken gradient");
          phi                        = prm.get_double ("Angle of internal friction");
          cohesion                   = prm.get_double ("Cohesion");
          reference_strain_rate      = prm.get_double ("Reference strain rate");
          minimum_strain_rate      = prm.get_double ("Minimum strain rate");
          activation_energy          = prm.get_double ("Activation energy");
          activation_volume          = prm.get_double ("Activation volume");
          phi *= numbers::PI/180.0;
          viscosity_scale            = prm.get_double ("Viscosity scale");
          min_viscosity              = prm.get_double ("Minimum viscosity");
          distance_melting_suppressed= prm.get_double("Distance melting suppressed");

          if (weaken_ridge && weaken_gradient <= 0.0)
            AssertThrow(false, ExcMessage("Error: Weaken gradient too small."));

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt simpler viscoelastic with Thermal viscosity exponent can not have reference_T=0."));

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          peridotite_melting_entropy_change
            = prm.get_double ("Peridotite melting entropy change");
          M_cpx           = prm.get_double ("Mass fraction cpx");

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          elastic_rheology.initialize_simulator (this->get_simulator());
          elastic_rheology.parse_parameters(prm);

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MeltSimplerViscoelastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      const unsigned int n_points = out.viscosities.size();

      if (this->get_parameters().use_operator_splitting && out.template get_additional_output<ReactionRateOutputs<dim> >() == nullptr)
        {
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }

       if ( out.template get_additional_output<AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
         {
           out.additional_outputs.push_back(
             std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));

         }

       Assert(!this->get_parameters().enable_additional_stokes_rhs
              ||
              out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
              == n_points, ExcInternalError());

      elastic_rheology.create_elastic_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltSimplerViscoelastic,
                                   "melt simpler viscoelastic",
                                   "A material model that takes the implentation of melt production "
                                   "from the melt simple script (dry peridotite solidus of \\cite{KSL2003}, "
                                   "expansion of melt fraction partial derivatives from Dannberg et al. (2016), "
                                   "and employs the viscoelastic rheology by tracking old stresses as compositional fields. ")
  }
}
