/******************************************************************************/
/*           GOLEM - Multiphysics of faulted geothermal reservoirs            */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/

#include "GolemFluidDensityT.h"
#include "GolemGlobals.h"

registerMooseObject("GolemApp", GolemFluidDensityT);

template <>
InputParameters
validParams<GolemFluidDensityT>()
{
  InputParameters params = validParams<GolemFluidDensity>();
  params.addClassDescription("fluid density as function of temperature");
  params.addRequiredParam<Real>("alpha", "the thermal expansion coefficient");
  return params;
}

GolemFluidDensityT::GolemFluidDensityT(const InputParameters & parameters)
  : GolemFluidDensity(parameters), _alpha(getParam<Real>("alpha"))
{
}

Real
GolemFluidDensityT::computeDensity(Real /*pressure*/, Real temperature, Real rho0) const
{
  Real alpha = _alpha;
  if (_has_scaled_properties)
  {
    alpha /= _scaling_uo->_s_expansivity;
  }
    
  return rho0*(1.0 - alpha*temperature); 
}

Real
GolemFluidDensityT::computedDensitydT(Real pressure, Real temperature, Real rho0) const
{
  if (_has_scaled_properties)
  {
    pressure *= _scaling_uo->_s_stress;
    temperature *= _scaling_uo->_s_temperature;
  }
  Real alpha = _alpha;
 
  if (_has_scaled_properties)
  {
    alpha /= _scaling_uo->_s_expansivity;
  }

    return  -1.0*rho0*alpha;
}

Real
GolemFluidDensityT::computedDensitydp(Real /*pressure*/, Real /*temperature*/) const
{
  return 0;
}
