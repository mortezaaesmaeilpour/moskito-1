#ifndef NEWTONRAPHSONTRAINING_2_H
#define NEWTONRAPHSONTRAINING_2_H


#include "Material.h"
#include "NewtonIteration.h"

class NewtonRaphsonTraining_2;

template <>
InputParameters validParams<NewtonRaphsonTraining_2>();

class NewtonRaphsonTraining_2 : public Material, public NewtonIteration
{
public:
  NewtonRaphsonTraining_2(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  virtual Real computeReferenceResidual(const Real iteration_value, const Real scalar) override;
  virtual Real computeResidual(const Real iteration_value, const Real scalar) override;
  virtual Real computeDerivative(const Real iteration_value, const Real scalar) override;
  virtual Real initialGuess(const Real trial_value);
  virtual Real calculate_Reynolds(Real rho, Real v, Real kin_vis, Real length);

  virtual Real minimumPermissibleValue(const Real trail_value) const override;
  virtual Real maximumPermissibleValue(const Real trail_value) const override
  {
    return 100.0;
  }
};
#endif
