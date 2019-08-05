#ifndef NEWTONRAPHSONTRAINING_1_H
#define NEWTONRAPHSONTRAINING_1_H


#include "Material.h"
#include "NewtonIteration.h"

class NewtonRaphsonTraining_1;

template <>
InputParameters validParams<NewtonRaphsonTraining_1>();

class NewtonRaphsonTraining_1 : public Material, public NewtonIteration
{
public:
  NewtonRaphsonTraining_1(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  virtual Real computeReferenceResidual(const Real iteration_value, const Real scalar) override;
  virtual Real computeResidual(const Real iteration_value, const Real scalar) override;
  virtual Real computeDerivative(const Real iteration_value, const Real scalar) override;
  virtual Real initialGuess(const Real trial_value);
};
#endif
