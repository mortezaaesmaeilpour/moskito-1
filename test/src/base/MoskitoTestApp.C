//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "MoskitoTestApp.h"
#include "MoskitoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<MoskitoTestApp>()
{
  InputParameters params = validParams<MoskitoApp>();
  return params;
}

MoskitoTestApp::MoskitoTestApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  MoskitoApp::registerObjectDepends(_factory);
  MoskitoApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  MoskitoApp::associateSyntaxDepends(_syntax, _action_factory);
  MoskitoApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  MoskitoApp::registerExecFlags(_factory);

  bool use_test_objs = getParam<bool>("allow_test_objects");
  if (use_test_objs)
  {
    MoskitoTestApp::registerObjects(_factory);
    MoskitoTestApp::associateSyntax(_syntax, _action_factory);
    MoskitoTestApp::registerExecFlags(_factory);
  }
}

MoskitoTestApp::~MoskitoTestApp() {}

void
MoskitoTestApp::registerApps()
{
  registerApp(MoskitoApp);
  registerApp(MoskitoTestApp);
}

void
MoskitoTestApp::registerObjects(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new test objects here! */
}

void
MoskitoTestApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
  /* Uncomment Syntax and ActionFactory parameters and register your new test objects here! */
}

void
MoskitoTestApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execute flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
MoskitoTestApp__registerApps()
{
  MoskitoTestApp::registerApps();
}

// External entry point for dynamic object registration
extern "C" void
MoskitoTestApp__registerObjects(Factory & factory)
{
  MoskitoTestApp::registerObjects(factory);
}

// External entry point for dynamic syntax association
extern "C" void
MoskitoTestApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  MoskitoTestApp::associateSyntax(syntax, action_factory);
}

// External entry point for dynamic execute flag loading
extern "C" void
MoskitoTestApp__registerExecFlags(Factory & factory)
{
  MoskitoTestApp::registerExecFlags(factory);
}
