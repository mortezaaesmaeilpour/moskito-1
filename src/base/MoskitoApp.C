#include "MoskitoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<MoskitoApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

MoskitoApp::MoskitoApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  MoskitoApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  MoskitoApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  MoskitoApp::registerExecFlags(_factory);
}

MoskitoApp::~MoskitoApp() {}

void
MoskitoApp::registerApps()
{
  registerApp(MoskitoApp);
}

void
MoskitoApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"MoskitoApp"});
}

void
MoskitoApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"MoskitoApp"});

  /* Uncomment Syntax parameter and register your new production objects here! */
}

void
MoskitoApp::registerObjectDepends(Factory & /*factory*/)
{
}

void
MoskitoApp::associateSyntaxDepends(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
MoskitoApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execution flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
MoskitoApp__registerApps()
{
  MoskitoApp::registerApps();
}

extern "C" void
MoskitoApp__registerObjects(Factory & factory)
{
  MoskitoApp::registerObjects(factory);
}

extern "C" void
MoskitoApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  MoskitoApp::associateSyntax(syntax, action_factory);
}

extern "C" void
MoskitoApp__registerExecFlags(Factory & factory)
{
  MoskitoApp::registerExecFlags(factory);
}
