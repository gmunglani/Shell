#ifndef FUNCTION_H_
#define FUNCTION_H_

#include "libmesh.h"
#include "../muparser/include/muParser.h"

class Function
{
public:
  Function() : value(NULL), is_const(false) {}

  Function(Real& val, const std::string& expr) : value(&val), is_const(false)
  {
    p.SetExpr(expr);
    p.DefineConst("pi", (double)pi);
  }

  inline void add_variable(const std::string& name, Real& var)
  {
    p.DefineVar(name, &var);
  }

  inline void init()
  {
    try
    {
      is_const = (p.GetUsedVar().size() == 0);
    }
    catch (mu::Parser::exception_type& e)
    {
      std::cerr << "ERROR: " << e.GetMsg() << ". Expression: " << e.GetExpr() << std::endl;
      exit(1);
    }
  }

  inline bool is_constant()
  {
    return is_const;
  }

  inline void evaluate()
  {
    libmesh_assert(value);
    (*value) = p.Eval();
  }
  
protected:
  mu::Parser p;
  Real* value;
  bool is_const;
};

#endif /* FUNCTION_H_ */
