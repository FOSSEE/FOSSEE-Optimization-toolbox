#include <wchar.h>
#include "FOSSEE_Optimization_Toolbox.hxx"
extern "C"
{
#include "FOSSEE_Optimization_Toolbox.h"
#include "addfunction.h"
}

#define MODULE_NAME L"FOSSEE_Optimization_Toolbox"

int FOSSEE_Optimization_Toolbox(wchar_t* _pwstFuncName)
{
    if(wcscmp(_pwstFuncName, L"linearprog") == 0){ addCFunction(L"linearprog", &sci_linearprog, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"rmps") == 0){ addCFunction(L"rmps", &sci_rmps, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"quadprog_CLP") == 0){ addCFunction(L"quadprog_CLP", &sci_quadprog_CLP, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"solveqp") == 0){ addCFunction(L"solveqp", &sci_solveqp, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"solveminuncp") == 0){ addCFunction(L"solveminuncp", &sci_solveminuncp, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"solveminbndp") == 0){ addCFunction(L"solveminbndp", &sci_solveminbndp, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"solveminconp") == 0){ addCFunction(L"solveminconp", &sci_solveminconp, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"matintlinprog") == 0){ addCFunction(L"matintlinprog", &sci_matintlinprog, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"mpsintlinprog") == 0){ addCFunction(L"mpsintlinprog", &sci_mpsintlinprog, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"solveintqp") == 0){ addCFunction(L"solveintqp", &sci_solveintqp, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"inter_fminunc") == 0){ addCFunction(L"inter_fminunc", &cpp_intfminunc, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"inter_fminbnd") == 0){ addCFunction(L"inter_fminbnd", &cpp_intfminbnd, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"inter_fmincon") == 0){ addCFunction(L"inter_fmincon", &cpp_intfmincon, MODULE_NAME); }
    if(wcscmp(_pwstFuncName, L"fotversion") == 0){ addCFunction(L"fotversion", &sci_fotversion, MODULE_NAME); }

    return 1;
}
