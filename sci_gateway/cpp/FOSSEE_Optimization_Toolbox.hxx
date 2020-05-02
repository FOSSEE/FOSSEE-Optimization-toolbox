#ifndef __FOSSEE_OPTIMIZATION_TOOLBOX_GW_HXX__
#define __FOSSEE_OPTIMIZATION_TOOLBOX_GW_HXX__

#ifdef _MSC_VER
#ifdef FOSSEE_OPTIMIZATION_TOOLBOX_GW_EXPORTS
#define FOSSEE_OPTIMIZATION_TOOLBOX_GW_IMPEXP __declspec(dllexport)
#else
#define FOSSEE_OPTIMIZATION_TOOLBOX_GW_IMPEXP __declspec(dllimport)
#endif
#else
#define FOSSEE_OPTIMIZATION_TOOLBOX_GW_IMPEXP
#endif

extern "C" FOSSEE_OPTIMIZATION_TOOLBOX_GW_IMPEXP int FOSSEE_Optimization_Toolbox(wchar_t* _pwstFuncName);
















#endif /* __FOSSEE_OPTIMIZATION_TOOLBOX_GW_HXX__ */
