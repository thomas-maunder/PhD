/* File: _interfacemodule.c
 * This file is auto-generated with f2py (version:1.21.4).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Wed Dec  1 02:56:24 2021
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *_interface_error;
static PyObject *_interface_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f##_,F##_)
#else
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f,F)
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj, 0);
    if (tmp) {
        PyErr_Clear();
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = _interface_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}

static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }
    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj,0);
    if (tmp) {
        PyErr_Clear();
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = _interface_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_WRAPPEDFUNC_US(mapping_,MAPPING_)(int*,int*,int*,double*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** mapping_ **********************************/
static char doc_f2py_rout__interface_mapping_[] = "\
mapping_(nr,nth,nph,yywt,nx,ny,nz,nnuc,x0,y0,z0,rho,xnu,r3c,r3l,r3r,dx,dy,dz,dphi,rhonew,xnunew,muc,mul,mur,phic,phil,phir,tol)\n\nWrapper for ``mapping_``.\
\n\nParameters\n----------\n"
"nr : input int\n"
"nth : input int\n"
"nph : input int\n"
"yywt : input rank-2 array('d') with bounds (f2py_yywt_d0,f2py_yywt_d1)\n"
"nx : input int\n"
"ny : input int\n"
"nz : input int\n"
"nnuc : input int\n"
"x0 : input float\n"
"y0 : input float\n"
"z0 : input float\n"
"rho : input rank-3 array('d') with bounds (f2py_rho_d0,f2py_rho_d1,f2py_rho_d2)\n"
"xnu : input rank-4 array('d') with bounds (f2py_xnu_d0,f2py_xnu_d1,f2py_xnu_d2,f2py_xnu_d3)\n"
"r3c : input rank-1 array('d') with bounds (f2py_r3c_d0)\n"
"r3l : input rank-1 array('d') with bounds (f2py_r3l_d0)\n"
"r3r : input rank-1 array('d') with bounds (f2py_r3r_d0)\n"
"dx : input float\n"
"dy : input float\n"
"dz : input float\n"
"dphi : input float\n"
"rhonew : in/output rank-3 array('d') with bounds (f2py_rhonew_d0,f2py_rhonew_d1,f2py_rhonew_d2)\n"
"xnunew : in/output rank-4 array('d') with bounds (f2py_xnunew_d0,f2py_xnunew_d1,f2py_xnunew_d2,f2py_xnunew_d3)\n"
"muc : input rank-1 array('d') with bounds (f2py_muc_d0)\n"
"mul : input rank-1 array('d') with bounds (f2py_mul_d0)\n"
"mur : input rank-1 array('d') with bounds (f2py_mur_d0)\n"
"phic : input rank-1 array('d') with bounds (f2py_phic_d0)\n"
"phil : input rank-1 array('d') with bounds (f2py_phil_d0)\n"
"phir : input rank-1 array('d') with bounds (f2py_phir_d0)\n"
"tol : input float";
/* extern void F_WRAPPEDFUNC_US(mapping_,MAPPING_)(int*,int*,int*,double*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*); */
static PyObject *f2py_rout__interface_mapping_(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,int*,double*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

  int nr = 0;
  PyObject *nr_capi = Py_None;
  int nth = 0;
  PyObject *nth_capi = Py_None;
  int nph = 0;
  PyObject *nph_capi = Py_None;
  double *yywt = NULL;
  npy_intp yywt_Dims[2] = {-1, -1};
  const int yywt_Rank = 2;
  PyArrayObject *capi_yywt_tmp = NULL;
  int capi_yywt_intent = 0;
  PyObject *yywt_capi = Py_None;
  int nx = 0;
  PyObject *nx_capi = Py_None;
  int ny = 0;
  PyObject *ny_capi = Py_None;
  int nz = 0;
  PyObject *nz_capi = Py_None;
  int nnuc = 0;
  PyObject *nnuc_capi = Py_None;
  double x0 = 0;
  PyObject *x0_capi = Py_None;
  double y0 = 0;
  PyObject *y0_capi = Py_None;
  double z0 = 0;
  PyObject *z0_capi = Py_None;
  double *rho = NULL;
  npy_intp rho_Dims[3] = {-1, -1, -1};
  const int rho_Rank = 3;
  PyArrayObject *capi_rho_tmp = NULL;
  int capi_rho_intent = 0;
  PyObject *rho_capi = Py_None;
  double *xnu = NULL;
  npy_intp xnu_Dims[4] = {-1, -1, -1, -1};
  const int xnu_Rank = 4;
  PyArrayObject *capi_xnu_tmp = NULL;
  int capi_xnu_intent = 0;
  PyObject *xnu_capi = Py_None;
  double *r3c = NULL;
  npy_intp r3c_Dims[1] = {-1};
  const int r3c_Rank = 1;
  PyArrayObject *capi_r3c_tmp = NULL;
  int capi_r3c_intent = 0;
  PyObject *r3c_capi = Py_None;
  double *r3l = NULL;
  npy_intp r3l_Dims[1] = {-1};
  const int r3l_Rank = 1;
  PyArrayObject *capi_r3l_tmp = NULL;
  int capi_r3l_intent = 0;
  PyObject *r3l_capi = Py_None;
  double *r3r = NULL;
  npy_intp r3r_Dims[1] = {-1};
  const int r3r_Rank = 1;
  PyArrayObject *capi_r3r_tmp = NULL;
  int capi_r3r_intent = 0;
  PyObject *r3r_capi = Py_None;
  double dx = 0;
  PyObject *dx_capi = Py_None;
  double dy = 0;
  PyObject *dy_capi = Py_None;
  double dz = 0;
  PyObject *dz_capi = Py_None;
  double dphi = 0;
  PyObject *dphi_capi = Py_None;
  double *rhonew = NULL;
  npy_intp rhonew_Dims[3] = {-1, -1, -1};
  const int rhonew_Rank = 3;
  PyArrayObject *capi_rhonew_tmp = NULL;
  int capi_rhonew_intent = 0;
  PyObject *rhonew_capi = Py_None;
  double *xnunew = NULL;
  npy_intp xnunew_Dims[4] = {-1, -1, -1, -1};
  const int xnunew_Rank = 4;
  PyArrayObject *capi_xnunew_tmp = NULL;
  int capi_xnunew_intent = 0;
  PyObject *xnunew_capi = Py_None;
  double *muc = NULL;
  npy_intp muc_Dims[1] = {-1};
  const int muc_Rank = 1;
  PyArrayObject *capi_muc_tmp = NULL;
  int capi_muc_intent = 0;
  PyObject *muc_capi = Py_None;
  double *mul = NULL;
  npy_intp mul_Dims[1] = {-1};
  const int mul_Rank = 1;
  PyArrayObject *capi_mul_tmp = NULL;
  int capi_mul_intent = 0;
  PyObject *mul_capi = Py_None;
  double *mur = NULL;
  npy_intp mur_Dims[1] = {-1};
  const int mur_Rank = 1;
  PyArrayObject *capi_mur_tmp = NULL;
  int capi_mur_intent = 0;
  PyObject *mur_capi = Py_None;
  double *phic = NULL;
  npy_intp phic_Dims[1] = {-1};
  const int phic_Rank = 1;
  PyArrayObject *capi_phic_tmp = NULL;
  int capi_phic_intent = 0;
  PyObject *phic_capi = Py_None;
  double *phil = NULL;
  npy_intp phil_Dims[1] = {-1};
  const int phil_Rank = 1;
  PyArrayObject *capi_phil_tmp = NULL;
  int capi_phil_intent = 0;
  PyObject *phil_capi = Py_None;
  double *phir = NULL;
  npy_intp phir_Dims[1] = {-1};
  const int phir_Rank = 1;
  PyArrayObject *capi_phir_tmp = NULL;
  int capi_phir_intent = 0;
  PyObject *phir_capi = Py_None;
  double tol = 0;
  PyObject *tol_capi = Py_None;
  int f2py_yywt_d0 = 0;
  int f2py_yywt_d1 = 0;
  int f2py_rho_d0 = 0;
  int f2py_rho_d1 = 0;
  int f2py_rho_d2 = 0;
  int f2py_xnu_d0 = 0;
  int f2py_xnu_d1 = 0;
  int f2py_xnu_d2 = 0;
  int f2py_xnu_d3 = 0;
  int f2py_r3c_d0 = 0;
  int f2py_r3l_d0 = 0;
  int f2py_r3r_d0 = 0;
  int f2py_rhonew_d0 = 0;
  int f2py_rhonew_d1 = 0;
  int f2py_rhonew_d2 = 0;
  int f2py_xnunew_d0 = 0;
  int f2py_xnunew_d1 = 0;
  int f2py_xnunew_d2 = 0;
  int f2py_xnunew_d3 = 0;
  int f2py_muc_d0 = 0;
  int f2py_mul_d0 = 0;
  int f2py_mur_d0 = 0;
  int f2py_phic_d0 = 0;
  int f2py_phil_d0 = 0;
  int f2py_phir_d0 = 0;
    static char *capi_kwlist[] = {"nr","nth","nph","yywt","nx","ny","nz","nnuc","x0","y0","z0","rho","xnu","r3c","r3l","r3r","dx","dy","dz","dphi","rhonew","xnunew","muc","mul","mur","phic","phil","phir","tol",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOOOOOOOOOOOOOOOOOOOOOOOOOO|:_interface.mapping_",\
        capi_kwlist,&nr_capi,&nth_capi,&nph_capi,&yywt_capi,&nx_capi,&ny_capi,&nz_capi,&nnuc_capi,&x0_capi,&y0_capi,&z0_capi,&rho_capi,&xnu_capi,&r3c_capi,&r3l_capi,&r3r_capi,&dx_capi,&dy_capi,&dz_capi,&dphi_capi,&rhonew_capi,&xnunew_capi,&muc_capi,&mul_capi,&mur_capi,&phic_capi,&phil_capi,&phir_capi,&tol_capi))
        return NULL;
/*frompyobj*/
  /* Processing variable nr */
    f2py_success = int_from_pyobj(&nr,nr_capi,"_interface.mapping_() 1st argument (nr) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nth */
    f2py_success = int_from_pyobj(&nth,nth_capi,"_interface.mapping_() 2nd argument (nth) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nph */
    f2py_success = int_from_pyobj(&nph,nph_capi,"_interface.mapping_() 3rd argument (nph) can't be converted to int");
  if (f2py_success) {
  /* Processing variable yywt */
  ;
  capi_yywt_intent |= F2PY_INTENT_IN;
  capi_yywt_tmp = array_from_pyobj(NPY_DOUBLE,yywt_Dims,yywt_Rank,capi_yywt_intent,yywt_capi);
  if (capi_yywt_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 4th argument `yywt' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    yywt = (double *)(PyArray_DATA(capi_yywt_tmp));

  /* Processing variable nx */
    f2py_success = int_from_pyobj(&nx,nx_capi,"_interface.mapping_() 5th argument (nx) can't be converted to int");
  if (f2py_success) {
  /* Processing variable ny */
    f2py_success = int_from_pyobj(&ny,ny_capi,"_interface.mapping_() 6th argument (ny) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nz */
    f2py_success = int_from_pyobj(&nz,nz_capi,"_interface.mapping_() 7th argument (nz) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nnuc */
    f2py_success = int_from_pyobj(&nnuc,nnuc_capi,"_interface.mapping_() 8th argument (nnuc) can't be converted to int");
  if (f2py_success) {
  /* Processing variable x0 */
    f2py_success = double_from_pyobj(&x0,x0_capi,"_interface.mapping_() 9th argument (x0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable y0 */
    f2py_success = double_from_pyobj(&y0,y0_capi,"_interface.mapping_() 10th argument (y0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable z0 */
    f2py_success = double_from_pyobj(&z0,z0_capi,"_interface.mapping_() 11st argument (z0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable rho */
  ;
  capi_rho_intent |= F2PY_INTENT_IN;
  capi_rho_tmp = array_from_pyobj(NPY_DOUBLE,rho_Dims,rho_Rank,capi_rho_intent,rho_capi);
  if (capi_rho_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 12nd argument `rho' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    rho = (double *)(PyArray_DATA(capi_rho_tmp));

  /* Processing variable xnu */
  ;
  capi_xnu_intent |= F2PY_INTENT_IN;
  capi_xnu_tmp = array_from_pyobj(NPY_DOUBLE,xnu_Dims,xnu_Rank,capi_xnu_intent,xnu_capi);
  if (capi_xnu_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 13rd argument `xnu' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    xnu = (double *)(PyArray_DATA(capi_xnu_tmp));

  /* Processing variable r3c */
  ;
  capi_r3c_intent |= F2PY_INTENT_IN;
  capi_r3c_tmp = array_from_pyobj(NPY_DOUBLE,r3c_Dims,r3c_Rank,capi_r3c_intent,r3c_capi);
  if (capi_r3c_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 14th argument `r3c' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    r3c = (double *)(PyArray_DATA(capi_r3c_tmp));

  /* Processing variable r3l */
  ;
  capi_r3l_intent |= F2PY_INTENT_IN;
  capi_r3l_tmp = array_from_pyobj(NPY_DOUBLE,r3l_Dims,r3l_Rank,capi_r3l_intent,r3l_capi);
  if (capi_r3l_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 15th argument `r3l' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    r3l = (double *)(PyArray_DATA(capi_r3l_tmp));

  /* Processing variable r3r */
  ;
  capi_r3r_intent |= F2PY_INTENT_IN;
  capi_r3r_tmp = array_from_pyobj(NPY_DOUBLE,r3r_Dims,r3r_Rank,capi_r3r_intent,r3r_capi);
  if (capi_r3r_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 16th argument `r3r' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    r3r = (double *)(PyArray_DATA(capi_r3r_tmp));

  /* Processing variable dx */
    f2py_success = double_from_pyobj(&dx,dx_capi,"_interface.mapping_() 17th argument (dx) can't be converted to double");
  if (f2py_success) {
  /* Processing variable dy */
    f2py_success = double_from_pyobj(&dy,dy_capi,"_interface.mapping_() 18th argument (dy) can't be converted to double");
  if (f2py_success) {
  /* Processing variable dz */
    f2py_success = double_from_pyobj(&dz,dz_capi,"_interface.mapping_() 19th argument (dz) can't be converted to double");
  if (f2py_success) {
  /* Processing variable dphi */
    f2py_success = double_from_pyobj(&dphi,dphi_capi,"_interface.mapping_() 20th argument (dphi) can't be converted to double");
  if (f2py_success) {
  /* Processing variable rhonew */
  ;
  capi_rhonew_intent |= F2PY_INTENT_INOUT;
  capi_rhonew_tmp = array_from_pyobj(NPY_DOUBLE,rhonew_Dims,rhonew_Rank,capi_rhonew_intent,rhonew_capi);
  if (capi_rhonew_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 21st argument `rhonew' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    rhonew = (double *)(PyArray_DATA(capi_rhonew_tmp));

  /* Processing variable xnunew */
  ;
  capi_xnunew_intent |= F2PY_INTENT_INOUT;
  capi_xnunew_tmp = array_from_pyobj(NPY_DOUBLE,xnunew_Dims,xnunew_Rank,capi_xnunew_intent,xnunew_capi);
  if (capi_xnunew_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 22nd argument `xnunew' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    xnunew = (double *)(PyArray_DATA(capi_xnunew_tmp));

  /* Processing variable muc */
  ;
  capi_muc_intent |= F2PY_INTENT_IN;
  capi_muc_tmp = array_from_pyobj(NPY_DOUBLE,muc_Dims,muc_Rank,capi_muc_intent,muc_capi);
  if (capi_muc_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 23rd argument `muc' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    muc = (double *)(PyArray_DATA(capi_muc_tmp));

  /* Processing variable mul */
  ;
  capi_mul_intent |= F2PY_INTENT_IN;
  capi_mul_tmp = array_from_pyobj(NPY_DOUBLE,mul_Dims,mul_Rank,capi_mul_intent,mul_capi);
  if (capi_mul_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 24th argument `mul' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    mul = (double *)(PyArray_DATA(capi_mul_tmp));

  /* Processing variable mur */
  ;
  capi_mur_intent |= F2PY_INTENT_IN;
  capi_mur_tmp = array_from_pyobj(NPY_DOUBLE,mur_Dims,mur_Rank,capi_mur_intent,mur_capi);
  if (capi_mur_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 25th argument `mur' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    mur = (double *)(PyArray_DATA(capi_mur_tmp));

  /* Processing variable phic */
  ;
  capi_phic_intent |= F2PY_INTENT_IN;
  capi_phic_tmp = array_from_pyobj(NPY_DOUBLE,phic_Dims,phic_Rank,capi_phic_intent,phic_capi);
  if (capi_phic_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 26th argument `phic' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    phic = (double *)(PyArray_DATA(capi_phic_tmp));

  /* Processing variable phil */
  ;
  capi_phil_intent |= F2PY_INTENT_IN;
  capi_phil_tmp = array_from_pyobj(NPY_DOUBLE,phil_Dims,phil_Rank,capi_phil_intent,phil_capi);
  if (capi_phil_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 27th argument `phil' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    phil = (double *)(PyArray_DATA(capi_phil_tmp));

  /* Processing variable phir */
  ;
  capi_phir_intent |= F2PY_INTENT_IN;
  capi_phir_tmp = array_from_pyobj(NPY_DOUBLE,phir_Dims,phir_Rank,capi_phir_intent,phir_capi);
  if (capi_phir_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : _interface_error,"failed in converting 28th argument `phir' of _interface.mapping_ to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    phir = (double *)(PyArray_DATA(capi_phir_tmp));

  /* Processing variable tol */
    f2py_success = double_from_pyobj(&tol,tol_capi,"_interface.mapping_() 29th argument (tol) can't be converted to double");
  if (f2py_success) {
  /* Processing variable f2py_yywt_d0 */
  f2py_yywt_d0 = shape(yywt, 0);
  /* Processing variable f2py_yywt_d1 */
  f2py_yywt_d1 = shape(yywt, 1);
  /* Processing variable f2py_rho_d0 */
  f2py_rho_d0 = shape(rho, 0);
  /* Processing variable f2py_rho_d1 */
  f2py_rho_d1 = shape(rho, 1);
  /* Processing variable f2py_rho_d2 */
  f2py_rho_d2 = shape(rho, 2);
  /* Processing variable f2py_xnu_d0 */
  f2py_xnu_d0 = shape(xnu, 0);
  /* Processing variable f2py_xnu_d1 */
  f2py_xnu_d1 = shape(xnu, 1);
  /* Processing variable f2py_xnu_d2 */
  f2py_xnu_d2 = shape(xnu, 2);
  /* Processing variable f2py_xnu_d3 */
  f2py_xnu_d3 = shape(xnu, 3);
  /* Processing variable f2py_r3c_d0 */
  f2py_r3c_d0 = shape(r3c, 0);
  /* Processing variable f2py_r3l_d0 */
  f2py_r3l_d0 = shape(r3l, 0);
  /* Processing variable f2py_r3r_d0 */
  f2py_r3r_d0 = shape(r3r, 0);
  /* Processing variable f2py_rhonew_d0 */
  f2py_rhonew_d0 = shape(rhonew, 0);
  /* Processing variable f2py_rhonew_d1 */
  f2py_rhonew_d1 = shape(rhonew, 1);
  /* Processing variable f2py_rhonew_d2 */
  f2py_rhonew_d2 = shape(rhonew, 2);
  /* Processing variable f2py_xnunew_d0 */
  f2py_xnunew_d0 = shape(xnunew, 0);
  /* Processing variable f2py_xnunew_d1 */
  f2py_xnunew_d1 = shape(xnunew, 1);
  /* Processing variable f2py_xnunew_d2 */
  f2py_xnunew_d2 = shape(xnunew, 2);
  /* Processing variable f2py_xnunew_d3 */
  f2py_xnunew_d3 = shape(xnunew, 3);
  /* Processing variable f2py_muc_d0 */
  f2py_muc_d0 = shape(muc, 0);
  /* Processing variable f2py_mul_d0 */
  f2py_mul_d0 = shape(mul, 0);
  /* Processing variable f2py_mur_d0 */
  f2py_mur_d0 = shape(mur, 0);
  /* Processing variable f2py_phic_d0 */
  f2py_phic_d0 = shape(phic, 0);
  /* Processing variable f2py_phil_d0 */
  f2py_phil_d0 = shape(phil, 0);
  /* Processing variable f2py_phir_d0 */
  f2py_phir_d0 = shape(phir, 0);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(&nr,&nth,&nph,yywt,&nx,&ny,&nz,&nnuc,&x0,&y0,&z0,rho,xnu,r3c,r3l,r3r,&dx,&dy,&dz,&dphi,rhonew,xnunew,muc,mul,mur,phic,phil,phir,&tol,&f2py_yywt_d0,&f2py_yywt_d1,&f2py_rho_d0,&f2py_rho_d1,&f2py_rho_d2,&f2py_xnu_d0,&f2py_xnu_d1,&f2py_xnu_d2,&f2py_xnu_d3,&f2py_r3c_d0,&f2py_r3l_d0,&f2py_r3r_d0,&f2py_rhonew_d0,&f2py_rhonew_d1,&f2py_rhonew_d2,&f2py_xnunew_d0,&f2py_xnunew_d1,&f2py_xnunew_d2,&f2py_xnunew_d3,&f2py_muc_d0,&f2py_mul_d0,&f2py_mur_d0,&f2py_phic_d0,&f2py_phil_d0,&f2py_phir_d0);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  /* End of cleaning variable f2py_phir_d0 */
  /* End of cleaning variable f2py_phil_d0 */
  /* End of cleaning variable f2py_phic_d0 */
  /* End of cleaning variable f2py_mur_d0 */
  /* End of cleaning variable f2py_mul_d0 */
  /* End of cleaning variable f2py_muc_d0 */
  /* End of cleaning variable f2py_xnunew_d3 */
  /* End of cleaning variable f2py_xnunew_d2 */
  /* End of cleaning variable f2py_xnunew_d1 */
  /* End of cleaning variable f2py_xnunew_d0 */
  /* End of cleaning variable f2py_rhonew_d2 */
  /* End of cleaning variable f2py_rhonew_d1 */
  /* End of cleaning variable f2py_rhonew_d0 */
  /* End of cleaning variable f2py_r3r_d0 */
  /* End of cleaning variable f2py_r3l_d0 */
  /* End of cleaning variable f2py_r3c_d0 */
  /* End of cleaning variable f2py_xnu_d3 */
  /* End of cleaning variable f2py_xnu_d2 */
  /* End of cleaning variable f2py_xnu_d1 */
  /* End of cleaning variable f2py_xnu_d0 */
  /* End of cleaning variable f2py_rho_d2 */
  /* End of cleaning variable f2py_rho_d1 */
  /* End of cleaning variable f2py_rho_d0 */
  /* End of cleaning variable f2py_yywt_d1 */
  /* End of cleaning variable f2py_yywt_d0 */
  } /*if (f2py_success) of tol*/
  /* End of cleaning variable tol */
  if((PyObject *)capi_phir_tmp!=phir_capi) {
    Py_XDECREF(capi_phir_tmp); }
  }  /*if (capi_phir_tmp == NULL) ... else of phir*/
  /* End of cleaning variable phir */
  if((PyObject *)capi_phil_tmp!=phil_capi) {
    Py_XDECREF(capi_phil_tmp); }
  }  /*if (capi_phil_tmp == NULL) ... else of phil*/
  /* End of cleaning variable phil */
  if((PyObject *)capi_phic_tmp!=phic_capi) {
    Py_XDECREF(capi_phic_tmp); }
  }  /*if (capi_phic_tmp == NULL) ... else of phic*/
  /* End of cleaning variable phic */
  if((PyObject *)capi_mur_tmp!=mur_capi) {
    Py_XDECREF(capi_mur_tmp); }
  }  /*if (capi_mur_tmp == NULL) ... else of mur*/
  /* End of cleaning variable mur */
  if((PyObject *)capi_mul_tmp!=mul_capi) {
    Py_XDECREF(capi_mul_tmp); }
  }  /*if (capi_mul_tmp == NULL) ... else of mul*/
  /* End of cleaning variable mul */
  if((PyObject *)capi_muc_tmp!=muc_capi) {
    Py_XDECREF(capi_muc_tmp); }
  }  /*if (capi_muc_tmp == NULL) ... else of muc*/
  /* End of cleaning variable muc */
  if((PyObject *)capi_xnunew_tmp!=xnunew_capi) {
    Py_XDECREF(capi_xnunew_tmp); }
  }  /*if (capi_xnunew_tmp == NULL) ... else of xnunew*/
  /* End of cleaning variable xnunew */
  if((PyObject *)capi_rhonew_tmp!=rhonew_capi) {
    Py_XDECREF(capi_rhonew_tmp); }
  }  /*if (capi_rhonew_tmp == NULL) ... else of rhonew*/
  /* End of cleaning variable rhonew */
  } /*if (f2py_success) of dphi*/
  /* End of cleaning variable dphi */
  } /*if (f2py_success) of dz*/
  /* End of cleaning variable dz */
  } /*if (f2py_success) of dy*/
  /* End of cleaning variable dy */
  } /*if (f2py_success) of dx*/
  /* End of cleaning variable dx */
  if((PyObject *)capi_r3r_tmp!=r3r_capi) {
    Py_XDECREF(capi_r3r_tmp); }
  }  /*if (capi_r3r_tmp == NULL) ... else of r3r*/
  /* End of cleaning variable r3r */
  if((PyObject *)capi_r3l_tmp!=r3l_capi) {
    Py_XDECREF(capi_r3l_tmp); }
  }  /*if (capi_r3l_tmp == NULL) ... else of r3l*/
  /* End of cleaning variable r3l */
  if((PyObject *)capi_r3c_tmp!=r3c_capi) {
    Py_XDECREF(capi_r3c_tmp); }
  }  /*if (capi_r3c_tmp == NULL) ... else of r3c*/
  /* End of cleaning variable r3c */
  if((PyObject *)capi_xnu_tmp!=xnu_capi) {
    Py_XDECREF(capi_xnu_tmp); }
  }  /*if (capi_xnu_tmp == NULL) ... else of xnu*/
  /* End of cleaning variable xnu */
  if((PyObject *)capi_rho_tmp!=rho_capi) {
    Py_XDECREF(capi_rho_tmp); }
  }  /*if (capi_rho_tmp == NULL) ... else of rho*/
  /* End of cleaning variable rho */
  } /*if (f2py_success) of z0*/
  /* End of cleaning variable z0 */
  } /*if (f2py_success) of y0*/
  /* End of cleaning variable y0 */
  } /*if (f2py_success) of x0*/
  /* End of cleaning variable x0 */
  } /*if (f2py_success) of nnuc*/
  /* End of cleaning variable nnuc */
  } /*if (f2py_success) of nz*/
  /* End of cleaning variable nz */
  } /*if (f2py_success) of ny*/
  /* End of cleaning variable ny */
  } /*if (f2py_success) of nx*/
  /* End of cleaning variable nx */
  if((PyObject *)capi_yywt_tmp!=yywt_capi) {
    Py_XDECREF(capi_yywt_tmp); }
  }  /*if (capi_yywt_tmp == NULL) ... else of yywt*/
  /* End of cleaning variable yywt */
  } /*if (f2py_success) of nph*/
  /* End of cleaning variable nph */
  } /*if (f2py_success) of nth*/
  /* End of cleaning variable nth */
  } /*if (f2py_success) of nr*/
  /* End of cleaning variable nr */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/****************************** end of mapping_ ******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"mapping_",-1,{{-1}},0,(char *)F_WRAPPEDFUNC_US(mapping_,MAPPING_),(f2py_init_func)f2py_rout__interface_mapping_,doc_f2py_rout__interface_mapping_},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_interface",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit__interface(void) {
  int i;
  PyObject *m,*d, *s, *tmp;
  m = _interface_module = PyModule_Create(&moduledef);
  Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module _interface (failed to import numpy)"); return m;}
  d = PyModule_GetDict(m);
  s = PyUnicode_FromString("1.21.4");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString(
    "This module '_interface' is auto-generated with f2py (version:1.21.4).\nFunctions:\n"
"  mapping_(nr,nth,nph,yywt,nx,ny,nz,nnuc,x0,y0,z0,rho,xnu,r3c,r3l,r3r,dx,dy,dz,dphi,rhonew,xnunew,muc,mul,mur,phic,phil,phir,tol)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString("1.21.4");
  PyDict_SetItemString(d, "__f2py_numpy_version__", s);
  Py_DECREF(s);
  _interface_error = PyErr_NewException ("_interface.error", NULL, NULL);
  /*
   * Store the error object inside the dict, so that it could get deallocated.
   * (in practice, this is a module, so it likely will not and cannot.)
   */
  PyDict_SetItemString(d, "__interface_error", _interface_error);
  Py_DECREF(_interface_error);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
    tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
    PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
    Py_DECREF(tmp);
  }

    {
      extern void F_FUNC_US(mapping_,MAPPING_)(void);
      PyObject* o = PyDict_GetItemString(d,"mapping_");
      tmp = F2PyCapsule_FromVoidPtr((void*)F_FUNC_US(mapping_,MAPPING_),NULL);
      PyObject_SetAttrString(o,"_cpointer", tmp);
      Py_DECREF(tmp);
      s = PyUnicode_FromString("mapping_");
      PyObject_SetAttrString(o,"__name__", s);
      Py_DECREF(s);
    }
    
/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"_interface");
#endif
  return m;
}
#ifdef __cplusplus
}
#endif
