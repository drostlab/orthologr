#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _orthologr_ambigousNucleotides(SEXP);
extern SEXP _orthologr_codonPrecondition(SEXP);
extern SEXP _orthologr_Different(SEXP, SEXP, SEXP, SEXP);
extern SEXP _orthologr_gestimator(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _orthologr_intToNuc(SEXP);
extern SEXP _orthologr_NotAGap(SEXP);
extern SEXP _orthologr_nucToInt(SEXP);
extern SEXP _orthologr_NumDiffs(SEXP, SEXP, SEXP, SEXP);
extern SEXP _orthologr_toChar(SEXP);
extern SEXP _orthologr_TranslateCodon(SEXP);
extern SEXP _orthologr_TsTv(SEXP, SEXP);
extern SEXP _orthologr_Universal(SEXP);

static const R_CallMethodDef CallEntries[] = {
        {"_orthologr_ambigousNucleotides", (DL_FUNC) &_orthologr_ambigousNucleotides, 1},
        {"_orthologr_codonPrecondition",   (DL_FUNC) &_orthologr_codonPrecondition,   1},
        {"_orthologr_Different",           (DL_FUNC) &_orthologr_Different,           4},
        {"_orthologr_gestimator",          (DL_FUNC) &_orthologr_gestimator,          5},
        {"_orthologr_intToNuc",            (DL_FUNC) &_orthologr_intToNuc,            1},
        {"_orthologr_NotAGap",             (DL_FUNC) &_orthologr_NotAGap,             1},
        {"_orthologr_nucToInt",            (DL_FUNC) &_orthologr_nucToInt,            1},
        {"_orthologr_NumDiffs",            (DL_FUNC) &_orthologr_NumDiffs,            4},
        {"_orthologr_toChar",              (DL_FUNC) &_orthologr_toChar,              1},
        {"_orthologr_TranslateCodon",      (DL_FUNC) &_orthologr_TranslateCodon,      1},
        {"_orthologr_TsTv",                (DL_FUNC) &_orthologr_TsTv,                2},
        {"_orthologr_Universal",           (DL_FUNC) &_orthologr_Universal,           1},
        {NULL, NULL, 0}
};

void R_init_orthologr(DllInfo *dll)
{
        R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
}