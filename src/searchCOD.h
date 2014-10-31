#ifndef _crossover_SEARCH_H
#define _crossover_SEARCH_H

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <R.h>
#include <Rmath.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

RcppExport SEXP searchCOD(SEXP sS, SEXP pS, SEXP vS, SEXP designS, SEXP linkMS, SEXP contrastS, SEXP modelS, SEXP vRepS, SEXP balanceSS, SEXP balancePS, SEXP verboseS, SEXP nS, SEXP jumpS, SEXP s2, SEXP checkES, SEXP randomSS, SEXP correlationS, SEXP interchangeS);

RcppExport SEXP rcd2R(SEXP designS, SEXP vS, SEXP modelS);
RcppExport SEXP rcdMatrix2R(SEXP designS, SEXP vS, SEXP modelS);
RcppExport SEXP infMatrix2R(SEXP designS, SEXP vS, SEXP modelS, SEXP xtxS);
RcppExport SEXP estimable2R(SEXP rcDesignS, SEXP vS, SEXP modelS, SEXP linkMS, SEXP contrastS, SEXP ZS, SEXP verboseS);
RcppExport SEXP getS12R(SEXP designS, SEXP vS, SEXP modelS, SEXP linkMS, SEXP contrastS, SEXP randomSS);
RcppExport SEXP getZ2R(SEXP sS, SEXP pS, SEXP randomSS);
RcppExport SEXP designMatrix2R(SEXP designS, SEXP vS, SEXP modelS);

arma::mat rcd(arma::mat design, int v, int model);
arma::mat rcdMatrix(arma::mat rcDesign, int v, int model);
arma::mat infMatrix(arma::mat rcDesign, int v, int model, bool xtx);
arma::mat designMatrix(arma::mat design, int v, int model);
double getS1(arma::mat rcDesign, int v, int model, arma::mat linkM, arma::mat tCC, bool randomS);
bool estimable(arma::mat rcDesign, int v, int model, arma::mat linkM, arma::mat C, arma::mat Z, int verbose);
arma::mat getZ(int s, int p, bool randomS);

#endif
