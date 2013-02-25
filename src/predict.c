#include "common.h"

/* predict the value of a gaussian node without parents. */
SEXP gpred(SEXP fitted, SEXP data, SEXP debug) {

int i = 0, *ndata = INTEGER(data), *debuglevel = LOGICAL(debug);
double *mean = NULL, *res = NULL;
SEXP result;

  /* get the (only) coefficient of the linear regression. */
  mean = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, *ndata));
  res = REAL(result);

  /* copy the mean in the return value. */
  for (i = 0; i < *ndata; i++)
    res[i] = *mean;

  if (*debuglevel > 0) {

    Rprintf("  > prediction for all observations is %lf.\n", *mean);

  }/*THEN*/

  UNPROTECT(1);

  return result;

}/*GPRED*/

/* predict the value of a gaussian node with one or more parents. */
SEXP cgpred(SEXP fitted, SEXP data, SEXP debug)  {

int i = 0, j = 0, ndata = LENGTH(VECTOR_ELT(data, 0)), ncols = LENGTH(data);
int *debuglevel = LOGICAL(debug);
double *res = NULL, *coefs = NULL;
double **columns = NULL;
SEXP result;

  /* get the coefficient of the linear regression. */
  coefs = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, ndata));
  res = REAL(result);

  /* dereference the columns of the data frame. */
  columns = (double **) alloc1dpointer(ncols);
  for (i = 0; i < ncols; i++)
    columns[i] = REAL(VECTOR_ELT(data, i));

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    res[i] = coefs[0];

    for (j = 0; j < ncols; j++)
      res[i] += columns[j][i] * coefs[j + 1];

    if (*debuglevel > 0) {

      Rprintf("  > prediction for observation %d is %lf with predictor:\n",
        i + 1, res[i]);

      Rprintf("    (%lf) + (%lf) * (%lf)", coefs[0], columns[0][i], coefs[1]);
      for (j = 1; j < ncols; j++)
        Rprintf(" + (%lf) * (%lf)", columns[j][i], coefs[j + 1]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CGPRED*/

/* predict the value of a discrete node without parents. */
SEXP dpred(SEXP fitted, SEXP data, SEXP debug) {

int i = 0, nmax = 0, ndata = LENGTH(data), length = 0;
int *res = NULL, *debuglevel = LOGICAL(debug), *iscratch = NULL, *maxima = NULL;
double *prob = NULL, *dscratch = NULL, *buf = NULL;
SEXP ptab, result, tr_levels = getAttrib(data, R_LevelsSymbol);

  /* get the probabilities of the multinomial distribution. */
  ptab = getListElement(fitted, "prob");
  length = LENGTH(ptab);
  prob = REAL(ptab);

  /* create the vector of indexes. */
  iscratch = alloc1dcont(length);
  for (i = 0; i < length; i++)
    iscratch[i] = i + 1;

  /* create a scratch copy of the array. */
  buf = alloc1dreal(length);
  dscratch = alloc1dreal(length);
  memcpy(dscratch, prob, length * sizeof(double));

  /* allocate the array for the indexes of the maxima. */
  maxima = alloc1dcont(length);

  /* find out the mode(s). */
  all_max(dscratch, length, maxima, &nmax, iscratch, buf);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  if (nmax == 1) {

    /* copy the index of the mode in the return value. */
    for (i = 0; i < ndata; i++)
      res[i] = maxima[0];

    if (*debuglevel > 0) {

      if (res[0] == NA_INTEGER)
        Rprintf("  > prediction for all observations is NA with probabilities:\n");
      else
        Rprintf("  > prediction for all observations is '%s' with probabilities:\n",
          CHAR(STRING_ELT(tr_levels, res[0] - 1)));

      Rprintf("  ");
      for (i = 0; i < LENGTH(ptab); i++)
        Rprintf("  %lf", prob[i]);
      Rprintf("\n");

    }/*THEN*/

  }/*THEN*/
  else {

    /* break ties: sample with replacement from all the maxima. */
    GetRNGstate();
    SampleReplace(ndata, nmax, res, maxima);
    PutRNGstate();

    if (*debuglevel > 0) {

      Rprintf("  > there are %d levels tied for prediction, applying tie breaking.\n", nmax);
      Rprintf("  > tied levels are:");
      for (i = 0; i < nmax; i++)
        Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[i] - 1)));
      Rprintf(".\n");

    }/*THEN*/

  }/*ELSE*/

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*DPRED*/

/* predict the value of a discrete node with one or more parents. */
SEXP cdpred(SEXP fitted, SEXP data, SEXP parents, SEXP debug) {

int i = 0, k = 0, ndata = LENGTH(data), nrows = 0, ncols = 0;
int *configs = INTEGER(parents), *debuglevel = LOGICAL(debug);
int *iscratch = NULL, *maxima = NULL, *nmax = NULL, *res = NULL;
double *prob = NULL, *dscratch = NULL, *buf = NULL;
SEXP temp, result, tr_levels = getAttrib(data, R_LevelsSymbol);

  /* get the probabilities of the multinomial distribution. */
  temp = getListElement(fitted, "prob");
  nrows = INT(getAttrib(temp, R_DimSymbol));
  ncols = LENGTH(temp) / nrows;
  prob = REAL(temp);

  /* create the vector of indexes. */
  iscratch = alloc1dcont(nrows);

  /* create a scratch copy of the array. */
  buf = alloc1dreal(nrows);
  dscratch = alloc1dreal(nrows * ncols);
  memcpy(dscratch, prob, nrows * ncols * sizeof(double));

  /* allocate the array for the indexes of the maxima. */
  maxima = alloc1dcont(nrows * ncols);

  /* allocate the maxima counters. */
  nmax = alloc1dcont(ncols);

  /* get the mode for each configuration. */
  for (i = 0; i < ncols; i++) {

    /* initialize the vector of indexes. */
    for (k = 0; k < nrows; k++)
      iscratch[k] = k + 1;

    /* find out the mode(s). */
    all_max(dscratch + i * nrows, nrows, maxima + i * nrows,
      nmax + i, iscratch, buf);

  }/*FOR*/

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  /* copy the index of the mode in the return value. */
  for (i = 0; i < ndata; i++) {

    if (nmax[configs[i] - 1] == 1) {

      res[i] = maxima[CMC(0, configs[i] - 1, nrows)];

      if (*debuglevel > 0) {

        if (res[i] == NA_INTEGER)
          Rprintf("  > prediction for observation %d is NA with probabilities:\n");
        else
          Rprintf("  > prediction for observation %d is '%s' with probabilities:\n",
            i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

        Rprintf("  ");
        for (int k = 0; k < nrows; k++)
          Rprintf("  %lf", (prob + nrows * (configs[i] - 1))[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(1, nmax[configs[i] - 1], res + i, maxima + (configs[i] - 1) * nrows);

      if (*debuglevel > 0) {

        Rprintf("  > there are %d levels tied for prediction of observation %d, applying tie breaking.\n",
          nmax[configs[i] - 1], i + 1);
        Rprintf("  > tied levels are:");
        for (k = 0; k < nmax[configs[i] - 1]; k++)
          Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[CMC(k, configs[i] - 1, nrows)] - 1)));
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  /* save the state of the random number generator. */
  PutRNGstate();

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*CDPRED*/

/* predict the value of the training variable in a naive Bayes or Tree-Augmented
 * naive Bayes classifier. */
SEXP naivepred(SEXP fitted, SEXP data, SEXP parents, SEXP training, SEXP prior,
    SEXP prob, SEXP debug) {

int i = 0, j = 0, k = 0, n = 0, nvars = LENGTH(fitted), nmax = 0, tr_nlevels = 0;
int *res = NULL, **ex = NULL, *ex_nlevels = NULL;
int idx = 0, *tr_id = INTEGER(training);
int *iscratch = NULL, *maxima = NULL, *prn = NULL, *debuglevel = LOGICAL(debug);
int *include_prob = LOGICAL(prob);
double **cpt = NULL, *pr = NULL, *scratch = NULL, *buf = NULL, *pt = NULL;
double sum = 0;
SEXP class, temp, tr, tr_levels, result, nodes, probtab, dimnames;

  /* cache the node labels. */
  nodes = getAttrib(fitted, R_NamesSymbol);

  /* cache the pointers to all the variables. */
  ex = (int **) alloc1dpointer(nvars);
  ex_nlevels = alloc1dcont(nvars);

  for (i = 0; i < nvars; i++) {

    temp = VECTOR_ELT(data, i);
    ex[i] = INTEGER(temp);
    ex_nlevels[i] = NLEVELS(temp);

  }/*FOR*/

  /* get the training variable and its levels. */
  n = LENGTH(VECTOR_ELT(data, 0));
  tr = getListElement(VECTOR_ELT(fitted, *tr_id - 1), "prob");
  tr_levels = VECTOR_ELT(getAttrib(tr, R_DimNamesSymbol), 0);
  tr_nlevels = LENGTH(tr_levels);
  /* get the prior distribution. */
  pr = REAL(prior);

  if (*debuglevel > 0) {

    Rprintf("* the prior distribution for the target variable is:\n");
    PrintValue(prior);

  }/*THEN*/

  /* allocate the scratch space used to compute posterior probabilities. */
  scratch = alloc1dreal(tr_nlevels);
  buf = alloc1dreal(tr_nlevels);

  /* cache the pointers to the conditional probability tables. */
  cpt = (double **) alloc1dpointer(nvars);

  for (i = 0; i < nvars; i++) 
    cpt[i] = REAL(getListElement(VECTOR_ELT(fitted, i), "prob"));

  /* dereference the parents' vector. */
  prn = INTEGER(parents);

  /* create the vector of indexes. */
  iscratch = alloc1dcont(tr_nlevels);

  /* allocate the array for the indexes of the maxima. */
  maxima = alloc1dcont(tr_nlevels);

  /* allocate the return value. */
  PROTECT(result = allocVector(INTSXP, n));
  res = INTEGER(result);

  /* allocate and initialize the table of the posterior probabilities. */
  if (*include_prob > 0) {

    PROTECT(probtab = allocMatrix(REALSXP, tr_nlevels, n));
    pt = REAL(probtab);
    memset(pt, '\0', n * tr_nlevels * sizeof(double));

  }/*THEN*/

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  /* for each observation... */
  for (i = 0; i < n; i++) {

    /* ... reset the scratch space and the indexes array... */
    for (k = 0; k < tr_nlevels; k++) {

      scratch[k] = log(pr[k]);
      iscratch[k] = k + 1;

    }/*FOR*/

    if (*debuglevel > 0)
      Rprintf("* predicting the value of observation %d.\n", i + 1);

    /* ... and for each conditional probability table... */
    for (j = 0; j < nvars; j++) {

      /* ... skip the training variable... */
      if (*tr_id == j + 1)
        continue;

      /* ... (this is the root node of the Chow-Liu tree) ... */
      if (prn[j] == NA_INTEGER) {

        /* ... and for each row of the conditional probability table... */
        for (k = 0; k < tr_nlevels; k++) {

          if (*debuglevel > 0) {

            Rprintf("  > node %s: picking cell %d (%d, %d) from the CPT (p = %lf).\n",
              NODE(j), CMC(ex[j][i] - 1, k, ex_nlevels[j]), ex[j][i], k + 1,
              cpt[j][CMC(ex[j][i] - 1, k, ex_nlevels[j])]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] += log(cpt[j][CMC(ex[j][i] - 1, k, ex_nlevels[j])]);

        }/*FOR*/

      }/*THEN*/
      else {

        /* ... and for each row of the conditional probability table... */
        for (k = 0; k < tr_nlevels; k++) {

          /* (the first dimension corresponds to the current node [X], the second
           * to the training node [Y], the third to the only parent of the current
           * node [Z]; CMC coordinates are computed as X + Y * NX + Z * NX * NY. */
          idx = (ex[j][i] - 1) + k * ex_nlevels[j] + 
                  (ex[prn[j] - 1][i] - 1) * ex_nlevels[j] * tr_nlevels;

          if (*debuglevel > 0) {

            Rprintf("  > node %s: picking cell %d (%d, %d, %d) from the CPT (p = %lf).\n",
              NODE(j), idx, ex[j][i], k + 1, ex[prn[j] - 1][i], cpt[j][idx]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] += log(cpt[j][idx]);

        }/*FOR*/

      }/*ELSE*/

    }/*FOR*/

    /* find out the mode(s). */
    all_max(scratch, tr_nlevels, maxima, &nmax, iscratch, buf);

    /* compute the posterior probabilities on the right scale, to attach them
     * to the return value. */
    if (*include_prob) {

      /* copy the log-probabilities from scratch. */
      memcpy(pt + i * tr_nlevels, scratch, tr_nlevels * sizeof(double));

      /* transform log-probabilitiees into plain probabilities. */
      for (k = 0, sum = 0; k < tr_nlevels; k++)
        sum += pt[i * tr_nlevels + k] = exp(pt[i * tr_nlevels + k] - scratch[maxima[0] - 1]);

      /* rescale them to sum up to 1. */
      for (k = 0; k < tr_nlevels; k++)
        pt[i * tr_nlevels + k] /= sum;

    }/*THEN*/

    if (nmax == 1) {

      res[i] = maxima[0];

      if (*debuglevel > 0) {

        Rprintf("  @ prediction for observation %d is '%s' with (log-)posterior:\n",
          i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

        Rprintf("  ");
        for (k = 0; k < tr_nlevels; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(1, nmax, res + i, maxima);

      if (*debuglevel > 0) {

        Rprintf("  @ there are %d levels tied for prediction of observation %d, applying tie breaking.\n", nmax, i + 1);

        Rprintf("  ");
        for (k = 0; k < tr_nlevels; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

        Rprintf("  @ tied levels are:");
        for (k = 0; k < nmax; k++)
          Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[k] - 1)));
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  /* save the state of the random number generator. */
  PutRNGstate();

  /* add back the attributes and the class to the return value. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, class);

  if (*include_prob > 0) {

    /* set the levels of the taregt variable as rownames. */
    PROTECT(dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, tr_levels);
    setAttrib(probtab, R_DimNamesSymbol, dimnames);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, install("prob"), probtab);

    UNPROTECT(4);

  }/*THEN*/
  else {

    UNPROTECT(2);

  }/*ELSE*/

  return result;

}/*NAIVEPRED*/

