
// monotone spline model, with unconstrained lambda
//  - can handle empty time periods

data {
  int p;                // total number of athletes
  int tmax;             // max time period
  int nT[tmax];         // number of athletes per time period
  int nT_unique[tmax];  // number of unique athletes per time period

  int sumnT;            // TODO: remove this, this is just sum(nt)
  int sumnT_unique;     // TODO: remove this, this is just sum(nT_unique)

  int B;                 // number of non-intercept monotone spline basis functions
  real sumB;             // initial range of data (not used)
  vector[B+1] init_lam;  // initial lambda parameters, including intercept
  matrix[sumnT, B] i_basis;
  matrix[sumnT, B] m_basis;
  real alpha;            // dispersion of lambda parameters

  // see https://mc-stan.org/docs/2_25/stan-users-guide/ragged-data-structs-section.html
  //  for info on working with "ragged" data structures
  vector[sumnT] y;          // all scores, in a single vector
  matrix[sumnT, p] Xbar;    // all model matrices, in a single matrix
  int ptcps[sumnT_unique];  // all participants, in a single vector
}
transformed data {
  // initial model values, always fixed
  real a0;
  // real b0;
  vector[p] m0;
  real V0_init;
  vector[p] V0;
  real lam0;

  // a0 = 0.01;
  // b0 = 0.01;
  a0 = 2.0;
  m0 = rep_vector(0, p);
  V0_init = 5.0;
  V0 = rep_vector(V0_init, p);
  lam0 = init_lam[1];
}

parameters {
  real<lower=0> w;
  vector<lower=0>[B] lam;
}

model {
  // observed outcomes
  vector[sumnT] ty;
  // useful position trackers
  int pos;
  int pos2;

  // model values
  vector[p] m[tmax];
  vector[p] V[tmax];   // only store diagonal of V matrices
  vector[tmax] a;
  vector[tmax] b;

  // priors on w, lam
  for (i in 1:B) {
    lam[i] ~ normal(init_lam[i+1], alpha);
  }
  w ~ normal(0, 1);
  // w ~ inv_gamma(2, 1);

  // transform data
  ty = lam0 + i_basis * lam;

  // priors on m, V, a, b
  m[1] = m0;
  V[1] = V0;
  a[1] = a0;
  //b[1] = b0;
  b[1] = variance(ty);

  // compute KF updates to get m, V, a, b for all time periods
  pos = 1;
  pos2 = 1;
  for (t in 1:tmax) {
    // get important values for given time period
    int n_t = nT[t];                                   // number of observations
    vector[n_t] ty_t = segment(ty, pos, n_t);          // observed scores
    matrix[n_t,p] X_t = block(Xbar, pos, 1, n_t, p);   // X matrix

    // get columns of participating athletes
    int psm = nT_unique[t];                            // number of athletes
    int ptcps_t[psm] = segment(ptcps, pos2, psm);      // participant indices

    // declare small posterior model values for given time period
    matrix[n_t,psm] Xsm_t = X_t[,ptcps_t];
    matrix[psm,psm] Vsm_t;
    matrix[psm,psm] Vinvsm_t;
    vector[psm] msm_t;

    // set up useful representation of V_t, t-dist scale matrix
    vector[p] Vdiag_t;
    matrix[n_t,n_t] scale_mat;

    // it's simpler to just store these instead of alwys typing out
    matrix[psm,psm] Vinvsm_pr = diag_matrix(1.0 ./ V[t][ptcps_t]);
    vector[psm] msm_pr = m[t][ptcps_t];

    // print("t: ", t)
    // update position trackers
    pos = pos + n_t;
    pos2 = pos2 + psm;

    // if time period is empty, simply move to next time period
    if (n_t == 0) {
      // add innovation variance but cap variances at V0_init
      Vdiag_t = V[t] + rep_vector(w, p);
      for(x in 1:p) {
        Vdiag_t[x] = fmin(Vdiag_t[x], V0_init);
      }

      V[t+1] = Vdiag_t;
      m[t+1] = m[t];
      a[t+1] = a[t];
      b[t+1] = b[t];
      continue;
    }

    // update model values
    Vinvsm_t = Vinvsm_pr + crossprod(Xsm_t);
    msm_t = mdivide_left_spd(
      Vinvsm_t,
      Vinvsm_pr * msm_pr + Xsm_t' * ty_t);

    // add posterior density! needs obs from time t, model from time t-1
    scale_mat = b[t]/a[t] *
      (diag_matrix(rep_vector(1.0, n_t)) +
        // quad_form(diag_matrix(V[t][ptcps_t]), Xsm_t'));  // slow for some reason
        Xsm_t * diag_matrix(V[t][ptcps_t]) * Xsm_t');
    target += multi_student_t_lpdf(
      ty_t |
      2*a[t],
      Xsm_t * msm_pr,
      scale_mat);

    // slow; compute full inverse of this non-diagonal matrix
    Vsm_t = inverse(Vinvsm_t);

    // store posterior values for this time period
    V[t][ptcps_t] = diagonal(Vsm_t);
    m[t][ptcps_t] = msm_t;
    a[t] = a[t] + n_t/2.0;
    b[t] = b[t] +
      (quad_form(Vinvsm_pr, msm_pr)
        + dot_self(ty_t)
        - quad_form(Vinvsm_t, msm_t)) / 2.0;

    if (t < tmax) {
      // cap variances at V0_init
      Vdiag_t = V[t] + rep_vector(w, p);
      for(x in 1:p) {
        Vdiag_t[x] = fmin(Vdiag_t[x], V0_init);
      }
      V[t+1] = Vdiag_t;
      m[t+1] = m[t];
      a[t+1] = a[t];
      b[t+1] = b[t];
    }
  }

  // add jacobian
  // NOTE: why don't we run into -inf issues???
  //  - Basis is a vector of 0s for the minimum element,
  //    so log(.) should be -Inf...
  target += sum(log(m_basis * lam));
}

