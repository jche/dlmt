// Generated by rstantools.  Do not edit by hand.

/*
    dlmt is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dlmt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dlmt.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_model_unc_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_model_unc");
    reader.add_event(176, 174, "end", "model_model_unc");
    return reader;
}
#include <stan_meta_header.hpp>
class model_model_unc
  : public stan::model::model_base_crtp<model_model_unc> {
private:
        int p;
        int tmax;
        std::vector<int> nT;
        std::vector<int> nT_unique;
        int sumnT;
        int sumnT_unique;
        int B;
        double sumB;
        vector_d init_lam;
        matrix_d i_basis;
        matrix_d m_basis;
        double alpha;
        vector_d y;
        matrix_d Xbar;
        std::vector<int> ptcps;
        double a0;
        vector_d m0;
        double V0_init;
        vector_d V0;
        double lam0;
public:
    model_model_unc(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_model_unc(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_model_unc_namespace::model_model_unc";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "tmax", "int", context__.to_vec());
            tmax = int(0);
            vals_i__ = context__.vals_i("tmax");
            pos__ = 0;
            tmax = vals_i__[pos__++];
            current_statement_begin__ = 8;
            validate_non_negative_index("nT", "tmax", tmax);
            context__.validate_dims("data initialization", "nT", "int", context__.to_vec(tmax));
            nT = std::vector<int>(tmax, int(0));
            vals_i__ = context__.vals_i("nT");
            pos__ = 0;
            size_t nT_k_0_max__ = tmax;
            for (size_t k_0__ = 0; k_0__ < nT_k_0_max__; ++k_0__) {
                nT[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("nT_unique", "tmax", tmax);
            context__.validate_dims("data initialization", "nT_unique", "int", context__.to_vec(tmax));
            nT_unique = std::vector<int>(tmax, int(0));
            vals_i__ = context__.vals_i("nT_unique");
            pos__ = 0;
            size_t nT_unique_k_0_max__ = tmax;
            for (size_t k_0__ = 0; k_0__ < nT_unique_k_0_max__; ++k_0__) {
                nT_unique[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "sumnT", "int", context__.to_vec());
            sumnT = int(0);
            vals_i__ = context__.vals_i("sumnT");
            pos__ = 0;
            sumnT = vals_i__[pos__++];
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "sumnT_unique", "int", context__.to_vec());
            sumnT_unique = int(0);
            vals_i__ = context__.vals_i("sumnT_unique");
            pos__ = 0;
            sumnT_unique = vals_i__[pos__++];
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "B", "int", context__.to_vec());
            B = int(0);
            vals_i__ = context__.vals_i("B");
            pos__ = 0;
            B = vals_i__[pos__++];
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "sumB", "double", context__.to_vec());
            sumB = double(0);
            vals_r__ = context__.vals_r("sumB");
            pos__ = 0;
            sumB = vals_r__[pos__++];
            current_statement_begin__ = 16;
            validate_non_negative_index("init_lam", "(B + 1)", (B + 1));
            context__.validate_dims("data initialization", "init_lam", "vector_d", context__.to_vec((B + 1)));
            init_lam = Eigen::Matrix<double, Eigen::Dynamic, 1>((B + 1));
            vals_r__ = context__.vals_r("init_lam");
            pos__ = 0;
            size_t init_lam_j_1_max__ = (B + 1);
            for (size_t j_1__ = 0; j_1__ < init_lam_j_1_max__; ++j_1__) {
                init_lam(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 17;
            validate_non_negative_index("i_basis", "sumnT", sumnT);
            validate_non_negative_index("i_basis", "B", B);
            context__.validate_dims("data initialization", "i_basis", "matrix_d", context__.to_vec(sumnT,B));
            i_basis = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(sumnT, B);
            vals_r__ = context__.vals_r("i_basis");
            pos__ = 0;
            size_t i_basis_j_2_max__ = B;
            size_t i_basis_j_1_max__ = sumnT;
            for (size_t j_2__ = 0; j_2__ < i_basis_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < i_basis_j_1_max__; ++j_1__) {
                    i_basis(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 18;
            validate_non_negative_index("m_basis", "sumnT", sumnT);
            validate_non_negative_index("m_basis", "B", B);
            context__.validate_dims("data initialization", "m_basis", "matrix_d", context__.to_vec(sumnT,B));
            m_basis = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(sumnT, B);
            vals_r__ = context__.vals_r("m_basis");
            pos__ = 0;
            size_t m_basis_j_2_max__ = B;
            size_t m_basis_j_1_max__ = sumnT;
            for (size_t j_2__ = 0; j_2__ < m_basis_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < m_basis_j_1_max__; ++j_1__) {
                    m_basis(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 19;
            context__.validate_dims("data initialization", "alpha", "double", context__.to_vec());
            alpha = double(0);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            alpha = vals_r__[pos__++];
            current_statement_begin__ = 23;
            validate_non_negative_index("y", "sumnT", sumnT);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(sumnT));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(sumnT);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = sumnT;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 24;
            validate_non_negative_index("Xbar", "sumnT", sumnT);
            validate_non_negative_index("Xbar", "p", p);
            context__.validate_dims("data initialization", "Xbar", "matrix_d", context__.to_vec(sumnT,p));
            Xbar = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(sumnT, p);
            vals_r__ = context__.vals_r("Xbar");
            pos__ = 0;
            size_t Xbar_j_2_max__ = p;
            size_t Xbar_j_1_max__ = sumnT;
            for (size_t j_2__ = 0; j_2__ < Xbar_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Xbar_j_1_max__; ++j_1__) {
                    Xbar(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 25;
            validate_non_negative_index("ptcps", "sumnT_unique", sumnT_unique);
            context__.validate_dims("data initialization", "ptcps", "int", context__.to_vec(sumnT_unique));
            ptcps = std::vector<int>(sumnT_unique, int(0));
            vals_i__ = context__.vals_i("ptcps");
            pos__ = 0;
            size_t ptcps_k_0_max__ = sumnT_unique;
            for (size_t k_0__ = 0; k_0__ < ptcps_k_0_max__; ++k_0__) {
                ptcps[k_0__] = vals_i__[pos__++];
            }
            // initialize transformed data variables
            current_statement_begin__ = 29;
            a0 = double(0);
            stan::math::fill(a0, DUMMY_VAR__);
            current_statement_begin__ = 31;
            validate_non_negative_index("m0", "p", p);
            m0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            stan::math::fill(m0, DUMMY_VAR__);
            current_statement_begin__ = 32;
            V0_init = double(0);
            stan::math::fill(V0_init, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("V0", "p", p);
            V0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            stan::math::fill(V0, DUMMY_VAR__);
            current_statement_begin__ = 34;
            lam0 = double(0);
            stan::math::fill(lam0, DUMMY_VAR__);
            // execute transformed data statements
            current_statement_begin__ = 38;
            stan::math::assign(a0, 2.0);
            current_statement_begin__ = 39;
            stan::math::assign(m0, rep_vector(0, p));
            current_statement_begin__ = 40;
            stan::math::assign(V0_init, 5.0);
            current_statement_begin__ = 41;
            stan::math::assign(V0, rep_vector(V0_init, p));
            current_statement_begin__ = 42;
            stan::math::assign(lam0, get_base1(init_lam, 1, "init_lam", 1));
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 46;
            num_params_r__ += 1;
            current_statement_begin__ = 47;
            validate_non_negative_index("lam", "B", B);
            num_params_r__ += B;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_model_unc() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 46;
        if (!(context__.contains_r("w")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable w missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("w");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "w", "double", context__.to_vec());
        double w(0);
        w = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, w);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable w: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 47;
        if (!(context__.contains_r("lam")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lam missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lam");
        pos__ = 0U;
        validate_non_negative_index("lam", "B", B);
        context__.validate_dims("parameter initialization", "lam", "vector_d", context__.to_vec(B));
        Eigen::Matrix<double, Eigen::Dynamic, 1> lam(B);
        size_t lam_j_1_max__ = B;
        for (size_t j_1__ = 0; j_1__ < lam_j_1_max__; ++j_1__) {
            lam(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, lam);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lam: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 46;
            local_scalar_t__ w;
            (void) w;  // dummy to suppress unused var warning
            if (jacobian__)
                w = in__.scalar_lb_constrain(0, lp__);
            else
                w = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 47;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lam;
            (void) lam;  // dummy to suppress unused var warning
            if (jacobian__)
                lam = in__.vector_lb_constrain(0, B, lp__);
            else
                lam = in__.vector_lb_constrain(0, B);
            // model body
            {
            current_statement_begin__ = 52;
            validate_non_negative_index("ty", "sumnT", sumnT);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ty(sumnT);
            stan::math::initialize(ty, DUMMY_VAR__);
            stan::math::fill(ty, DUMMY_VAR__);
            current_statement_begin__ = 54;
            int pos(0);
            (void) pos;  // dummy to suppress unused var warning
            stan::math::fill(pos, std::numeric_limits<int>::min());
            current_statement_begin__ = 55;
            int pos2(0);
            (void) pos2;  // dummy to suppress unused var warning
            stan::math::fill(pos2, std::numeric_limits<int>::min());
            current_statement_begin__ = 58;
            validate_non_negative_index("m", "p", p);
            validate_non_negative_index("m", "tmax", tmax);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>  > m(tmax, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(p));
            stan::math::initialize(m, DUMMY_VAR__);
            stan::math::fill(m, DUMMY_VAR__);
            current_statement_begin__ = 59;
            validate_non_negative_index("V", "p", p);
            validate_non_negative_index("V", "tmax", tmax);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>  > V(tmax, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(p));
            stan::math::initialize(V, DUMMY_VAR__);
            stan::math::fill(V, DUMMY_VAR__);
            current_statement_begin__ = 60;
            validate_non_negative_index("a", "tmax", tmax);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a(tmax);
            stan::math::initialize(a, DUMMY_VAR__);
            stan::math::fill(a, DUMMY_VAR__);
            current_statement_begin__ = 61;
            validate_non_negative_index("b", "tmax", tmax);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> b(tmax);
            stan::math::initialize(b, DUMMY_VAR__);
            stan::math::fill(b, DUMMY_VAR__);
            current_statement_begin__ = 64;
            for (int i = 1; i <= B; ++i) {
                current_statement_begin__ = 65;
                lp_accum__.add(normal_log<propto__>(get_base1(lam, i, "lam", 1), get_base1(init_lam, (i + 1), "init_lam", 1), alpha));
            }
            current_statement_begin__ = 67;
            lp_accum__.add(normal_log<propto__>(w, 0, 1));
            current_statement_begin__ = 71;
            stan::math::assign(ty, add(lam0, multiply(i_basis, lam)));
            current_statement_begin__ = 74;
            stan::model::assign(m, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        m0, 
                        "assigning variable m");
            current_statement_begin__ = 75;
            stan::model::assign(V, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        V0, 
                        "assigning variable V");
            current_statement_begin__ = 76;
            stan::model::assign(a, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        a0, 
                        "assigning variable a");
            current_statement_begin__ = 78;
            stan::model::assign(b, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        variance(ty), 
                        "assigning variable b");
            current_statement_begin__ = 81;
            stan::math::assign(pos, 1);
            current_statement_begin__ = 82;
            stan::math::assign(pos2, 1);
            current_statement_begin__ = 83;
            for (int t = 1; t <= tmax; ++t) {
                {
                current_statement_begin__ = 85;
                int n_t(0);
                (void) n_t;  // dummy to suppress unused var warning
                stan::math::fill(n_t, std::numeric_limits<int>::min());
                stan::math::assign(n_t,get_base1(nT, t, "nT", 1));
                current_statement_begin__ = 86;
                validate_non_negative_index("ty_t", "n_t", n_t);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ty_t(n_t);
                stan::math::initialize(ty_t, DUMMY_VAR__);
                stan::math::fill(ty_t, DUMMY_VAR__);
                stan::math::assign(ty_t,segment(ty, pos, n_t));
                current_statement_begin__ = 87;
                validate_non_negative_index("X_t", "n_t", n_t);
                validate_non_negative_index("X_t", "p", p);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> X_t(n_t, p);
                stan::math::initialize(X_t, DUMMY_VAR__);
                stan::math::fill(X_t, DUMMY_VAR__);
                stan::math::assign(X_t,block(Xbar, pos, 1, n_t, p));
                current_statement_begin__ = 90;
                int psm(0);
                (void) psm;  // dummy to suppress unused var warning
                stan::math::fill(psm, std::numeric_limits<int>::min());
                stan::math::assign(psm,get_base1(nT_unique, t, "nT_unique", 1));
                current_statement_begin__ = 91;
                validate_non_negative_index("ptcps_t", "psm", psm);
                std::vector<int  > ptcps_t(psm, int(0));
                stan::math::fill(ptcps_t, std::numeric_limits<int>::min());
                stan::math::assign(ptcps_t,segment(ptcps, pos2, psm));
                current_statement_begin__ = 94;
                validate_non_negative_index("Xsm_t", "n_t", n_t);
                validate_non_negative_index("Xsm_t", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Xsm_t(n_t, psm);
                stan::math::initialize(Xsm_t, DUMMY_VAR__);
                stan::math::fill(Xsm_t, DUMMY_VAR__);
                stan::math::assign(Xsm_t,stan::model::rvalue(X_t, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list())), "X_t"));
                current_statement_begin__ = 95;
                validate_non_negative_index("Vsm_t", "psm", psm);
                validate_non_negative_index("Vsm_t", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Vsm_t(psm, psm);
                stan::math::initialize(Vsm_t, DUMMY_VAR__);
                stan::math::fill(Vsm_t, DUMMY_VAR__);
                current_statement_begin__ = 96;
                validate_non_negative_index("Vinvsm_t", "psm", psm);
                validate_non_negative_index("Vinvsm_t", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Vinvsm_t(psm, psm);
                stan::math::initialize(Vinvsm_t, DUMMY_VAR__);
                stan::math::fill(Vinvsm_t, DUMMY_VAR__);
                current_statement_begin__ = 97;
                validate_non_negative_index("msm_t", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> msm_t(psm);
                stan::math::initialize(msm_t, DUMMY_VAR__);
                stan::math::fill(msm_t, DUMMY_VAR__);
                current_statement_begin__ = 100;
                validate_non_negative_index("Vdiag_t", "p", p);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Vdiag_t(p);
                stan::math::initialize(Vdiag_t, DUMMY_VAR__);
                stan::math::fill(Vdiag_t, DUMMY_VAR__);
                current_statement_begin__ = 101;
                validate_non_negative_index("scale_mat", "n_t", n_t);
                validate_non_negative_index("scale_mat", "n_t", n_t);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> scale_mat(n_t, n_t);
                stan::math::initialize(scale_mat, DUMMY_VAR__);
                stan::math::fill(scale_mat, DUMMY_VAR__);
                current_statement_begin__ = 104;
                validate_non_negative_index("Vinvsm_pr", "psm", psm);
                validate_non_negative_index("Vinvsm_pr", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Vinvsm_pr(psm, psm);
                stan::math::initialize(Vinvsm_pr, DUMMY_VAR__);
                stan::math::fill(Vinvsm_pr, DUMMY_VAR__);
                stan::math::assign(Vinvsm_pr,diag_matrix(elt_divide(1.0, stan::model::rvalue(get_base1(V, t, "V", 1), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list()), "V[t]"))));
                current_statement_begin__ = 105;
                validate_non_negative_index("msm_pr", "psm", psm);
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> msm_pr(psm);
                stan::math::initialize(msm_pr, DUMMY_VAR__);
                stan::math::fill(msm_pr, DUMMY_VAR__);
                stan::math::assign(msm_pr,stan::model::rvalue(get_base1(m, t, "m", 1), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list()), "m[t]"));
                current_statement_begin__ = 109;
                stan::math::assign(pos, (pos + n_t));
                current_statement_begin__ = 110;
                stan::math::assign(pos2, (pos2 + psm));
                current_statement_begin__ = 113;
                if (as_bool(logical_eq(n_t, 0))) {
                    current_statement_begin__ = 115;
                    stan::math::assign(Vdiag_t, add(get_base1(V, t, "V", 1), rep_vector(w, p)));
                    current_statement_begin__ = 116;
                    for (int x = 1; x <= p; ++x) {
                        current_statement_begin__ = 117;
                        stan::model::assign(Vdiag_t, 
                                    stan::model::cons_list(stan::model::index_uni(x), stan::model::nil_index_list()), 
                                    stan::math::fmin(get_base1(Vdiag_t, x, "Vdiag_t", 1), V0_init), 
                                    "assigning variable Vdiag_t");
                    }
                    current_statement_begin__ = 120;
                    stan::model::assign(V, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                Vdiag_t, 
                                "assigning variable V");
                    current_statement_begin__ = 121;
                    stan::model::assign(m, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(m, t, "m", 1), 
                                "assigning variable m");
                    current_statement_begin__ = 122;
                    stan::model::assign(a, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(a, t, "a", 1), 
                                "assigning variable a");
                    current_statement_begin__ = 123;
                    stan::model::assign(b, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(b, t, "b", 1), 
                                "assigning variable b");
                    current_statement_begin__ = 124;
                    continue;
                }
                current_statement_begin__ = 128;
                stan::math::assign(Vinvsm_t, add(Vinvsm_pr, crossprod(Xsm_t)));
                current_statement_begin__ = 129;
                stan::math::assign(msm_t, mdivide_left_spd(Vinvsm_t, add(multiply(Vinvsm_pr, msm_pr), multiply(transpose(Xsm_t), ty_t))));
                current_statement_begin__ = 134;
                stan::math::assign(scale_mat, multiply((get_base1(b, t, "b", 1) / get_base1(a, t, "a", 1)), add(diag_matrix(rep_vector(1.0, n_t)), multiply(multiply(Xsm_t, diag_matrix(stan::model::rvalue(get_base1(V, t, "V", 1), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list()), "V[t]"))), transpose(Xsm_t)))));
                current_statement_begin__ = 138;
                lp_accum__.add(multi_student_t_log(ty_t, (2 * get_base1(a, t, "a", 1)), multiply(Xsm_t, msm_pr), scale_mat));
                current_statement_begin__ = 145;
                stan::math::assign(Vsm_t, inverse(Vinvsm_t));
                current_statement_begin__ = 148;
                stan::model::assign(V, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list())), 
                            diagonal(Vsm_t), 
                            "assigning variable V");
                current_statement_begin__ = 149;
                stan::model::assign(m, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::cons_list(stan::model::index_multi(ptcps_t), stan::model::nil_index_list())), 
                            msm_t, 
                            "assigning variable m");
                current_statement_begin__ = 150;
                stan::model::assign(a, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::nil_index_list()), 
                            (get_base1(a, t, "a", 1) + (n_t / 2.0)), 
                            "assigning variable a");
                current_statement_begin__ = 151;
                stan::model::assign(b, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::nil_index_list()), 
                            (get_base1(b, t, "b", 1) + (((quad_form(Vinvsm_pr, msm_pr) + dot_self(ty_t)) - quad_form(Vinvsm_t, msm_t)) / 2.0)), 
                            "assigning variable b");
                current_statement_begin__ = 156;
                if (as_bool(logical_lt(t, tmax))) {
                    current_statement_begin__ = 158;
                    stan::math::assign(Vdiag_t, add(get_base1(V, t, "V", 1), rep_vector(w, p)));
                    current_statement_begin__ = 159;
                    for (int x = 1; x <= p; ++x) {
                        current_statement_begin__ = 160;
                        stan::model::assign(Vdiag_t, 
                                    stan::model::cons_list(stan::model::index_uni(x), stan::model::nil_index_list()), 
                                    stan::math::fmin(get_base1(Vdiag_t, x, "Vdiag_t", 1), V0_init), 
                                    "assigning variable Vdiag_t");
                    }
                    current_statement_begin__ = 162;
                    stan::model::assign(V, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                Vdiag_t, 
                                "assigning variable V");
                    current_statement_begin__ = 163;
                    stan::model::assign(m, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(m, t, "m", 1), 
                                "assigning variable m");
                    current_statement_begin__ = 164;
                    stan::model::assign(a, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(a, t, "a", 1), 
                                "assigning variable a");
                    current_statement_begin__ = 165;
                    stan::model::assign(b, 
                                stan::model::cons_list(stan::model::index_uni((t + 1)), stan::model::nil_index_list()), 
                                get_base1(b, t, "b", 1), 
                                "assigning variable b");
                }
                }
            }
            current_statement_begin__ = 173;
            lp_accum__.add(sum(stan::math::log(multiply(m_basis, lam))));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("w");
        names__.push_back("lam");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(B);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_model_unc_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double w = in__.scalar_lb_constrain(0);
        vars__.push_back(w);
        Eigen::Matrix<double, Eigen::Dynamic, 1> lam = in__.vector_lb_constrain(0, B);
        size_t lam_j_1_max__ = B;
        for (size_t j_1__ = 0; j_1__ < lam_j_1_max__; ++j_1__) {
            vars__.push_back(lam(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_model_unc";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "w";
        param_names__.push_back(param_name_stream__.str());
        size_t lam_j_1_max__ = B;
        for (size_t j_1__ = 0; j_1__ < lam_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lam" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "w";
        param_names__.push_back(param_name_stream__.str());
        size_t lam_j_1_max__ = B;
        for (size_t j_1__ = 0; j_1__ < lam_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lam" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_model_unc_namespace::model_model_unc stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif