/*
File $Id: QuadProg++.cc 232 2007-06-21 12:29:00Z digasper $

 Author: Luca Di Gaspero
 DIEGM - University of Udine, Italy
 luca.digaspero@uniud.it
 http://www.diegm.uniud.it/digaspero/

 This software may be modified and distributed under the terms
 of the MIT license.  See the LICENSE file for details.

 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "QuadProg++.hh"



#ifdef QUADPROGPP_ENABLE_EIGEN

#include "EigenHelpers.hh"

#else

#include "ArrayHelpers.hh"

#endif



#ifdef QUADPROGPP_ENABLE_TRACING

#define QPPP_TRACE_MESSAGE(msg)                         std::cout << msg
#define QPPP_TRACE_MATRIX(name, matrix, rows, cols)     print_matrix(name, matrix, rows, cols)
#define QPPP_TRACE_VECTOR(name, matrix, rows)           print_vector(name, matrix, rows)

#else

#define QPPP_TRACE_MESSAGE(msg)
#define QPPP_TRACE_MATRIX(name, matrix, rows, cols)
#define QPPP_TRACE_VECTOR(name, matrix, rows)

#endif



namespace QuadProgpp
{

template<typename T>
void print_vector(const char* name, const std::vector<T> & v, int n = -1)
{
    std::ostringstream s;
    std::string t;
    if (n == -1)
        n = v.size();

    s << name << ": " << std::endl << " ";
    for (int i = 0; i < n; i++)
    {
        s << v[i] << ", ";
    }
    t = s.str();
    t = t.substr(0, t.size() - 2); // To remove the trailing space and comma

    std::cout << t << std::endl;
}


//-----------------------------------------------------------------------
// Utility function for computing the euclidean distance between two numbers
//-----------------------------------------------------------------------
inline double distance(double a, double b)
{
    register double a1, b1, t;
    a1 = std::fabs(a);
    b1 = std::fabs(b);
    if (a1 > b1)
    {
        t = (b1 / a1);
        return a1 * std::sqrt(1.0 + t * t);
    }
    else if (b1 > a1)
    {
        t = (a1 / b1);
        return b1 * std::sqrt(1.0 + t * t);
    }
    return a1 * std::sqrt(2.0);
}


class ConstraintStatus
{
    public:
        enum Status
        {
            INACTIVE = 0,
            ACTIVE = 1
        };
};


class Solver::Implementation
{
    public:
        int     iter;
        double  f_value;


    private:
        double  inf;
        double  epsilon;

        QPPP_MATRIX(double) R;
        QPPP_MATRIX(double) J;
        QPPP_VECTOR(double) ci_violations;
        QPPP_VECTOR(double) z;
        QPPP_VECTOR(double) r;
        QPPP_VECTOR(double) d;
        QPPP_VECTOR(double) np;
        QPPP_VECTOR(double) u;
        QPPP_VECTOR(double) x_old;
        QPPP_VECTOR(double) u_old;
        QPPP_VECTOR(int) A;
        QPPP_VECTOR(int) A_old;
        QPPP_VECTOR(bool) iaexcl;

        std::vector<ConstraintStatus::Status> ci_status;

        CholeskyDecomposition<double>   chol;


    private:
        bool add_constraint(QPPP_VECTOR(double)& d,
                            double& R_norm,
                            const int iq,
                            const int num_var)
        {
            QPPP_TRACE_MESSAGE("Add constraint " << iq << '/');

            register int i, j, k;
            double cc, ss, h, t1, t2, xny;

            /* we have to find the Givens rotation which will reduce the element
              d[j] to zero.
              if it is already zero we don't have to do anything, except of
              decreasing j */
            for (j = num_var - 1; j >= iq + 1; j--)
            {
                /* The Givens rotation is done with the matrix (cc cs, cs -cc).
                If cc is one, then element (j) of d is zero compared with element
                (j - 1). Hence we don't have to do anything.
                If cc is zero, then we just have to switch column (j) and column (j - 1)
                of J. Since we only switch columns in J, we have to be careful how we
                update d depending on the sign of gs.
                Otherwise we have to apply the Givens rotation to these columns.
                The i - 1 element of d has to be updated to h. */
                cc = d[j - 1];
                ss = d[j];
                h = distance(cc, ss);
                if (std::fabs(h) < epsilon) // h == 0
                    continue;
                d[j] = 0.0;
                ss = ss / h;
                cc = cc / h;
                if (cc < 0.0)
                {
                    cc = -cc;
                    ss = -ss;
                    d[j - 1] = -h;
                }
                else
                    d[j - 1] = h;
                xny = ss / (1.0 + cc);
                for (k = 0; k < num_var; k++)
                {
                    t1 = J(k,j - 1);
                    t2 = J(k, j);
                    J(k, j - 1) = t1 * cc + t2 * ss;
                    J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
                }
            }
            /* To update R we have to put the iq+1 components of the d vector
              into column iq of R
              */
            for (i = 0; i < iq+1; i++)
                R(i, iq) = d[i];

            QPPP_TRACE_MESSAGE(iq << std::endl);
            QPPP_TRACE_MATRIX("R", R, iq+1, iq+1);
            QPPP_TRACE_MATRIX("J", J, J.rows(), J.cols());
            QPPP_TRACE_VECTOR("d", d, d.size());

            if (std::fabs(d[iq]) <= epsilon * R_norm)
            {
                // problem degenerate
                return false;
            }
            R_norm = std::max<double>(R_norm, std::fabs(d[iq]));
            return true;
        }


        void delete_constraint( QPPP_VECTOR(double)& u,
                                int iq,
                                int l,
                                const int num_var,
                                const int num_ce)
        {
            QPPP_TRACE_MESSAGE("Delete constraint " << l << ' ' << iq);


            register int i, j, k, qq = -1; // just to prevent warnings from smart compilers
            double cc, ss, h, xny, t1, t2;

            /* Find the index qq for active constraint l to be removed */
            for (i = 0; i < iq - num_ce; i++)
                if (A[i] == l)
                {
                    qq = i;
                    break;
                }

            /* remove the constraint from the active set and the duals */
            for (i = qq; i < iq - 1; i++)
            {
                A[i-num_ce] = A[i + 1 - num_ce];
                u[i] = u[i + 1];
                for (j = 0; j < num_var; j++)
                    R(j, i) = R(j, i + 1);
            }

            A[iq - 1 - num_ce] = A[iq - num_ce];
            u[iq - 1] = u[iq];
            /* constraint has been fully removed */

            QPPP_TRACE_MESSAGE('/' << iq << std::endl);

            for (j = qq; j < iq-1; j++)
            {
                cc = R(j, j);
                ss = R(j + 1, j);
                h = distance(cc, ss);
                if (std::fabs(h) < epsilon) // h == 0
                    continue;
                cc = cc / h;
                ss = ss / h;
                R(j + 1, j) = 0.0;
                if (cc < 0.0)
                {
                    R(j, j) = -h;
                    cc = -cc;
                    ss = -ss;
                }
                else
                    R(j, j) = h;

                xny = ss / (1.0 + cc);
                for (k = j + 1; k < iq-1; k++)
                {
                    t1 = R(j, k);
                    t2 = R(j + 1, k);
                    R(j, k) = t1 * cc + t2 * ss;
                    R(j + 1,k) = xny * (t1 + R(j, k)) - t2;
                }
                for (k = 0; k < num_var; k++)
                {
                    t1 = J(k, j);
                    t2 = J(k,j + 1);
                    J(k, j) = t1 * cc + t2 * ss;
                    J(k,j + 1) = xny * (J(k, j) + t1) - t2;
                }
            }
        }


    public:
        Implementation()
        {
            if (std::numeric_limits<double>::has_infinity)
            {
                inf = std::numeric_limits<double>::infinity();
            }
            else
            {
                inf = 1.0E300;
            }

            epsilon = std::numeric_limits<double>::epsilon();
        }


        // The Solving function, implementing the Goldfarb-Idnani method
        Status::Value solve(QPPP_MATRIX(double)& G, QPPP_VECTOR(double)& g0,
                     const QPPP_MATRIX(double)& CE, const QPPP_VECTOR(double)& ce0,
                     const QPPP_MATRIX(double)& CI, const QPPP_VECTOR(double)& ci0,
                     QPPP_VECTOR(double)& x)
        {
            int num_var = G.cols();
            int num_ce = CE.cols();
            int num_ci = CI.cols();

            if (G.rows() != num_var)
            {
                std::ostringstream msg;
                msg << "The matrix G is not a squared matrix (" << G.rows() << " x " << G.cols() << ")";
                throw std::logic_error(msg.str());
            }
            if (CE.rows() != num_var)
            {
                std::ostringstream msg;
                msg << "The matrix CE is incompatible (incorrect number of rows " << CE.rows() << " , expecting " << num_var << ")";
                throw std::logic_error(msg.str());
            }
            if (ce0.size() != num_ce)
            {
                std::ostringstream msg;
                msg << "The vector ce0 is incompatible (incorrect dimension " << ce0.size() << ", expecting " << num_ce << ")";
                throw std::logic_error(msg.str());
            }
            if (CI.rows() != num_var)
            {
                std::ostringstream msg;
                msg << "The matrix CI is incompatible (incorrect number of rows " << CI.rows() << " , expecting " << num_var << ")";
                throw std::logic_error(msg.str());
            }
            if (ci0.size() != num_ci)
            {
                std::ostringstream msg;
                msg << "The vector ci0 is incompatible (incorrect dimension " << ci0.size() << ", expecting " << num_ci << ")";
                throw std::logic_error(msg.str());
            }
            x.resize(num_var);


            QPPP_TRACE_MESSAGE(std::endl << "Starting solve_quadprog" << std::endl);
            QPPP_TRACE_MATRIX("G", G, G.rows(), G.cols());
            QPPP_TRACE_VECTOR("g0", g0, g0.size());
            QPPP_TRACE_MATRIX("CE", CE, CE.rows(), CE.cols());
            QPPP_TRACE_VECTOR("ce0", ce0, ce0.size());
            QPPP_TRACE_MATRIX("CI", CI, CI.rows(), CI.cols());
            QPPP_TRACE_VECTOR("ci0", ci0, ci0.size());


            /*
             * Preprocessing phase
             */

            double trace_G, trace_J;
            /* compute the trace of the original matrix G */
            trace_G = G.trace();
            /* decompose the matrix G in the form L^T L */
            chol.compute(G);
            /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
            J.resize(num_var, num_var);
            chol.invert_upper(G,J);
            trace_J = J.trace();

            QPPP_TRACE_MATRIX("L", G, G.rows(), G.cols());
            QPPP_TRACE_MATRIX("J", J, J.rows(), J.cols());



            /* trace_G * trace_J is an estimate for cond(G) */

            /*
             * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
             * this is a feasible point in the dual space
             * x = G^-1 * g0
             */
            chol.solve(G, x, g0);
            for (int i = 0; i < num_var; ++i)
                x[i] = -x[i];
            /* and compute the current solution value */
            f_value = 0.5 * g0.dot(x);


            QPPP_TRACE_MESSAGE("Unconstrained solution: " << f_value << std::endl);
            QPPP_TRACE_VECTOR("x", x, x.size());


            /* Add equality constraints to the working set A */
            z.resize(num_var);
            r.resize(num_ce + num_ci);
            d.resize(num_var);
            np.resize(num_var);
            u.resize(num_ce + num_ci);
            R.resize(num_var, num_var);
            double R_norm = 1.0; /* this variable will hold the norm of the matrix R */

            for (int i = 0; i < num_ce; ++i)
            {
                for (int j = 0; j < num_var; ++j)
                    np[j] = CE(j, i);
                compute_d(d, J, np);
                update_z(z, J, d, i);
                update_r(R, r, d, i);


                QPPP_TRACE_MATRIX("R", R, num_var, i);
                QPPP_TRACE_VECTOR("z", z, num_var);
                QPPP_TRACE_VECTOR("r", r, i);
                QPPP_TRACE_VECTOR("d", d, num_var);


                /*
                 * compute full step length t2: i.e., the minimum step in
                 * primal space s.t. the contraint becomes feasible
                 */
                double t2 = 0.0;
                if (z.dot(z) > epsilon) // i.e. z != 0
                {
                    double z_dot_np = z.dot(np);

                    t2 = (-np.dot(x) - ce0[i]) / z_dot_np;

                    /* set x = x + t2 * z */
                    for (int k = 0; k < num_var; ++k)
                        x[k] += t2 * z[k];

                    /* set u = u+ */
                    for (int k = 0; k < i; ++k)
                        u[k] -= t2 * r[k];

                    /* compute the new solution value */
                    f_value += 0.5 * t2 * t2 * z_dot_np;
                }
                u[i] = t2;

                if (!add_constraint(d, R_norm, i, num_var))
                {
                    // Equality constraints are linearly dependent
                    throw std::runtime_error("Constraints are linearly dependent");
                    return (Status::FAILURE);
                }
            }


            /* set ci_status = K \ A */
            ci_status.assign(num_ci, ConstraintStatus::INACTIVE);

            A.resize(num_ci);
            A_old.resize(num_ci);
            x_old.resize(num_var);
            u_old.resize(num_ci + num_ce);

            ci_violations.resize(num_ci);
            iaexcl.resize(num_ci);
            iter = 0;
            /*
             * t is the step lenght, which is the minimum of the partial step
             * length t1 and the full step length t2
             */
            double t, t1, t2;

            int iq = num_ce;

        l1:
            ++iter;
            QPPP_TRACE_VECTOR("x", x, x.size());

            /* step 1: choose a violated constraint */
            for (int i = 0; i < iq-num_ce; ++i)
            {
                ci_status[ A[i] ] = ConstraintStatus::ACTIVE;
            }

            /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
            double ss = 0.0;
            int ip = 0; /* ip will be the index of the chosen violated constraint */


            // ci_violations = CI^T*x + ci0
            multiply_and_add(ci_violations,CI,x,ci0);


            QPPP_TRACE_VECTOR("s", ci_violations, num_ci);


            double psi = 0.0; /* this value will contain the sum of all infeasibilities */
            for (int i = 0; i < num_ci; ++i)
            {
                iaexcl[i] = true;
                psi += std::min(0.0, ci_violations[i]);
            }

            if (std::fabs(psi) <= num_ci * epsilon * trace_G * trace_J * 100.0)
            {
                /* numerically there are not infeasibilities anymore */
                return (Status::OK);
            }

            /* save old values for u and A */
            for (int i = 0; i < iq; ++i)
            {
                u_old[i] = u[i];
            }
            for (int i = 0; i < iq-num_ce; ++i)
            {
                A_old[i] = A[i];
            }
            /* and for x */
            x_old = x;

        l2: /* Step 2: check for feasibility and determine a new S-pair */
            for (int i = 0; i < num_ci; ++i)
            {
                if ((ci_violations[i] < ss) && (ci_status[i] == ConstraintStatus::INACTIVE) && (iaexcl[i]))
                {
                    ss = ci_violations[i];
                    ip = i;
                }
            }
            if (ss >= 0.0)
            {
                return (Status::OK);
            }

            /* set np = n[ip] */
            for (int i = 0; i < num_var; ++i)
                np[i] = CI(i, ip);
            /* set u = [u 0]^T */
            u[iq] = 0.0;
            /* add ip to the active set A */
            A[iq-num_ce] = ip;


            QPPP_TRACE_MESSAGE("Trying with constraint " << ip << std::endl);
            QPPP_TRACE_VECTOR("np", np, np.size());


        l2a:/* Step 2a: determine step direction */
            /* compute z = H np: the step direction in the primal space (through J, see the paper) */
            compute_d(d, J, np);
            update_z(z, J, d, iq);
            /* compute N* np (if q > 0): the negative of the step direction in the dual space */
            update_r(R, r, d, iq);

            double z_dot_np = z.dot(np);

            QPPP_TRACE_MESSAGE("Step direction z" << std::endl);
            QPPP_TRACE_VECTOR("z", z, z.size());
            QPPP_TRACE_VECTOR("r", r, iq + 1);
            QPPP_TRACE_VECTOR("u", u, iq + 1);
            QPPP_TRACE_VECTOR("d", d, d.size());
            QPPP_TRACE_VECTOR("A", A, iq + 1);


            /* Step 2b: compute step length */
            /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
            t1 = inf; /* +inf */
            /* find the index l s.t. it reaches the minimum of u+[x] / r */
            int l = 0;
            for (int k = num_ce; k < iq; ++k)
            {
                if (r[k] > 0.0)
                {
                    if (u[k] / r[k] < t1)
                    {
                        t1 = u[k] / r[k];
                        l = A[k-num_ce];
                    }
                }
            }
            /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
            if (z.dot(z) > epsilon) // i.e. z != 0
            {
                t2 = -ci_violations[ip] / z_dot_np;
                if (t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
                    t2 = inf;
            }
            else
                t2 = inf; /* +inf */

            /* the step is chosen as the minimum of t1 and t2 */
            t = std::min(t1, t2);


            QPPP_TRACE_MESSAGE("Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ");


            /* Step 2c: determine new S-pair and take step: */

            /* case (i): no step in primal or dual space */
            if (t >= inf)
            {
                /* QPP is infeasible */
                // FIXME: unbounded to raise
                return (Status::FAILURE);
            }
            /* case (ii): step in dual space */
            if (t2 >= inf)
            {
                /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
                for (int k = 0; k < iq; ++k)
                    u[k] -= t * r[k];
                u[iq] += t;
                ci_status[l] = ConstraintStatus::INACTIVE;
                delete_constraint(u, iq, l, num_var, num_ce);
                --iq;

                QPPP_TRACE_MESSAGE(" in dual space: " << f_value << std::endl);
                QPPP_TRACE_VECTOR("x", x, x.size());
                QPPP_TRACE_VECTOR("A", A, iq + 1);

                goto l2a;
            }

            /* case (iii): step in primal and dual space */

            /* set x = x + t * z */
            for (int k = 0; k < num_var; ++k)
                x[k] += t * z[k];
            /* update the solution value */
            f_value += t * z_dot_np * (0.5 * t + u[iq]);
            /* u = u + t * [-r 1] */
            for (int k = 0; k < iq; ++k)
                u[k] -= t * r[k];
            u[iq] += t;

            QPPP_TRACE_MESSAGE(" in both spaces: " << f_value << std::endl);
            QPPP_TRACE_VECTOR("x", x, x.size());
            QPPP_TRACE_VECTOR("u", u, iq + 1);
            QPPP_TRACE_VECTOR("r", r, iq + 1);
            QPPP_TRACE_VECTOR("A", A, iq + 1);


            if (std::fabs(t - t2) < epsilon)
            {
                QPPP_TRACE_MESSAGE("Full step has taken " << t << std::endl);
                QPPP_TRACE_VECTOR("x", x, x.size());

                /* full step has taken */
                /* add constraint ip to the active set*/
                if (add_constraint(d, R_norm, iq, num_var))
                {
                    ++iq;
                    ci_status[ip] = ConstraintStatus::ACTIVE;
                }
                else
                {
                    iaexcl[ip] = false;
                    delete_constraint(u, iq, ip, num_var, num_ce);

                    QPPP_TRACE_MATRIX("R", R, R.rows(), R.cols());
                    QPPP_TRACE_VECTOR("A", A, iq);
                    QPPP_TRACE_VECTOR("ci_status", ci_status, ci_status.size());

                    ci_status.assign(num_ci, ConstraintStatus::INACTIVE);
                    for (int i = num_ce; i < iq; ++i)
                    {
                        A[i-num_ce] = A_old[i-num_ce];
                        u[i] = u_old[i];
                        ci_status[A[i]] = ConstraintStatus::ACTIVE;
                    }
                    x = x_old;
                    goto l2; /* go to step 2 */
                }

                QPPP_TRACE_MATRIX("R", R, R.rows(), R.cols());
                QPPP_TRACE_VECTOR("A", A, iq);
                QPPP_TRACE_VECTOR("ci_status", ci_status, ci_status.size());

                goto l1;
            }

            /* a patial step has taken */
            QPPP_TRACE_MESSAGE("Partial step has taken " << t << std::endl);
            QPPP_TRACE_VECTOR("x", x, x.size());

            /* drop constraint l */
            ci_status[l] = ConstraintStatus::INACTIVE;
            delete_constraint(u, iq, l, num_var, num_ce);
            --iq;

            QPPP_TRACE_MATRIX("R", R, R.rows(), R.cols());
            QPPP_TRACE_VECTOR("A", A, iq);

            /* update s[ip] = CI * x + ci0 */
            multiply_and_add_i(ci_violations, CI, x, ci0, ip);

            QPPP_TRACE_VECTOR("s", ci_violations, num_ci);

            goto l2a;
        }
};


Solver::Solver()
{
    impl = new Implementation();
}

Solver::~Solver()
{
    delete impl;
}

Status::Value   Solver::solve(
        QPPP_MATRIX(double)& G,
        QPPP_VECTOR(double)& g0,
        const QPPP_MATRIX(double)& CE,
        const QPPP_VECTOR(double)& ce0,
        const QPPP_MATRIX(double)& CI,
        const QPPP_VECTOR(double)& ci0,
        QPPP_VECTOR(double)& x)
{
    return (impl->solve(G, g0, CE, ce0, CI, ci0, x));
}
} // namespace QuadProgpp
