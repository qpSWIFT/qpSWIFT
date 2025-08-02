
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>  // Add this line for std::vector conversion
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
#include <iostream>
#include <queue>
#include <mutex>
#include <future>
#include <condition_variable>
#include <functional>
#include "qpSWIFT.h"
#include "GlobalOptions.h" // Include GlobalOptions for the MAXIT, RELTOL, etc. macros

namespace py = pybind11;
using namespace Eigen;
using SpMat = Eigen::SparseMatrix<qp_real, Eigen::ColMajor, qp_int>;


class ThreadPool {
public:
    ThreadPool(size_t threads)
        :   stop(false)
    {
        for (size_t i = 0; i < threads; ++i)
            workers.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
    }

    // enqueue a job, get a future for its return value
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using Ret = typename std::result_of<F(Args...)>::type;
        auto task = std::make_shared<std::packaged_task<Ret()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<Ret> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop) throw std::runtime_error("enqueue on stopped ThreadPool");
            tasks.emplace([task](){ (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread>        workers;
    std::queue<std::function<void()>> tasks;
    std::mutex                      queue_mutex;
    std::condition_variable         condition;
    bool                            stop;
};

// Structure to hold solution data
struct __attribute__((visibility("hidden"))) QPSolution {
//    py::array_t<double> x;
    std::vector<double> raw_x;
    int exit_flag;
    int iterations;
    double setup_time;
    double solve_time;
    double obj_value;
};


// Structure to hold solver options
struct QPOptions {
    qp_int maxiter;
    qp_real abstol;
    qp_real reltol;
    qp_real sigma;
    qp_int verbose;
    
    // Constructor with default values from GlobalOptions.h
    QPOptions() : 
        maxiter(MAXIT),
        abstol(ABSTOL),
        reltol(RELTOL),
        sigma(SIGMA),
        verbose(VERBOSE) {}
};

// Utility functions
template <typename VectorType>
bool has_infinities(const VectorType& vec) {
    for (Eigen::Index i = 0; i < vec.size(); i++) {
        if (!std::isfinite(vec[i])) return true;
    }
    return false;
}

template <typename VectorType>
bool all_infinite(const VectorType& vec, bool lower_bound = true) {
    double inf_val = lower_bound ? -INFINITY : INFINITY;
    for (Eigen::Index i = 0; i < vec.size(); i++) {
        if (vec[i] != inf_val) return false;
    }
    return true;
}

// Process bounds and create qpSWIFT-compatible constraint matrices
template <typename MatrixType, typename VectorType1, typename VectorType2>
void process_range_constraints(
    const MatrixType& A,
    const VectorType1& lb,
    const VectorType2& ub,
    MatrixXd& G,
    VectorXd& h,
    bool is_box_constraint = false
) {
    int n_constraints = A.rows();
    int n_vars = is_box_constraint ? A.cols() : A.cols();  // Fixed: Always use cols for variable count
    int active_constraints = 0;

    // First, count how many actual constraints we have
    std::vector<int> lb_active_idx, ub_active_idx;

    for (int i = 0; i < n_constraints; i++) {
        if (std::isfinite(lb[i])) {
            lb_active_idx.push_back(i);
            active_constraints++;
        }
        if (std::isfinite(ub[i])) {
            ub_active_idx.push_back(i);
            active_constraints++;
        }
    }

    // Prepare G and h
    G.resize(active_constraints, n_vars);
    h.resize(active_constraints);

    // Fill G and h
    int constraint_idx = 0;

    // Add lower bound constraints: -A*x <= -lb
    for (int i : lb_active_idx) {
        if (is_box_constraint) {
            G.row(constraint_idx) = MatrixXd::Zero(1, n_vars);
            G(constraint_idx, i) = -1.0;
        } else {
            // For non-box constraints, copy the entire row
            for (int j = 0; j < n_vars; j++) {
                G(constraint_idx, j) = -A(i, j);
            }
        }
        h[constraint_idx] = -lb[i];
        constraint_idx++;
    }

    // Add upper bound constraints: A*x <= ub
    for (int i : ub_active_idx) {
        if (is_box_constraint) {
            G.row(constraint_idx) = MatrixXd::Zero(1, n_vars);
            G(constraint_idx, i) = 1.0;
        } else {
            // For non-box constraints, copy the entire row
            for (int j = 0; j < n_vars; j++) {
                G(constraint_idx, j) = A(i, j);
            }
        }
        h[constraint_idx] = ub[i];
        constraint_idx++;
    }
}

// Apply solver options to a QP struct
void apply_options(QP* myQP, const QPOptions& options) {
    if (myQP && myQP->options) {
        myQP->options->maxit = options.maxiter;
        myQP->options->abstol = options.abstol;
        myQP->options->reltol = options.reltol;
        myQP->options->sigma = options.sigma;
        myQP->options->verbose = options.verbose;
        if (options.verbose){
            std::cout << "\n****qpSWIFT : Sparse Quadratic Programming Solver****\n\n" << std::endl;
            std::cout << "================Settings Applied======================\n" << std::endl;
            std::cout << "Maximum Iterations :" << myQP->options->maxit << std::endl;
            std::cout << "ABSTOL             :" << myQP->options->abstol << std::endl;
            std::cout << "RELTOL             :" << myQP->options->reltol << std::endl;
            std::cout << "SIGMA              :" << myQP->options->sigma << std::endl;
            std::cout << "VERBOSE            :" << myQP->options->verbose << std::endl;
            std::cout << "Permutation vector : AMD Solver\n\n" << std::endl;
        }
    }
}

// Main solver function for dense matrices
QPSolution solve_qp_dense(
    const MatrixXd& H,
    const VectorXd& g,
    const VectorXd& lb,
    const VectorXd& ub,
    const MatrixXd& E,
    const VectorXd& b,
    const MatrixXd& A,
    const VectorXd& lbA,
    const VectorXd& ubA,
    const QPOptions& options
) {
    // Initialize solution
    QPSolution solution;
    qp_int n = H.cols();  // Number of variables

    // Prepare constraint matrices
    MatrixXd G_box, G_linear;
    VectorXd h_box, h_linear;

    // Process box constraints lb <= x <= ub
    MatrixXd identity = MatrixXd::Identity(n, n);

    // Handle empty vectors properly
    VectorXd lb_effective = lb;
    VectorXd ub_effective = ub;

    if (lb.size() == 0) {
        lb_effective = VectorXd::Constant(n, -INFINITY);
    } else if (lb.size() != n) {
        std::cout << "WARNING: lb size mismatch! Expected " << n << ", got " << lb.size() << std::endl;
        lb_effective = VectorXd::Constant(n, -INFINITY);
    }

    if (ub.size() == 0) {
        ub_effective = VectorXd::Constant(n, INFINITY);
    } else if (ub.size() != n) {
        std::cout << "WARNING: ub size mismatch! Expected " << n << ", got " << ub.size() << std::endl;
        ub_effective = VectorXd::Constant(n, INFINITY);
    }

    process_range_constraints(identity, lb_effective, ub_effective, G_box, h_box, true);

    // Process linear constraints lbA <= Ax <= ubA
    if (A.rows() > 0) {
        VectorXd lbA_effective = lbA;
        VectorXd ubA_effective = ubA;

        if (lbA.size() == 0) {
            lbA_effective = VectorXd::Constant(A.rows(), -INFINITY);
        } else if (lbA.size() != A.rows()) {
            std::cout << "WARNING: lbA size mismatch! Expected " << A.rows() << ", got " << lbA.size() << std::endl;
            lbA_effective = VectorXd::Constant(A.rows(), -INFINITY);
        }

        if (ubA.size() == 0) {
            ubA_effective = VectorXd::Constant(A.rows(), INFINITY);
        } else if (ubA.size() != A.rows()) {
            std::cout << "WARNING: ubA size mismatch! Expected " << A.rows() << ", got " << ubA.size() << std::endl;
            ubA_effective = VectorXd::Constant(A.rows(), INFINITY);
        }

        process_range_constraints(A, lbA_effective, ubA_effective, G_linear, h_linear);
    }

    // Combine constraints
    int total_ineq_constraints = G_box.rows() + G_linear.rows();
    MatrixXd G(total_ineq_constraints, n);
    VectorXd h(total_ineq_constraints);

    if (G_box.rows() > 0) {
        G.topRows(G_box.rows()) = G_box;
        h.head(G_box.rows()) = h_box;
    }

    if (G_linear.rows() > 0) {
        G.bottomRows(G_linear.rows()) = G_linear;
        h.tail(G_linear.rows()) = h_linear;
    }

    // Check equality constraints
    qp_int p = 0;
    MatrixXd E_effective;
    VectorXd b_effective;

    if (E.rows() > 0) {
        if (E.rows() != b.size()) {
            std::cout << "WARNING: E and b size mismatch! E rows: " << E.rows() << ", b size: " << b.size() << std::endl;
            E_effective = MatrixXd(0, n);
            b_effective = VectorXd(0);
        } else {
            E_effective = E;
            b_effective = b;
            p = E.rows();
        }
    } else {
        E_effective = MatrixXd(0, n);
        b_effective = VectorXd(0);
    }

    // Prepare qpSWIFT input
    qp_int m = G.rows();  // Number of inequality constraints

    // Create mutable copies for qpSWIFT (since it requires non-const pointers)
    MatrixXd Hessian_copy = H;
    MatrixXd G_copy = G;
    MatrixXd E_copy = E_effective;
    VectorXd g_copy = g;
    VectorXd h_copy = h;
    VectorXd b_copy = b_effective;

//    std::cout << "Setting up QP problem (dense)...=========================================" << std::endl;
//    std::cout << "Hessian matrix H:\n" << Hessian_copy << std::endl;
//    std::cout << "Gradient vector g:\n" << g_copy.transpose() << std::endl;
//    std::cout << "Equality matrix E:\n" << E_copy << std::endl;
//    std::cout << "Equality bounds b:\n" << b_copy.transpose() << std::endl;
//    std::cout << "Inequality matrix G:\n" << G_copy << std::endl;
//    std::cout << "Inequality bounds h:\n" << h_copy.transpose() << std::endl;

    QP *myQP = QP_SETUP_dense(
        n, m, p,
        Hessian_copy.data(),  // H matrix
        E_copy.data(),        // E matrix (equality constraints)
        G_copy.data(),        // G matrix (inequality constraints)
        g_copy.data(),        // g vector (linear cost)
        h_copy.data(),        // h vector (inequality bounds)
        b_copy.data(),        // b vector (equality bounds)
        NULL,                 // No permutation
        COLUMN_MAJOR_ORDERING    // Row major ordering
    );

    if (myQP == NULL) {
        std::cout << "ERROR: QP_SETUP_dense failed!" << std::endl;
        solution.exit_flag = -1;
        solution.iterations = 0;
        solution.setup_time = 0;
        solution.solve_time = 0;
        solution.obj_value = 0;
        return solution;
    }
    
    // Apply user options
    apply_options(myQP, options);

    qp_int result = QP_SOLVE(myQP);

    solution.raw_x = std::vector<double>(n);
    for (int i = 0; i < n; i++) {
        solution.raw_x[i] = myQP->x[i];
    }

    // Extract statistics
    solution.exit_flag = result;
    solution.iterations = myQP->stats->IterationCount;
    solution.setup_time = myQP->stats->tsetup;
    solution.solve_time = myQP->stats->tsolve;
    solution.obj_value = myQP->stats->fval;

    // Clean up
    QP_CLEANUP_dense(myQP);

    return solution;
}

QPSolution solve_qp_sparse(
    const SpMat& H,
    const VectorXd& g,
    const VectorXd& lb,
    const VectorXd& ub,
    const SpMat& E,
    const VectorXd& b,
    const SpMat& A,
    const VectorXd& lbA,
    const VectorXd& ubA,
    const QPOptions& options
) {
    // Initialize solution
    QPSolution solution;
    qp_int n = H.cols();  // Number of variables

    // Make sure all input matrices are in CSC (column-major) format
    SpMat H_csc = H;
    SpMat E_csc = E;
    SpMat A_csc = A;

    // Create constraint matrices
    std::vector<Eigen::Triplet<double>> G_triplets;
    std::vector<double> h_values;

    // Process box constraints lb <= x <= ub
    VectorXd lb_effective = lb;
    VectorXd ub_effective = ub;

    if (lb.size() == 0) {
        lb_effective = VectorXd::Constant(n, -INFINITY);
    } else if (lb.size() != n) {
        std::cout << "WARNING: lb size mismatch! Expected " << n << ", got " << lb.size() << std::endl;
        lb_effective = VectorXd::Constant(n, -INFINITY);
    }

    if (ub.size() == 0) {
        ub_effective = VectorXd::Constant(n, INFINITY);
    } else if (ub.size() != n) {
        std::cout << "WARNING: ub size mismatch! Expected " << n << ", got " << ub.size() << std::endl;
        ub_effective = VectorXd::Constant(n, INFINITY);
    }

    // Add box constraints
    int constraint_idx = 0;
    for (int i = 0; i < n; i++) {
        if (std::isfinite(lb_effective[i])) {
            G_triplets.push_back(Eigen::Triplet<double>(constraint_idx, i, -1.0));
            h_values.push_back(-lb_effective[i]);
            constraint_idx++;
        }
    }

    for (int i = 0; i < n; i++) {
        if (std::isfinite(ub_effective[i])) {
            G_triplets.push_back(Eigen::Triplet<double>(constraint_idx, i, 1.0));
            h_values.push_back(ub_effective[i]);
            constraint_idx++;
        }
    }

    // Process linear constraints lbA <= Ax <= ubA
    if (A_csc.rows() > 0) {
        // Prepare effective bounds
        VectorXd lbA_effective = lbA;
        VectorXd ubA_effective = ubA;

        if (lbA.size() == 0) {
            lbA_effective = VectorXd::Constant(A_csc.rows(), -INFINITY);
        } else if (lbA.size() != A_csc.rows()) {
            std::cout << "WARNING: lbA size mismatch! Expected " << A_csc.rows() << ", got " << lbA.size() << std::endl;
            lbA_effective = VectorXd::Constant(A_csc.rows(), -INFINITY);
        }

        if (ubA.size() == 0) {
            ubA_effective = VectorXd::Constant(A_csc.rows(), INFINITY);
        } else if (ubA.size() != A_csc.rows()) {
            std::cout << "WARNING: ubA size mismatch! Expected " << A_csc.rows() << ", got " << ubA.size() << std::endl;
            ubA_effective = VectorXd::Constant(A_csc.rows(), INFINITY);
        }

        // Add lower bound constraints: -A*x <= -lbA
        if (A_csc.rows() > 0) {
            // Create temporary CSR matrix for efficient row access
            Eigen::SparseMatrix<double, Eigen::RowMajor> A_csr(A_csc);

            // Add lower bound constraints: -A*x <= -lbA
            for (int i = 0; i < A_csr.rows(); i++) {
                if (std::isfinite(lbA_effective[i])) {
                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_csr, i); it; ++it) {
                        G_triplets.push_back(Eigen::Triplet<double>(constraint_idx, it.col(), -it.value()));
                    }
                    h_values.push_back(-lbA_effective[i]);
                    constraint_idx++;
                }
            }

            // Add upper bound constraints: A*x <= ubA
            for (int i = 0; i < A_csr.rows(); i++) {
                if (std::isfinite(ubA_effective[i])) {
                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_csr, i); it; ++it) {
                        G_triplets.push_back(Eigen::Triplet<double>(constraint_idx, it.col(), it.value()));
                    }
                    h_values.push_back(ubA_effective[i]);
                    constraint_idx++;
                }
            }
        }
    }

    // Create the G matrix and h vector
    int m = constraint_idx;  // Total number of inequality constraints
    SpMat G(m, n);

    G.setFromTriplets(G_triplets.begin(), G_triplets.end());
    G.makeCompressed();

    // Process equality constraints
    qp_int p = 0;
    SpMat E_effective(0, n);
    VectorXd b_effective;

    if (E_csc.rows() > 0) {
        if (E_csc.rows() != b.size()) {
            std::cout << "WARNING: E and b size mismatch! E rows: " << E_csc.rows() << ", b size: " << b.size() << std::endl;
        } else {
            E_effective = E_csc;
            b_effective = b;
            p = E_csc.rows();
        }
    } else {
        b_effective.resize(0);
    }

    // Print problem details
//    std::cout << "Setting up QP problem (sparse)...==================================" << std::endl;
//    std::cout << "Hessian matrix H:\n" << H << std::endl;
//    std::cout << "Gradient vector g:\n" << g.transpose() << std::endl;
//    std::cout << "Equality matrix E:\n" << E_effective << std::endl;
//    std::cout << "Equality bounds b:\n" << b_effective.transpose() << std::endl;
//    std::cout << "Inequality matrix G:\n" << G << std::endl;
//    std::cout << "Inequality bounds h:\n" << h.transpose() << std::endl;

    // Create mutable copy of gradient
    VectorXd g_copy = g;

    QP *myQP = QP_SETUP(
        n, m, p,
        H_csc.outerIndexPtr(), H_csc.innerIndexPtr(), H_csc.valuePtr(),
        E_effective.outerIndexPtr(), E_effective.innerIndexPtr(), E_effective.valuePtr(),
        G.outerIndexPtr(), G.innerIndexPtr(), G.valuePtr(),
        g_copy.data(),
        h_values.data(),
        b_effective.data(),
        0.0,
        NULL  // No permutation
    );

    if (myQP == NULL) {
        std::cout << "ERROR: QP_SETUP failed!" << std::endl;
        solution.exit_flag = -1;
        solution.iterations = 0;
        solution.setup_time = 0;
        solution.solve_time = 0;
        solution.obj_value = 0;
        return solution;
    }

    apply_options(myQP, options);

    qp_int result = QP_SOLVE(myQP);

    solution.raw_x = std::vector<double>(n);
    for (int i = 0; i < n; i++) {
        solution.raw_x[i] = myQP->x[i];
    }

    solution.exit_flag = result;
    solution.iterations = myQP->stats->IterationCount;
    solution.setup_time = myQP->stats->tsetup;
    solution.solve_time = myQP->stats->tsolve;
    solution.obj_value = myQP->stats->fval;

    // Clean up
    QP_CLEANUP(myQP);

    return solution;
}

// Helper to convert Python dict to QPOptions
QPOptions dict_to_options(const py::dict& options_dict) {
    QPOptions options;
    
    if (options_dict.contains("MAXITER"))
        options.maxiter = options_dict["MAXITER"].cast<qp_int>();
    else options.maxiter = MAXIT;
    
    if (options_dict.contains("ABSTOL"))
        options.abstol = options_dict["ABSTOL"].cast<qp_real>();
    else options.abstol = ABSTOL;
    
    if (options_dict.contains("RELTOL"))
        options.reltol = options_dict["RELTOL"].cast<qp_real>();
    else options.reltol = RELTOL;
    
    if (options_dict.contains("SIGMA"))
        options.sigma = options_dict["SIGMA"].cast<qp_real>();
    else options.sigma = SIGMA;
    
    if (options_dict.contains("VERBOSE"))
        options.verbose = options_dict["VERBOSE"].cast<qp_int>();
    else options.verbose = VERBOSE;
    
    return options;
}

static ThreadPool pool(std::max<unsigned>(1u, std::thread::hardware_concurrency() - 1));

// Python module
PYBIND11_MODULE(qpSWIFT_sparse_bindings, m) {
    m.doc() = "Custom sparse python bindings for the qpSWIFT solver";

    py::class_<QPSolution>(m, "QPSolution")
        .def_readonly("exit_flag", &QPSolution::exit_flag, "Exit flag (0: success, >0: error)")
        .def_readonly("iterations", &QPSolution::iterations, "Number of iterations")
        .def_readonly("setup_time", &QPSolution::setup_time, "Setup time in seconds")
        .def_readonly("solve_time", &QPSolution::solve_time, "Solve time in seconds")
        .def_readonly("obj_value", &QPSolution::obj_value, "Objective function value")
        .def_property_readonly("x", [](QPSolution &self) {
            ssize_t N = self.raw_x.size();
            return py::array_t<double>(
              { N },
              { (ssize_t)sizeof(double) },
              self.raw_x.data(),
              py::cast(&self)
            );
          });

    // Dense version with all constraints and options
    m.def("solve",
        [](const MatrixXd& H, const VectorXd& g,
           const VectorXd& lb, const VectorXd& ub,
           const MatrixXd& E, const VectorXd& b,
           const MatrixXd& A, const VectorXd& lbA, const VectorXd& ubA,
           const py::dict& options_dict) {
            QPOptions options = dict_to_options(options_dict);
            {
                py::gil_scoped_release release;
                QPSolution sol = solve_qp_dense(H, g, lb, ub, E, b, A, lbA, ubA, options);
                return sol;
            }
        },
        py::arg("H"), py::arg("g"),
        py::arg("lb") = VectorXd(), py::arg("ub") = VectorXd(),
        py::arg("E") = MatrixXd(0, 0), py::arg("b") = VectorXd(),
        py::arg("A") = MatrixXd(0, 0), py::arg("lbA") = VectorXd(), py::arg("ubA") = VectorXd(),
        py::arg("options") = py::dict(),
        "Solve a QP with dense matrices"
    );

    // Sparse version with all constraints and options
    m.def("solve_sparse",
        [](const SpMat& H, const VectorXd& g,
           const VectorXd& lb, const VectorXd& ub,
           const SpMat& E, const VectorXd& b,
           const SpMat& A, const VectorXd& lbA, const VectorXd& ubA,
           const py::dict& options_dict) {
            QPOptions options = dict_to_options(options_dict);
            {
                py::gil_scoped_release release;
                QPSolution sol = solve_qp_sparse(H, g, lb, ub, E, b, A, lbA, ubA, options);
                return sol;
            }
        },
        py::arg("H"), py::arg("g"),
        py::arg("lb") = VectorXd(), py::arg("ub") = VectorXd(),
        py::arg("E") = SpMat(), py::arg("b") = VectorXd(),
        py::arg("A") = SpMat(), py::arg("lbA") = VectorXd(), py::arg("ubA") = VectorXd(),
        py::arg("options") = py::dict(),
        "Solve a QP with sparse matrices"
    );


    m.def("run_sparse",
        [](const VectorXd& P, const VectorXd& c,
           const VectorXd& h, const SpMat& G,
           const VectorXd& b, const SpMat& A,
           const py::dict& opts) {
            QPOptions options = dict_to_options(opts);
            QPSolution sol;
            {
                py::gil_scoped_release release;

                qp_int n = P.size();
                SpMat H_diag(n, n);

                 // Create triplet list for diagonal elements
                std::vector<Eigen::Triplet<double>> triplets;
                triplets.reserve(n);

                for (int i = 0; i < n; i++) {
                    triplets.push_back(Eigen::Triplet<double>(i, i, P(i)));
                }

                H_diag.setFromTriplets(triplets.begin(), triplets.end());
                H_diag.makeCompressed();

                sol = solve_qp_sparse(H_diag,
                                                 c,
                                                 VectorXd::Constant(P.size(), -INFINITY),
                                                 VectorXd::Constant(P.size(), INFINITY),
                                                 A, b,
                                                 G, VectorXd::Constant(h.size(), -INFINITY),
                                                 h, options);

            }
            // Convert to Python dict format
            py::dict result;
            py::dict basic_info;
            basic_info["ExitFlag"] = sol.exit_flag;
            result["basicInfo"] = basic_info;

            py::array_t<double> x_array = py::cast(sol.raw_x);
            result["sol"] = x_array;


            return result;
        },
        py::arg("P"), py::arg("c"),
        py::arg("h") = VectorXd(), py::arg("G") = SpMat(),
        py::arg("b") = VectorXd(), py::arg("A") = SpMat(),
        py::arg("opts") = py::dict(),
        "Solve a QP with sparse matrices with old interface"
    );

    // Sparse version with all constraints and options
    m.def("run_sparse_with_box_constraints",
        [](const VectorXd& P, const VectorXd& c, const VectorXd& h_box, const VectorXd& h,
           const SpMat& G_box, const SpMat& G,
           const VectorXd& b, const SpMat& A,
           const py::dict& opts) {
            QPOptions options = dict_to_options(opts);
            QPSolution sol;
            {
                py::gil_scoped_release release;

                qp_int n = P.size();
                SpMat H_diag(n, n);

                 // Create triplet list for diagonal elements
                std::vector<Eigen::Triplet<double>> triplets;
                triplets.reserve(n);

                for (int i = 0; i < n; i++) {
                    triplets.push_back(Eigen::Triplet<double>(i, i, P(i)));
                }

                H_diag.setFromTriplets(triplets.begin(), triplets.end());
                H_diag.makeCompressed();

                // Stack G_box and G vertically
                SpMat G_combined;
                if (G_box.rows() > 0 && G.rows() > 0) {
                    // Both matrices have data - stack them vertically
                    G_combined.resize(G_box.rows() + G.rows(), n);
                    std::vector<Eigen::Triplet<double>> combined_triplets;

                    // Add G_box triplets
                    for (int k = 0; k < G_box.outerSize(); ++k) {
                        for (SpMat::InnerIterator it(G_box, k); it; ++it) {
                            combined_triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
                        }
                    }

                    // Add G triplets with row offset
                    for (int k = 0; k < G.outerSize(); ++k) {
                        for (SpMat::InnerIterator it(G, k); it; ++it) {
                            combined_triplets.push_back(Eigen::Triplet<double>(it.row() + G_box.rows(), it.col(), it.value()));
                        }
                    }

                    G_combined.setFromTriplets(combined_triplets.begin(), combined_triplets.end());
                    G_combined.makeCompressed();
                } else if (G_box.rows() > 0) {
                    // Only G_box has data
                    G_combined = G_box;
                } else if (G.rows() > 0) {
                    // Only G has data
                    G_combined = G;
                } else {
                    // Neither has data - create empty matrix
                    G_combined.resize(0, n);
                }

                // Combine h_box and h vectors
                VectorXd h_combined;
                if (h_box.size() > 0 && h.size() > 0) {
                    // Both vectors have data - concatenate them
                    h_combined.resize(h_box.size() + h.size());
                    h_combined.head(h_box.size()) = h_box;
                    h_combined.tail(h.size()) = h;
                } else if (h_box.size() > 0) {
                    // Only h_box has data
                    h_combined = h_box;
                } else if (h.size() > 0) {
                    // Only h has data
                    h_combined = h;
                } else {
                    // Neither has data - create empty vector
                    h_combined.resize(0);
                }

                sol = solve_qp_sparse(H_diag,
                                                 c,
                                                 VectorXd::Constant(P.size(), -INFINITY),
                                                 VectorXd::Constant(P.size(), INFINITY),
                                                 A, b,
                                                 G_combined, VectorXd::Constant(h_combined.size(), -INFINITY),
                                                 h_combined, options);

            }
            // Convert to Python dict format
            py::dict result;
            py::dict basic_info;
            basic_info["ExitFlag"] = sol.exit_flag;
            result["basicInfo"] = basic_info;

            py::array_t<double> x_array = py::cast(sol.raw_x);
            result["sol"] = x_array;


            return result;
        },
        py::arg("P"), py::arg("c"), py::arg("h_box") = VectorXd(), py::arg("h") = VectorXd(),
        py::arg("G_box") = SpMat(), py::arg("G") = SpMat(),
        py::arg("b") = VectorXd(), py::arg("A") = SpMat(),
        py::arg("opts") = py::dict(),
        "Solve a QP with sparse matrices with old interface"
    );

    m.def("solve_sparse_H_diag",
        [](const VectorXd& H, const VectorXd& g,
           const VectorXd& lb, const VectorXd& ub,
           const SpMat& E, const VectorXd& b,
           const SpMat& A, const VectorXd& lbA, const VectorXd& ubA,
           const py::dict& options_dict) {
            QPOptions options = dict_to_options(options_dict);
            {
                py::gil_scoped_release release;

                qp_int n = H.size();
                SpMat H_diag(n, n);

                 // Create triplet list for diagonal elements
                std::vector<Eigen::Triplet<double>> triplets;
                triplets.reserve(n);

                for (int i = 0; i < n; i++) {
                    triplets.push_back(Eigen::Triplet<double>(i, i, H(i)));
                }

                H_diag.setFromTriplets(triplets.begin(), triplets.end());
                H_diag.makeCompressed();


                QPSolution sol = solve_qp_sparse(H_diag, g, lb, ub, E, b, A, lbA, ubA, options);
                return sol;
            }
        },
        py::arg("H"), py::arg("g"),
        py::arg("lb") = VectorXd(), py::arg("ub") = VectorXd(),
        py::arg("E") = SpMat(), py::arg("b") = VectorXd(),
        py::arg("A") = SpMat(), py::arg("lbA") = VectorXd(), py::arg("ubA") = VectorXd(),
        py::arg("options") = py::dict(),
        "Solve a QP with a diagonal Hessian matrix (provided as a vector of diagonal elements)"
    );

    m.def("solve_sparse_H_diag_batch",
        [](const std::vector<py::tuple>& qp_data,
           const py::dict& options_dict)
        {
            // parse options once
            QPOptions options = dict_to_options(options_dict);

            size_t N = qp_data.size();
            // POD for each QP
            struct QP {
                    VectorXd H, g, lb, ub, b, lbA, ubA;
                    SpMat E, A;
                };
            std::vector<QP> data;
            data.reserve(N);
            // unpack under GIL
            for (auto &t : qp_data) {
                    if (t.size() != 9)
                    throw std::invalid_argument("Each tuple needs 9 elements");
                data.push_back(QP{
                        t[0].cast<VectorXd>(),
                        t[1].cast<VectorXd>(),
                        t[2].cast<VectorXd>(),
                        t[3].cast<VectorXd>(),
                        t[5].cast<VectorXd>(),
                        t[7].cast<VectorXd>(),
                        t[8].cast<VectorXd>(),
                        t[4].cast<SpMat>(),
                        t[6].cast<SpMat>()
                    });
            }

            std::vector<QPSolution> solutions(N);
            std::vector<std::future<void>> futures;
            futures.reserve(N);

            {
                py::gil_scoped_release release;
                // enqueue one job per QP
                for (size_t i = 0; i < N; ++i) {
                    futures.emplace_back(pool.enqueue([&,i] {
                            const auto& q = data[i];
                            // build diagonal H
                            qp_int n = q.H.size();
                            SpMat Hdiag(n, n);
                            std::vector<Eigen::Triplet<double>> T;
                            T.reserve(n);
                            for (qp_int k = 0; k < n; ++k)
                                T.emplace_back(k,k,q.H(k));
                            Hdiag.setFromTriplets(T.begin(), T.end());
                            Hdiag.makeCompressed();

                            solutions[i] = solve_qp_sparse(
                                Hdiag, q.g,
                                q.lb, q.ub,
                                q.E,  q.b,
                                q.A,
                                q.lbA, q.ubA,
                                options
                            );
                        })
                    );
                    }
                for (auto &f : futures) f.get();
            }

            return solutions;
            },
        py::arg("qp_data"),
        py::arg("options") = py::dict(),
        "Solve a batch of diagonal‐H QPs in parallel from a list of 9‐tuples "
        "(H, g, lb, ub, E, b, A, lbA, ubA) with a shared options dict."
    );
}