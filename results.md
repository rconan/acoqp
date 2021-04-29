# Active Optics Quadratic Programming Algorithm

# Example case:

```
Number of lenslets:13824
Valid lenslets:7360
rho_3:0.1
Number of lenslets:13824
Number of Tu cols:271
-----------------------------------------------------------------
           OSQP v0.6.2  -  Operator Splitting QP Solver
              (c) Bartolomeo Stellato,  Goran Banjac
        University of Oxford  -  Stanford University 2021
-----------------------------------------------------------------
problem:  variables n = 271, constraints m = 1228
          nnz(P) + nnz(A) = 369644
settings: linear system solver = qdldl,
          eps_abs = 1.0e-08, eps_rel = 1.0e-06,
          eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
          rho = 1.00e-01 (adaptive),
          sigma = 1.00e-06, alpha = 1.60, max_iter = 135500
          check_termination: on (interval 25),
          scaling: on, scaled_termination: off
          warm start: on, polish: off, time_limit: off

iter   objective    pri res    dua res    rho        time
   1  -2.6383e-08   9.66e-20   1.50e-03   1.00e-01   1.27e-01s
  75  -3.6471e-08   7.06e-19   1.79e-11   1.00e-06   2.84e-01s

status:               solved
number of iterations: 75
optimal objective:    -0.0000
run time:             2.84e-01s
optimal rho estimate: 1.00e-06

-1.5061e-7
1.2179e-7
-1.5978e-7
1.6365e-6
-2.3834e-6
-5.8098e-7
-2.1422e-7
```
