# Born_approximation
 
Born approximation programs in Julia

Please check that you have installed the required packages and revise the paths
in the programs. 


Some programs require functions defined in funciones_born.jl. 
Check that the path to this file is correct in the programs

The list of programs for the Born approximation are:
- Born_ap_potential_stepfun.jl     (piecewise potential)
- Born_ap_potential_stepfun_bump.jl  (piecewise bump potential)
- Born_ap_potential_smooth.jl  (smooth potential)

The list of programs for the iterative algorithm are:
- Iterative_method_smooth.jl (compute the different iterations for smooth potential)
- Iterative_method_step.jl (compute the different iterations for step potential)
- Figure_iter.jl (From the data files created by the previous two programs 
  generates graphics with an error comparison)
