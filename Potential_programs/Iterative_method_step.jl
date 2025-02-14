# Iterative algorithm to approximate the potential
## functions are defined in funciones_born.jl
# The following packages must be installed
using ArbNumerics, Flux, SpecialFunctions, FFTW
using JLD, Plots, Measures, LaTeXStrings
# Include the directory where functions_born.jl is located
include("C:/Users/Carlos/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_1/funciones_born.jl")
using .funciones_born
setextrabits(0)
setworkingprecision(ArbFloat, bits = 1024) # digits = 
# Include your directiory where figures will be saved
cd("C:/Users/carlos/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_4_iter")

#############################
## Definimos la función potencial continua: fg(r) (artículo Knudsen)
h=0.001;
r=0:h:1; # radial variable
d1=1/3;
d2=2/3;
fg=0.5*any.(x->(x<d2),r)
fg=fg .+ 0.5*any.(x->(x<d1),r)

## Compute q(r) piecewise constant: rj, ga
N_s=3;        # number of steps
rj=(1:N_s)/N_s; # discrete radius
ga=[1; 0.5; 0]

## Compute eigenvalues
N_e=400;
I,la=funciones_born.general_eig_fun(rj,ga,N_e) # lambda_k-k
N_m=N_e; # number of moments (used later)

## Compute t_exp
xi_max=160;
xi=collect(0:ArbFloat(pi)/10:xi_max);
ftq_exp = zeros(length(xi))*ArbFloat(0)
for k=1:N_e       
    global ftq_exp += (-1)^(k-1)*(xi./2).^(2*(k-1))./(gamma(ArbFloat(k))*(gamma(ArbFloat(k+1/2)))).*(la[k]);          
end
ftq_exp=2*(ArbFloat(pi))^(3/2)*ftq_exp;

## Invert Fourier transform to obtain q_exp, 
fr=fg
delq=qp[2]-qp[1];
ftq_exp2=zeros(length(qp));
ftq_exp2[1:length(ftq_exp)]=ftq_exp;
q_exp,rp=funciones_born.ift_radial_fun(qp,delq,ftq_exp2);

# restrict q_exp to (0,1)
r_real=rp[any.(x->(x<=1.01),rp)]
q_exp_real=q_exp[any.(x->(x<=1.01),rp)]

## Compute q(r) piecewise constant: rj, ga
N_s=60;        # number of steps

## q_1(=q_exp) -- > q_2

function iter_fun(q_1,r,N_s,N_e,N_m,xi,qp,delq,ftq_exp)
    # 1. transform q_1 in a step function
    rj=(1:N_s)/N_s; # discrete radius
    valu=[findmin(abs.(r .- rj[jj]))[2] for jj in 1:N_s]
    q1_step = q_1[valu] # values at the steps
    # 2. Compute eigenvalues and momenta
    I1,la_q1=funciones_born.general_eig_fun(rj,q1_step,N_e) # lambda_k-k
    Im1,mom_q1=funciones_born.general_mom_fun(rj,q1_step,N_m)
    # 3. Compute t_2-t_exp Arb_precision
    ftq_2 = zeros(length(xi))*ArbFloat(0)
    for k=1:N_e       
        ftq_2 += (-1)^(k-1)*(xi/2).^(2*(k-1))/(gamma(ArbFloat(k))*gamma(ArbFloat(k+1/2)))*(la_q1[k]-mom_q1[k]);          
    end
    ftq_2=2*(ArbFloat(pi))^(3/2)*ftq_2;
    # 4. compute q_2
    ftq_exp22 = zeros(length(qp))
    ftq_exp22[1:length(ftq_exp)] = -ftq_2 + ftq_exp
    q_2,rp=funciones_born.ift_radial_fun(qp,delq,ftq_exp22);
    return q_2,ftq_exp22
end

# plot results
p4= plot(r,fr,label = L"q",linewidth=2,linestyle = :dash,
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     xlabel = "r", bottom_margin=5mm, right_margin=10mm)
 plot!(rp,q_exp,xlims=(0,1),label = L"q_{exp}",legendfontsize=16,linewidth=2)

n_iter= 6
err_lin=zeros(n_iter+1)
err_l2=zeros(n_iter+1)

err_l2[1] = sqrt(h * sum((q_exp[1:length(r)]-fr).^2))
err_lin[1],ind_i = findmax(abs.(q_exp[1:length(r)]-fr))

q_1 = q_exp_real
for kk=1:n_iter
    global q_1,ftq_exp22=iter_fun(q_1,r_real,N_s,N_e,N_m,xi,qp,delq,ftq_exp)
    if kk <=2
      p4=plot!(rp,q_1,xlims=(0,1),label = "\$ q^$kk \$",linewidth=2,legendfontsize=16)
    end
    err_l2[kk+1] = sqrt(h * sum((q_1[1:length(r)]-fr).^2))
    err_lin[kk+1],ind_i = findmax(abs.(q_1[1:length(r)]-fr))
end
display(p4)
savefig("ejem_pot_iter_d.png")

# save results
FileIO.save("results_iter_d.jld2","err_l2",err_l2,"err_lin",err_lin)



