# Compute Born approximation of a potential
## functions are defined in funciones_born.jl
# The following packages must be installed
using ArbNumerics, Flux, SpecialFunctions, FFTW
using JLD, Plots, Measures, LaTeXStrings

# Include the directory where functions_born.jl is located
include("C:/Users/Carlo/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_1/funciones_born.jl")
using .funciones_born
setextrabits(0)
setworkingprecision(ArbFloat, bits = 1024) # digits = 
# Include your directiory where figures will be saved
cd("C:/Users/Carlo/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_1")

## Define potential function 
h=0.001;
r=0:h:1; # radial variable
# steps at d1 and d2. q takes values
# q(x)= 2a, if x in (0,d1)
# q(x)= a, if x in (d1,d2) 
# q(x)= 0, if x in (d2,1)
d1=1/3;
d2=2/3;
a=1.0
fg=a*any.(x->(x<d2),r)
fg=fg .+ a*any.(x->(x<d1),r)

## Compute q(r) piecewise constant: rj, ga
N_s=3;          # number of steps
rj=(1:N_s)/N_s; # discrete radius
ga=[2*a; a; 0]

## Compute eigenvalues
N_e=1000;
I,la=general_eig_fun(rj,ga,N_e) # lambda_k-k

## Compute t_exp
xi_max=460;
xi=collect(0:ArbFloat(pi)/10:xi_max);
ftq_exp = zeros(length(xi))*ArbFloat(0)
for k=1:N_e       
    global ftq_exp += (-1)^(k-1)*(xi./2).^(2*(k-1))./(gamma(ArbFloat(k))*(gamma(ArbFloat(k+1/2)))).*(la[k]);          
end
ftq_exp=2*(ArbFloat(pi))^(3/2)*ftq_exp;

## Compute Fourier transform of q
fr=fg
r_1=0:h:10
fr_1=zeros(length(r_1))
fr_1[1:length(fr)]=fr;
Fq,qp=ft_radial_fun(r_1,h,fr_1);
Fq[1]=Fq[2]; # avoids singularity

## Invert Fourier transform to obtain q_exp, 
delq=qp[2]-qp[1];
ftq_exp2=zeros(length(qp));
ftq_exp2[1:length(ftq_exp)]=ftq_exp;
q_exp,rp=ift_radial_fun(qp,delq,ftq_exp2);

## Plot results
p3=plot(qp,Fq,label = "FT \$ q \$", linewidth=2,linestyle = :dash)
plot!(qp,ftq_exp2,xlims=(0,50),label = "FT \$q_{exp}\$",linewidth=2,
xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
bottom_margin=5mm, right_margin=10mm)
savefig("exp1_pot_disc_FT.png")

p4=plot(r,fr,label = L"q",linestyle = :dash,linewidth=2,
xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16)
plot!(rp,q_exp,xlims=(0,1),label = L"q_{exp}",linewidth=2,bottom_margin=5mm, right_margin=10mm,
      xlabel = "r",ylims=(-1,2.4),
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16)
display(p4)
savefig("exp1_pot_disc.png")
 
