# Compute Born approximation of a potential
## functions are defined in funciones_born.jl
# The following packages must be installed
using ArbNumerics, Flux, SpecialFunctions, FFTW
using JLD, Plots
# Include the directory where functions_born.jl is located
include("C:/Users/Carlo/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_1/funciones_born.jl")
using .funciones_born
setextrabits(0)
setworkingprecision(ArbFloat, bits = 256) # digits = 
# Include your directiory where figures will be saved
cd("C:/Users/carlo/ownCloud/trabajo/papers/Cristobal/julia/potencial/example_2_smooth")

## Define potential function 
h=0.0001;
r=0:h:1; # radial variable
d=0.9;
Psi=exp.(-(r.^2)./(((r.^2).-d^2).^2)).*(any.(x->(x<d),r))
alpha=0.3;
fg=(((alpha*Psi.+1).^2).-1);
fg=fg*1

## Compute q(r) piecewise constant: rj, ga
N_s=1000;        # number of steps
rj=(1:N_s)/N_s; # discrete radius
ga=zeros(N_s)
valu=[findmin(abs.(r .- rj[jj]))[2] for jj in 1:N_s]
ga = fg[valu]

## Compute eigenvalues
N_e=400;
I,la=general_eig_fun(rj,ga,N_e) # lambda_k-k

## Compute t_exp
xi_max=160;
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
plot!(qp,ftq_exp2,xlims=(0,50),label = "FT \$q_{exp}\$",xlabel = L"\xi",linewidth=2,
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     bottom_margin=5mm, right_margin=10mm)
savefig("exp2_pot_smooth_FT.png")

p4=plot(r,fr,label = L"q",linestyle = :dash,linewidth=2,
xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16)
plot!(rp,q_exp,label = L"q_{exp}", xlabel = "r",linewidth=2,xlims=(0,1),
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     bottom_margin=5mm, right_margin=10mm,ylims=(-0.6,0.8))
display(p4)

savefig("exp2_pot_smooth.png")


