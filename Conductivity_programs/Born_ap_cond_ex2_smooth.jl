# Compute Born approximation of a conductivity (Calderon pb)
## functions are defined at the beginning
# The following packages must be installed
using ArbNumerics, Flux, SpecialFunctions, FFTW
using JLD, Plots, Measures, LaTeXStrings
setextrabits(0)
setworkingprecision(ArbFloat, bits = 1024) # digits = 
###################
###### Functions ############

function ft_radial_fun(r,h,Fr)
    # Compute the Fourier transform of a 3d radial function
    # create a function for positive r
    n = length(r)-1;
    delr =h;
    r = (-n:n)*delr;
    rpos = r[any.(x->x>=0,r)]  # r(r>=0);
    Frpos = Fr;    # F(r), r >=0
    rFrpos = rpos.*Frpos;
    # create antisymmetric function with r=0 at the center
    # see plot 2
    rFr = [-reverse(rFrpos[2:end]); rFrpos];
    
    # transform to q, multiply by delr to approximate the Riemmann integral 
    N = length(r);
    qFq = -2*pi*imag(fftshift(fft(ifftshift(rFr))))*delr;
    # The formula has 4*pi but we are summing two times since we consider the
    # odd part of qFq
    delq = 2*pi/(N*delr);   # golden rule for fft
    q = (-n:n)*delq;
    Fq_sim=qFq./q;
    Fq=Fq_sim[Int(((N+1)/2)):N];
    qp=(0:n)*delq;
    return Fq,qp
    end
    
    #################
    
    function ift_radial_fun(qp,delq,Fq)
    # Compute the inverse Fourier transform of a 3d radial function
    # create a function for positive r
    n = length(qp)-1;
    qp = (-n:n)*delq;
    qpos = qp[any.(x->(x>=0),qp)]  #qp(qp>=0);
    Fqpos = Fq;    # F(r), r >=0
    qFqpos = qpos.*Fqpos;
    # create antisymmetric function with r=0 at the center
    # see plot 2
    qFq = [-reverse(qFqpos[2:end]); qFqpos];
    
    N = length(qp);
    rFr = (1/(2*pi)^2)*N*imag(fftshift(ifft(ifftshift(qFq))))*delq;
    
    delr = 2*pi/(N*delq);   # golden rule for fft
    r = (-n:n)*delr;
    
    Fr_sim=rFr./r;
    Fr=Fr_sim[Int((N+1)/2):N];
    Fr[1]=Fr[2];
    rp=(0:n)*delr;
    return Fr,rp
    end
    
    ########### Funciones Juan

    function trapz(r,f)
        N=length(f)
        int=(f[1]+f[2])/2*(r[2]-r[1])
        for jj = 2:(N-1)
            int = int + (f[jj]+f[jj+1])/2*(r[jj+1]-r[jj]) 
        end
        return int
    end

    
    # Eigenvalues of the DtN d=3: radial case
    function general_eig_fun(rj,ga,N_e)
        N=length(ga)
        rj=ArbFloat.(rj) 
        ga=ArbFloat.(ga)
        la=ArbFloat.(zeros(N_e))
        
        v=[1;0]
        for k=0:N_e-1
            for jj=1:N-1
                M_l = [rj[jj]^k rj[jj]^(-(k+1)) ; ga[jj]*k*rj[jj]^(k-1) -ga[jj]*(k+1)*rj[jj]^(-(k+2))]
                M_r = [rj[jj]^k rj[jj]^(-(k+1)) ; ga[jj+1]*k*rj[jj]^(k-1) -ga[jj+1]*(k+1)*rj[jj]^(-(k+2))]         
                v = M_r\(M_l*v);
            end
           la[k+1]=ga[N]*(k*v[1]-(k+1)*v[2])/(v[1]+v[2])-k;   
         end
        I=0:(N_e-1)
        la[isnan.(la)] .= 0.0
        return I,la
    end
    
    # momenta associated with a potential with n steps
    function general_mom_fun(r,gm,m)
    
        n=length(r);
        mom=ArbFloat.(zeros(m))
        r=ArbFloat.(r)
        gm=ArbFloat.(gm)
        for kk=0:m-1
            mom[kk+1]=gm[1]*r[1]^(2*kk+3)/(2*kk+3);
            for jj=2:n
                mom[kk+1] += gm[jj]*(r[jj]^(2*kk+3)/(2*kk+3)-r[jj-1]^(2*kk+3)/(2*kk+3));
            end
        end
        I=0:m-1;
        return I,mom
    end
    
    # momenta associated with a potential with n steps
    function general_mom_fun_ne(r,gm,m)
    
    mom = zeros(m)    
    for jj=1:m
        mom[jj]=trapz(r,r.^(2*(jj-1)+2).*gm);
    end
    I=0:m-1
    return I,mom
    end

############################# End Functions

## Define the conductivity
h=0.00001;
r=0:h:1; # radial variable
d=0.9;
Psi=exp.(-(r.^2)./(((r.^2).-d^2).^2)).*(any.(x->(x<d),r))
alpha=0.3 # Test different values: 0.3, 1, -1;
fg=(alpha*Psi.+1).^2;

## Compute ga(r) piecewise constant: rj, ga
N_s=15000;        # number of steps
#N_s=5;
rj=(1:N_s)/N_s; # discrete radius
ga=zeros(N_s)
valu=[findmin(abs.(r .- rj[jj]))[2] for jj in 1:N_s]
ga = fg[valu]

## Compute eigenvalues
N_e=600;
I,la=general_eig_fun(rj,ga,N_e) # lambda_k-k

## Compute t_exp
xi_max=260;
xi=collect(0:ArbFloat(pi)/10:xi_max);
ftq_exp = zeros(length(xi))*ArbFloat(0)
for k=2:N_e  # la[1]=0       
    global ftq_exp += (-1)^(k-1)*(xi./2).^(2*(k-2))./(gamma(ArbFloat(k))*(gamma(ArbFloat(k+1/2)))).*(la[k]);
    #global ftq_exp += (-1)^(k-1)*(xi./2).^(2*(k-1))./(gamma(ArbFloat(k))*(gamma(ArbFloat(k+1/2)))).*(la[k]);          
end
ftq_exp=-(ArbFloat(pi))^(3/2)*ftq_exp;

## Compute Fourier transform of q
fr=fg
r_1=0:h:10
fr_1=zeros(length(r_1))
fr_1[1:length(fr)]=(fr-ones(length(fr)));
Fq,qp=ft_radial_fun(r_1,h,fr_1);
Fq[1]=Fq[2]; # avoids singularity


p1 = plot(I,log10.(abs.(la)))
#plot!(Im,log10.(abs.(mom)))
p2=plot(qp,Fq,label = L"FT\left(\gamma-1\right)", linewidth=2,linestyle = :dash,legendfontsize=16,legend=:bottomright)
p2 = plot!(xi,ftq_exp,label = L"FT(\gamma_{exp}-1)", linewidth=2,xlims=(0,50),xlabel = L"\xi",
          xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
          bottom_margin=5mm, right_margin=10mm)

display(p2)          
cd("C:/Users/carlos/ownCloud/trabajo/papers/Cristobal/julia/Conductividad/dibujos")
savefig("exp3_cond_FTm1.png")

## Invert Fourier transform to obtain ga_exp, 

delq=qp[2]-qp[1];
ftq_exp2=zeros(length(qp));
ftq_exp2[1:length(ftq_exp)]=ftq_exp;
ga_exp,rp=ift_radial_fun(qp,delq,ftq_exp2);


p4= plot(r,fr,label = L"\gamma", linewidth=2,linestyle = :dash,xlims=(0,1),xlabel = L"r",
        xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
        bottom_margin=5mm, right_margin=10mm,legend=:bottomright)
plot!(rp,ga_exp.+ones(length(ga_exp)),linewidth=2,xlims=(0,1),label = L"\gamma_{exp}",
        xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
        bottom_margin=5mm, right_margin=10mm)

display(p4)
savefig("exp3_condm1.png")




