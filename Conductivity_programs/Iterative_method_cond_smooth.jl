# Iterative method to obtain a conductivity from the DtN map (Calderon pb)
## functions are defined at the beginning
# The following packages must be installed
using ArbNumerics, Flux, SpecialFunctions, FFTW
using FileIO, JLD2, Plots, Measures, LaTeXStrings
setextrabits(0)
setworkingprecision(ArbFloat, bits = 2048) # digits = 
###################
###### Functions ############

function ft_radial_fun(r,h,Fr)

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
    
    ########### 

    function trapz(r,f)
        N=length(f)
        if N==1
            int = 0
        elseif N==2
        int=(f[1]+f[2])/2*(r[2]-r[1])
        else
        int=(f[1]+f[2])/2*(r[2]-r[1])    
        for jj = 2:(N-1)
            int = int + (f[jj]+f[jj+1])/2*(r[jj+1]-r[jj]) 
        end
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
    
    function general_mom_fun_ga(r,gm,m)
    
        n=length(r);
        mom=ArbFloat.(zeros(m))
        r=ArbFloat.(r)
        s_gm=sqrt.(ArbFloat.(gm))
        for kk=0:m-1
            mom[kk+1]=log(s_gm[1])*r[1]^(2*kk+1);
            for jj=2:n
                mom[kk+1] += log(s_gm[jj])*(r[jj]^(2*kk+1)-r[jj-1]^(2*kk+1));
            end
            mom[kk+1]=mom[kk+1]*2*kk
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
#############################

## Define the conductivity
h=0.0001;
r=0:h:1; # radial variable
d=0.9;
Psi=exp.(-(r.^2)./(((r.^2).-d^2).^2)).*(any.(x->(x<d),r))
alpha=0.3;
fg=(alpha*Psi.+1).^2;

## Compute ga(r) piecewise constant approximation: rj, ga
N_s=10000;        # number of steps
#N_s=5;
rj=(1:N_s)/N_s; # discrete radius
ga=zeros(N_s)
valu=[findmin(abs.(r .- rj[jj]))[2] for jj in 1:N_s]
ga = fg[valu]

## Compute eigenvalues
N_e=400;
I,la=general_eig_fun(rj,ga,N_e) # lambda_k-k
#

## Compute t_exp
xi_max=160;
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
fr_1[1:length(fr)]=fr-ones(length(fr));
Fq,qp=ft_radial_fun(r_1,h,fr_1);
Fq[1]=Fq[2]; # avoids singularity

p2 = plot(xi,ftq_exp)
plot!(qp,Fq,xlims=(0,50))

## Invert Fourier transform to obtain ga_exp, 
delq=qp[2]-qp[1];
ftq_exp2=zeros(length(qp));
ftq_exp2[1:length(ftq_exp)]=ftq_exp;
ga_exp_m1,rp=ift_radial_fun(qp,delq,ftq_exp2);

# iterative point fixed for gamma
function iter_fun(q_1,valu,rj,N_e,N_m,xi,qp,delq,ftq_exp)
    # 1. transform q_1 in a step function
    q1_step=q_1[valu]
    # 2. Compute eigenvalues and momenta
    I1,la_q1=general_eig_fun(rj,q1_step .+ 1,N_e+1) # lambda_k-k Necesito el k+1!!
    Im1,mom_q1=general_mom_fun(rj,q1_step,N_e)
    # 3. Modify eigenvalues to compute correction
    lam_q1=la_q1 ./ (ArbFloat(2)*(I1) .* (I1 .+ ArbFloat(1/2)))
    # 4. Compute t_2-t_exp Arb_precision
    ftq_2 = zeros(length(xi))*ArbFloat(0)
    for k=1:N_e      
        ftq_2 += (-1)^(k-1)*(xi./2).^(2*(k-1))./(gamma(ArbFloat(k))*(gamma(ArbFloat(k+1/2)))).*(mom_q1[k]-lam_q1[k+1]); 
    end
    ftq_2=2*(ArbFloat(pi))^(3/2)*ftq_2;
    # 5. compute q_2
    ftq_exp22 = zeros(length(qp))
    ftq_exp22[1:length(ftq_exp)] = ftq_2 .+ ftq_exp
    q_2,rp=ift_radial_fun(qp,delq,ftq_exp22);
    return q_2,ftq_exp22
end

# Plot results
p4 = plot(r,fg,label = L"\gamma",linewidth=2,
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     xlabel = "r", bottom_margin=5mm, right_margin=10mm)
p4 = plot!(rp,ga_exp_m1 .+ 1,xlims=(0,1),label = L"\gamma^{exp}",legendfontsize=16,linewidth=2)
ga_1 = ga_exp_m1

################
p4= plot(r,fg,label = L"\gamma",linewidth=2,linestyle = :dash,
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     xlabel = "r", bottom_margin=5mm, right_margin=10mm)
 plot!(rp,ga_exp_m1 .+ 1,xlims=(0,1),label = L"\gamma^{exp}",legendfontsize=16,linewidth=2)

n_iter= 6
err_lin=zeros(n_iter+1)
err_l2=zeros(n_iter+1)

err_l2[1] = sqrt(h * sum((ga_exp_m1[1:length(r)] .+ 1 -fg).^2))
err_lin[1],ind_i = findmax(abs.(ga_exp_m1[1:length(r)] .+1 -fg))

ga_1 = ga_exp_m1
for kk=1:n_iter
    global ga_1
    ga_1,ftq_exp22=iter_fun(ga_1,valu,rj,N_e,N_e,xi,qp,delq,ftq_exp)
    if kk <=2
      p4=plot!(rp,ga_1 .+ 1,xlims=(0,1),label = "\$\\gamma^$kk \$",linewidth=2,legendfontsize=16)
    end
    err_l2[kk+1] = sqrt(h * sum((ga_1[1:length(r)] .+ 1 -fr).^2))
    err_lin[kk+1],ind_i = findmax(abs.(ga_1[1:length(r)] .+ 1 -fr))
end
display(p4)
savefig("ejem_ga_iter_s.png")

# Save results                                                                                
FileIO.save("results_iter_s.jld2","err_l2",err_l2,"err_lin",err_lin)