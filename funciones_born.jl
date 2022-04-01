module funciones_born
export trapz, ift_radial_fun, ft_radial_fun, general_eig_fun, general_mom_fun, general_mom_fun_ne
using ArbNumerics, SpecialFunctions, FFTW

function trapz(r,f)
    N=length(f)
    int=(f[1]+f[2])/2*(r[2]-r[1])
    for jj = 2:(N-1)
        int = int + (f[jj]+f[jj+1])/2*(r[jj+1]-r[jj]) 
    end
    return int
end

##########################################

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

###################    

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
######################

 function trapz(r,f)
    N=length(f)
    int=(f[1]+f[2])/2*(r[2]-r[1])
    for jj = 2:(N-1)
        int = int + (f[jj]+f[jj+1])/2*(r[jj+1]-r[jj]) 
    end
    return int
end

########### Funciones Juan

# eigenvalues associated with a potential with n steps
function general_eig_fun(rk,q,N_e)
    m=length(q);
    rk=ArbFloat.(rk)
    #We define the  vector that  they have as its components the numbers  sqrt(abs(-gm(j))) 
    alpha=ArbFloat.(abs.(q).^(1/2));
    
    #We define the spherival Bessel  functions j_l(k*r), and y_l(k*r).
    #k and r are mute variables, they are not indicating the vectorswe have
    #define before
    j= (order,k,r) -> sqrt(pi/(2*k*r))*besselj(ArbFloat(order+1/2),ArbFloat(k*r));
    y= (order,k,r) -> sqrt(pi/(2*k*r))*bessely(ArbFloat(order+1/2),ArbFloat(k*r));
    j0= (order,k,r)-> r.^order;
    bi= (order,k,r)-> sqrt(pi/(2*k*r))*besseli(ArbFloat(order+1/2),ArbFloat(k*r));
    bk=(order,k,r)-> (-1)^(order)*sqrt(pi/(2*k*r))*besselk(ArbFloat(order+1/2),ArbFloat(k*r));
    y0=(order,k,r)-> r.^(-(order+1));

    jp=(order,k,r)-> k*(j(order-1,k,r)-(order+1)./(r*k).*j(order,k,r));
    yp=(order,k,r)-> k*(y(order-1,k,r)-(order+1)./(r*k).*y(order,k,r));
    j0p=(order,k,r)-> order*r.^(order-1);
    bip=(order,k,r)-> k*(bi(order-1,k,r)-(order+1)./(r*k).*bi(order,k,r));
    bkp=(order,k,r)-> k*(bk(order-1,k,r)-(order+1)./(r*k).*bk(order,k,r));
    y0p=(order,k,r) -> -(order+1)*r.^(-(order+2));

    #We solve the system according to the theorical scheme we have developed
    I=collect(0:(N_e-1));
    lambda = ArbFloat.(zeros(length(I)))
    for l=0:N_e-1
        A=1;
        B=0;
        X=[A;B];
        for ss=1:m-1
        
            if q[ss]<0
                z=j(l,ArbFloat(alpha[ss]),ArbFloat(rk[ss]));
                w=y(l,alpha[ss],rk[ss]);
                zp=jp(l,alpha[ss],rk[ss]);
                wp=yp(l,alpha[ss],rk[ss]);
            elseif q[ss]>0
                z=bi(l,alpha[ss],rk[ss]);
                w=bk(l,alpha[ss],rk[ss]);
                zp=bip(l,alpha[ss],rk[ss]);
                wp=bkp(l,alpha[ss],rk[ss]);
            else
                z=j0(l,alpha[ss],rk[ss]);
                w=y0(l,alpha[ss],rk[ss]);
                zp=j0p(l,alpha[ss],rk[ss]);
                wp=y0p(l,alpha[ss],rk[ss]); 
            end
            MM=[z w;zp wp];
            X=MM*X;
            if q[ss+1]<0
                z=j(l,alpha[ss+1],rk[ss]);
                w=y(l,alpha[ss+1],rk[ss]);
                zp=jp(l,alpha[ss+1],rk[ss]);
                wp=yp(l,alpha[ss+1],rk[ss]);
            elseif q[ss+1]>0
                z=bi(l,alpha[ss+1],rk[ss]);
                w=bk(l,alpha[ss+1],rk[ss]);
                zp=bip(l,alpha[ss+1],rk[ss]);
                wp=bkp(l,alpha[ss+1],rk[ss]);
            else
                z=j0(l,alpha[ss+1],rk[ss]);
                w=y0(l,alpha[ss+1],rk[ss]);
                zp=j0p(l,alpha[ss+1],rk[ss]);
                wp=y0p(l,alpha[ss+1],rk[ss]); 
            end
            MM=1/(z*wp-zp*w)*[wp -w; -zp z];
            X=MM*X;
        end

        if q[m]<0
            z=j(l,alpha[m],rk[m]);
            w=y(l,alpha[m],rk[m]);
            zp=jp(l,alpha[m],rk[m]);
            wp=yp(l,alpha[m],rk[m]);
        elseif q[m]>0
            z=bi(l,alpha[m],rk[m]);
            w=bk(l,alpha[m],rk[m]);
            zp=bip(l,alpha[m],rk[m]);
            wp=bkp(l,alpha[m],rk[m]);
        else
            z=j0(l,alpha[m],rk[m]);
            w=y0(l,alpha[m],rk[m]);
            zp=j0p(l,alpha[m],rk[m]);
            wp=y0p(l,alpha[m],rk[m]); 
        end
        MM=[z w;zp wp];
        X=MM*X;
        lambda[l+1]=X[2]/X[1];
        #display(typeof(lambda[l+1]))
    end

    lambda=lambda .- I;
    
    return I,lambda
end

# momenta associated with a potential with n steps
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

end