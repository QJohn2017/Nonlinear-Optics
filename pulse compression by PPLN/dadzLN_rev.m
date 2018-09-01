function   [AFF, ASHG, n_FF, n_SHG] = dadzLN_rev(AFF,ASHG,lam,lamSHG,deff,n2_NL,Anl,Lnl,Z2,dzKTP)

global c eps_0 t w hbar;

um=1e-6;
w_FF = 2*pi*c/lam;
w_SHG = 2*pi*c/lamSHG;
tem=70;
% kapa=1E1;

n_FF= n_MgLN(lam/um,tem);
dk = -28e3;   %k_SHG - 2*k_FF-kG;
%[dk,phase_mis] = cascaded_quadratic_nonlinearity(Lnl,Z2,dzKTP);
n_SHG = lamSHG/2/pi*(dk+2*2*pi*n_FF/lam);
n2_eff=0;
-(4*pi/(c*eps_0))*(Lnl/lam)*(deff^2/(n_SHG*n_FF*n_FF))/(dk*Lnl);
        %%% efficiency prediction
        LL=sqrt(n_FF*n_SHG)*c/(2*w_FF*deff*max(abs(AFF)));
        
        f_FF = (conj(AFF).*ASHG)*exp(-j*dk*Z2);
        g_FF = Anl*(n_FF*eps_0*c/2)*abs(AFF).^2 .*AFF;
        %w_FF*deff/n_FF/c/Anl;
        kFF_1 = (-j*w_FF*deff/n_FF/c)*1*f_FF + (-j*w_FF*(n2_NL+n2_eff)/c/Anl)* g_FF; %
        f_SHG = (AFF.^2)*exp(j*dk*Z2);
        g_SHG = Anl*(n_SHG*eps_0*c/2)*abs(ASHG).^2 .*ASHG;
        kSHG_1 = (-j*w_SHG*deff/n_SHG/c/2)*f_SHG + (-j*w_SHG*(n2_NL+n2_eff)/c/Anl)* g_SHG ;  %w_SHG*deff/n_SHG/c/2
        clear f_FF f_SHG g_FF g_SHG;

        %A_half2 = A + k1*dz/2;
        AFF_half2 = AFF + kFF_1*dzKTP/2;
        ASHG_half2 = ASHG + kSHG_1*dzKTP/2;

        %k2 = dAdz(A_half2,gamma,w0,w); 
        f_FF = (conj(AFF_half2).*ASHG_half2)*exp(-j*dk*(Z2+dzKTP/2));
        g_FF = Anl*(n_FF*eps_0*c/2)*abs(AFF_half2).^2 .*AFF_half2;
        kFF_2 = (-j*w_FF*deff/n_FF/c)*1* f_FF  + (-j*w_FF*(n2_NL+n2_eff)/c/Anl)* g_FF ;
        f_SHG = (AFF_half2.^2)*exp(j*dk*(Z2+dzKTP/2));
        g_SHG = Anl*(n_SHG*eps_0*c/2)*abs(ASHG_half2).^2 .*ASHG_half2;
        kSHG_2 = (-j*w_SHG*deff/n_SHG/c/2)* f_SHG  + (-j*w_SHG*(n2_NL+n2_eff)/c/Anl)* g_SHG;
        clear f_FF f_SHG g_FF g_SHG;
        clear AFF_half2 ASHG_half2;
   
        %A_half3 = A + k2*dz/2;
        AFF_half3 = AFF + kFF_2*dzKTP/2;
        ASHG_half3 = ASHG + kSHG_2*dzKTP/2;
                
        %k3 = dAdz(A_half3,gamma,w0,w); 
        f_FF = (conj(AFF_half3).*ASHG_half3)*exp(-j*dk*(Z2+dzKTP/2));
        g_FF = Anl*(n_FF*eps_0*c/2)*abs(AFF_half3).^2 .*AFF_half3;
        kFF_3 = (-j*w_FF*deff/n_FF/c)*1* f_FF  + (-j*w_FF*(n2_NL+n2_eff)/c/Anl)* g_FF ;
        f_SHG = (AFF_half3.^2)*exp(j*dk*(Z2+dzKTP/2));
        g_SHG = Anl*(n_SHG*eps_0*c/2)*abs(ASHG_half3).^2 .*ASHG_half3;
        kSHG_3 = (-j*w_SHG*deff/n_SHG/c/2)* f_SHG + (-j*w_SHG*(n2_NL+n2_eff)/c/Anl)* g_SHG ;
        clear f_FF f_SHG g_FF g_SHG;
        clear AFF_half3 ASHG_half3;
                
        %A_full = A + k3*dz;
        AFF_full = AFF + kFF_3*dzKTP;
        ASHG_full = ASHG + kSHG_3*dzKTP;
        
        %k4 = dAdz(A_full,gamma,w0,w); 
        f_FF = (conj(AFF_full).*ASHG_full)*exp(-j*dk*(Z2+dzKTP));
        g_FF = Anl*(n_FF*eps_0*c/2)*abs(AFF_full).^2 .*AFF_full;
        kFF_4 = (-j*w_FF*deff/n_FF/c)*1* f_FF + (-j*w_FF*(n2_NL+n2_eff)/c/Anl)* g_FF ;
        f_SHG = (AFF_full.^2)*exp(j*dk*(Z2+dzKTP));
        g_SHG = Anl*(n_SHG*eps_0*c/2)*abs(ASHG_full).^2 .*ASHG_full;
        kSHG_4 = (-j*w_SHG*deff/n_SHG/c/2)* f_SHG + (-j*w_SHG*(n2_NL+n2_eff)/c/Anl)* g_SHG ;
        clear f_FF f_SHG g_FF g_SHG;
        clear AFF_full ASHG_full;
        
        %A = A + dz*(k1 + 2*k2 + 2*k3 + k4)/6;
        AFF = AFF + dzKTP*(kFF_1 + 2*kFF_2 + 2*kFF_3 + kFF_4)/6;
        ASHG = ASHG + dzKTP*(kSHG_1 + 2*kSHG_2 + 2*kSHG_3 + kSHG_4)/6;




end