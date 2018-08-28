%% parameter setting
    clear all;
    c = 299792458;
    tem = 50;
    freq = 0.8e14:0.1e12:3.2e14;
    um = 1e-6;
    wavelength = c./freq;
    n=n_MgLN(wavelength/um,tem);

%% propagation parameters

    TE1.beta = n.*(2*pi)./wavelength*1e-3;  % rad/mm
    TE1.beta = TE1.beta';
    TE1.omega = 2*pi.*freq*1e-15;            % rad/fs
    TE1.omega = TE1.omega';
    TE1.wavelength = 1e9.*wavelength;        %nm
    TE1.wavelength = TE1.wavelength';
    N = length ( TE1.beta);
    h = mean(diff(TE1.omega));
    
%% GV,GVD,TOD calculation
    FD = zeros(N);  % first derivative
    for ind = 3 : N-2
        FD(ind,ind-2) = 1/12;
        FD(ind,ind-1) = -2/3;
        FD(ind,ind+1) = 2/3;
        FD(ind,ind+2) = -1/12;
    end
    for ind = 1 : 2
        FD(ind,ind) = -25/12;
        FD(ind,ind+1) = 4;
        FD(ind,ind+2) = -3;
        FD(ind,ind+3) = 4/3;
        FD(ind,ind+4) = -1/4;
    end
    for ind = N-1 : N
        FD(ind,ind) = 25/12;
        FD(ind,ind-1) = -4;
        FD(ind,ind-2) = 3;
        FD(ind,ind-3) = -4/3;
        FD(ind,ind-4) = 1/4;
    end
    TE1.GV = FD*TE1.beta/h;
%     TM1.GV = FD*TM1.beta/h;
    subplot(1,3,1);plot(TE1.wavelength,TE1.GV);xlabel('wavelength (nm)');ylabel('GV (fs/mm)');
    grid on;
    
    SD = zeros(N);  % second derivative
    for ind = 3 : N-2
        SD(ind,ind-2) = -1/12;
        SD(ind,ind-1) = 4/3;
        SD(ind,ind)   = -5/2;
        SD(ind,ind+1) = 4/3;
        SD(ind,ind+2) = -1/12;
    end
    for ind = 1 : 2
        SD(ind,ind) = 15/4;
        SD(ind,ind+1) = -77/6;
        SD(ind,ind+2) = 107/6;
        SD(ind,ind+3) = -13;
        SD(ind,ind+4) = 61/12;
        SD(ind,ind+5) = -5/6;    
    end
    for ind = N-1 : N
        SD(ind,ind) = 15/4;
        SD(ind,ind-1) = -77/6;
        SD(ind,ind-2) = 107/6;
        SD(ind,ind-3) = -13;
        SD(ind,ind-4) = 61/12;
        SD(ind,ind-5) = -5/6;
    end
    TE1.GVD = SD*TE1.beta/h/h;
%     TM1.GVD = SD*TM1.beta/h/h;
    subplot(1,3,2);plot(TE1.wavelength,TE1.GVD);xlabel('wavelength (nm)');ylabel('GVD (fs^2/mm)');
    grid on;
    
    TD = zeros(N);  % third derivative
    for ind = 4 : N-3
        TD(ind,ind-3) = 1/8;
        TD(ind,ind-2) = -1;
        TD(ind,ind-1) = 13/8;
        TD(ind,ind+1) = -13/8;
        TD(ind,ind+2) = 1;
        TD(ind,ind+3) = -1/8;
    end
    for ind = 1 : 3
        TD(ind,ind) = -49/8;
        TD(ind,ind+1) = 29;
        TD(ind,ind+2) = -461/8;
        TD(ind,ind+3) = 62;
        TD(ind,ind+4) = -307/8;
        TD(ind,ind+5) = 13;
        TD(ind,ind+6) = -15/8;
    end
    for ind = N-2 : N
        TD(ind,ind) = 49/8;
        TD(ind,ind-1) = -29;
        TD(ind,ind-2) = 461/8;
        TD(ind,ind-3) = -62;
        TD(ind,ind-4) = 307/8;
        TD(ind,ind-5) = -13;
        TD(ind,ind-6) = 15/8;
    end
    TE1.TOD = TD*TE1.beta/h/h/h;
%     TM1.TOD = TD*TM1.beta/h/h/h;
     subplot(1,3,3);plot(TE1.wavelength,TE1.TOD);xlabel('wavelength (nm)');ylabel('TOD (fs^3/mm)');
     grid on;
     
    wav = TE1.wavelength;
    GV(:,1) = TE1.GV;
    GVD(:,1) = TE1.GVD;
    TOD(:,1) = TE1.TOD;