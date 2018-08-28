c = 299792458;

fname = dir('beta.dat');

for index = 1:length(fname)
    W = load(fname(index).name);
    DTE1 = W(1:1:end,:);
%     DTM1 = W(1:2:end,:);
%% parameter setting
    TE1.beta = DTE1(:,3)*1E-3;           % rad/mm
%     TM1.beta = DTM1(:,2)*1E-3;           % rad/mm

    TE1.omega = 2*pi*DTE1(:,1)*1E-15;    % rad/fs
%     TM1.omega = 2*pi*DTM1(:,1)*1E-15;    % rad/fs
    h = mean(diff(TE1.omega));
%     h = mean(diff(TM1.omega));
    TE1.wavelength = 1E9*c./DTE1(:,1);   % nm
%     TM1.wavelength = 1E9*c./DTM1(:,1);   % nm

    N = length(TE1.beta);
%     N = length(TM1.beta);
%% first derivative
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
 %% second derivative    
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

%% third derivative     
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


%% output     
    wav = TE1.wavelength;
    GV(:,index) = TE1.GV;
    GVD(:,index) = TE1.GVD;
    TOD(:,index) = TE1.TOD;
%     wav = TM1.wavelength;
%     GV(:,index) = TM1.GV;
%     GVD(:,index) = TM1.GVD;
%     TOD(:,index) = TM1.TOD;

% figure
% y1 = diff(GVD);
% x1= TE1.wavelength(1:300)
% plot(x1,y1)
end