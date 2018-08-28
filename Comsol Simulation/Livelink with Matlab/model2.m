%% Initialization
clear all;
c = 299792458;
model = mphload('designWG');

%% parameter setting

    Y=2.4E-6;
    model.param.set('slot_y',Y);
    X=1.5E-6;
    model.param.set('slot_x',X);
 for i=1:6   
    H=(3.7+0.3i)*1E-6
    model.param.set('height',H);
%% running the comsol model
    model.study('std1').run
    param=model.result.numerical('gev1').computeResult()

%% getting the mode pattern
    figure
    mphgeom(model);
    figure
    mphplot(model,'pg1');
%% dispersion calculation--setting beta and omega
    W = param{1,1};
    TE1.beta = W(:,1)*1E-3;           % rad/mm
%     TM1.beta = DTM1(:,2)*1E-3;           % rad/mm

    TE1.omega = 2*pi*W(:,2)*1E-15;    % rad/fs
%     TM1.omega = 2*pi*DTM1(:,1)*1E-15;    % rad/fs
    h = mean(diff(TE1.omega));
%     h = mean(diff(TM1.omega));
    TE1.wavelength = 1E9*c./W(:,2);   % nm
%     TM1.wavelength = 1E9*c./DTM1(:,1);   % nm

    N = length(TE1.beta);
%     N = length(TM1.beta);
%% dispersion calculation--calculating the first derivative
    figure
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
     subplot(1,3,1);plot(TE1.wavelength,TE1.GV);xlabel('wavelength (nm)');ylabel('GV (fs/mm)');grid on;
     
%% dispersion calculation--calculating the second derivative     
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
     subplot(1,3,2);plot(TE1.wavelength,TE1.GVD);xlabel('wavelength (nm)');ylabel('GVD (fs^2/mm)');grid on;
 
%% dispersion calculation--calculating the third derivative     
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
     subplot(1,3,3);plot(TE1.wavelength,TE1.TOD);xlabel('wavelength (nm)');ylabel('TOD (fs^3/mm)');grid on;

%% getting the result--GV GVD TOD     
    wav = TE1.wavelength;
    GV(:,1) = TE1.GV;
    GVD(:,1) = TE1.GVD;
    TOD(:,1) = TE1.TOD;
%     wav = TM1.wavelength;
%     GV(:,index) = TM1.GV;
%     GVD(:,index) = TM1.GVD;
%     TOD(:,index) = TM1.TOD;

    param_H(i)= H;
    TOD_Target(i) = TOD(18);        %GVD@2760nm
end
%% parameter-variable curve plot
    figure
    plot(param_H,TOD_Target);
    xlabel('height (m)');
    ylabel('TOD (fs^3/mm)');
    title('TOD-height');
    grid on;