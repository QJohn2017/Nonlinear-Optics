clear,clc
global c eps_0 t w hbar;
c =  2.99792458E8;      % speed of light
eps_0 = 8.854E-12;
lam = 2760E-9;          % m
lamSHG = lam/2;
w0 = 2*pi*c/lam;
t00 = 500E-15;
I00 = 3000;  % peak power units:W
t0 = 500E-15;
I0 = 3000;  % peak power units:W

%%%%%%%%%% SA behavior of NLM (LN crystal) %%%%%%%%%
deff = 14E-12;          % m/V
Lnl = 14E-3;             % length of NL crystal, m
Anl = pi*(2.591e-6)*(2.406e-6);       % mode area at NL crystal   W2=6um W1=5.04um
Rl = 0.75;              % linear reflectivity
n2_NL = 0;%2.1E-20;

%%%%%%%%%%% Sim parameters %%%%%%%%%%%
Nrt = 1;
Nw = 4001;
Nzfi = 300;
Nzcr= 2000;
%dw = 10*BW/(Nw-1);
%w = ( -(Nw-1)/2 : (Nw-1)/2 )*dw;
%t = linspace(-pi/dw,pi/dw,Nw);
%dt = mean(diff(t));
t = linspace(-3e-12,3e-12,Nw);
dt = t(2) - t(1);
w = 2*pi*linspace(-1/2/dt,1/2/dt,Nw);
dw = w(2)-w(1);
w = ifftshift(w);
t = ifftshift(t);
zKTP = linspace(0,Lnl,Nzcr);  % space coordinate in KTP
dzKTP = zKTP(2)-zKTP(1);

%%%%%%%%%%%% Initial pulse %%%%%%%%%%%
um=1e-6;
tem=70;
n_FF= n_MgLN(lam/um,tem);
FF=cos(0*pi/t0.*t);
AFF0 = FF.* sqrt(2*I0/Anl/(n_FF*c*eps_0)).*exp(-(sqrt(2*log(2))*t/t0).^2); %*(sech(1.22*t/t0)).^2; % units V/m
E= (Anl*n_FF*eps_0*c/2)*sum(abs(AFF0).^2)*dt  % initial energy
ASHG0 =0;

%%%%%%%%%%%%%%% GVM and GDD of PPLN%%%%%%%%%%%%
beta1_rel = -5.25e-11; %2.18/c-2.1845/c;p
%beta1_rel = -2e-12;
beta1_offset = 0;%11.4e-12;%5.25E-12;
beta2_FF =-483E-27;  
beta3_FF = 2658E-42;
beta2_SHG =116.4E-27;
beta3_SHG =315.5E-42;
D_FF = exp(-j*(beta1_offset*w + beta2_FF/2*w.^2 + beta3_FF/6*w.^3)*dzKTP/2);
D_SHG = exp(-1i*((beta1_rel + beta1_offset)*w + beta2_SHG/2*w.^2 + beta3_SHG/6*w.^3)*dzKTP/2);


supergaussian = exp(-(t./6e-12).^6);
AFF0 = AFF0 .*supergaussian;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Iff0 = (Anl*n_FF*eps_0*c/2)*abs(AFF0).^2;
AFF=AFF0;
ASHG=ASHG0; 
%%%%%%%%%%%passing PPLN
      T=1;
      Total_Energy_ff=zeros(Nzcr,1);
      Total_Energy_shg=zeros(Nzcr,1);
      w_pulse = zeros(Nzcr,1);
      L = zeros(Nzcr,1);
      fwhm_1=zeros(Nrt);
      for Z2 = zKTP
           % Propagation (1st half)
     sAFF = fft(AFF);
     sSHG = fft(ASHG);
     AFF = ifft(D_FF.*sAFF);
     ASHG = ifft(D_SHG.*sSHG);
           % nonlinear step using Runga-Kutta 4th order 
     [AFF, ASHG, n_FF, n_SHG] = dadzLN_rev(AFF,ASHG,lam,lamSHG,deff,n2_NL,Anl,Lnl,Z2,dzKTP);
           % Propagation (2st half)
     sAFF = fft(AFF);
     sSHG = fft(ASHG);
     AFF = ifft(D_FF.*sAFF);
     ASHG = ifft(D_SHG.*sSHG);   
     
     Ishgincry(1,T) = max(abs(ASHG)).^2/max(abs(AFF0)).^2;
     
     RNL=[];
     AFF_1(T,:)=AFF;
     ASHG_1(T,:)=ASHG;
     Iff_1(T,:) = (Anl*n_FF*eps_0*c/2)*abs(AFF_1(T,:)).^2;
     Ishg_1(T,:) = (Anl*n_SHG*eps_0*c/2)*abs(ASHG_1(T,:)).^2;
     Total_Energy_ff(T) = sum(abs(Iff_1(T,:))*dt);
     Total_Energy_shg(T) = sum(abs(Ishg_1(T,:))*dt);
    
     Intensity1=Iff_1;
     
     ind = T;
     AA(ind) = max(Intensity1(ind,:));
     fwhm_=find(abs(Intensity1(ind,:))>abs(max(Intensity1(ind,:))/2));
     fwhm_1(ind)=length(fwhm_)+1;
     fwhm_1(ind)=fwhm_1(ind)/(length(t)+1)*(6*1e3);   
     w_pulse(T) = fwhm_1(T);        % every pulse duration
     L(T) = T.*dzKTP;
     
     peakpower(T) = max(Iff_1(T,:)); %Total_Energy1(T)/w_pulse(T);    %peak power
     peakpower_shg(T) = max(Ishg_1(T,:)); %Total_Energy1(T)/w_pulse(T);    %peak power
     T=T+1;
      end     

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  figure(1); plot(zKTP',Ishgincry);
  Iff1(1,:) = (Anl*n_FF*eps_0*c/2)*abs(AFF).^2;
  Ishg1(1,:) = (Anl*n_SHG*eps_0*c/2)*abs(ASHG).^2;
 
fwhm=[];
RNL=[];
Intensity=Iff1;
fwhm1=zeros(Nrt);
for ind = 1:Nrt
    AA(ind) = max(Intensity(ind,:));
    fwhm=find(abs(Intensity(ind,:))>abs(max(Intensity(ind,:))/2));
    fwhm1(ind)=length(fwhm)+1;
    fwhm1(ind)=fwhm1(ind)/(length(t)+1)*(6*1e3);
    
end

fwhm1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% pulse
figure(2)
po=round(14e-3/Lnl*2000);
AFF00 = FF.*sqrt(2*I00/Anl/(n_FF*c*eps_0)).*exp(-(sqrt(2*log(2))*t/t00).^2);  % units V/m
Iff00 = (Anl*n_FF*eps_0*c/2)*abs(AFF00).^2;
plot(t*1e15,Iff00(1,:)./max(1),'k', t*1e15,Iff_1(po,:)./max(1),'r','LineWidth',2)

set(gca,'XTick',-900:300:900,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'YTick',0:1e4:4e4,'FontSize',15,'fontname','Arial')   % y axis range
xlabel('Pulse Duration(fs)')
ylabel('Power (W)')
axis([-1000 1000 0 4e4])
legend('Input','Output')

%%%%%%%%% spectrum
figure(3)
po=round(14e-3/Lnl*2000);
fre=(w+w0);
lamm=2*pi*c./fre;
plot(lamm*1e6,abs(fft(AFF_1(po,:))),'b','LineWidth',2); hold on;
fre=(w+2*w0);
lamm=2*pi*c./fre;
plot(lamm*1e6,abs(fft(ASHG_1(po,:))),'r','LineWidth',2);
set(gca,'XTick',0:1:8,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'YTick',0:0.5e10:4e10,'FontSize',15,'fontname','Arial')   % y axis range
xlabel('wavelength (\mum)')
ylabel('Intensity')
legend('FF','SHG')
axis([0 8 0 4e10]) 

%plot(zNrd,AA)
%plot(t,Ii(Nrt,:))
% yyaxis right
% plot(t,Ii(Nrt,:))
% yyaxis left
% plot(t,Ii(1,:))

%%%%%%%%%% conversion efficiency
FF0 = max(abs(Iff0(1,:)).^2);
FF1 = max(abs(Iff1(1,:)).^2);
SHG1 = max(abs(Ishg1(1,:)).^2);
SHG_conversion=(FF0-FF1)/FF0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% beam quality judgement
beam=[];
beam = Iff_1(po,:);

Total_Energy = sum(abs(beam(1,:))*dt);

RR=0;
LL=0;
for TT = 1:1:length(t)
    if beam(1,length(t)-TT+1)>(max(beam)/2)
        RR=length(t)-TT+1; 
        break
    end 
end

for TT = 1:1:length(t)       
    if beam(1,RR-TT)<=(max(beam)/2)
        LL=RR-TT;
        break
    end
end

beam_Energy= sum(abs(beam(1,LL:RR))*dt);
Q=beam_Energy/Total_Energy

ita=fwhm1*1e-15*max(Iff1)/(I0*t0)

figure(4)
subplot(1,2,1)
[AX,H1,H2] = plotyy(L*1e3,w_pulse,L*1e3,peakpower,'plot'); 
set(get(AX(1),'Ylabel'),'String','Pulse Duration(fs)'); 
set(get(AX(2),'Ylabel'),'String','Peak Power(W)'); 
set(AX(1),'FontSize',15,'fontname','Arial'); 
set(AX(2),'FontSize',15,'fontname','Arial');
set(AX(1),'yTick',[0:100:1000]); 
set(AX(2),'yTick',[0:1e4:4e4]);
xlabel('Crystal length(mm)'); 
%title('Pulse durantion and Peak power versus length of crystal'); 
set(H1,'LineWidth',3);
set(H2,'LineWIdth',3);
subplot(1,2,2)
plot(L*1e3,peakpower,L*1e3,peakpower_shg,'LineWidth',3);
set(gca,'YTick',0:0.5e4:4e4,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'XTick',0:5:20,'FontSize',15,'fontname','Arial')   % y axis range
ylabel('Peak Power(W)'); 
xlabel('Crystal length(mm)');
legend('FF','SHG')

figure(5)
Ifff_1(1001:2000,:)=Iff_1(1:1000,:);
Ifff_1(1:1000,:)=Iff_1(1001:2000,:);
imagesc((zKTP*1000),fftshift(t*1e15),fftshift(Ifff_1'))
set(gca,'YTick',-500:100:500,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'XTick',0:2:20,'FontSize',15,'fontname','Arial')   % y axis range
ylabel('Time domain(fs)')
xlabel('Crystal length(mm)')
axis([0 20 -500 500])
caxis([0 1e4])
colormap jet

figure(6)
Isshg_1(1001:2000,:)=Ishg_1(1:1000,:);
Isshg_1(1:1000,:)=Ishg_1(1001:2000,:);
imagesc((zKTP*1000),fftshift(t*1e15),fftshift(Isshg_1'))
set(gca,'YTick',-500:100:500,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'XTick',0:2:20,'FontSize',15,'fontname','Arial')   % y axis range
ylabel('Time domain(fs)')
xlabel('Crystal length(mm)')
axis([0 20 -500 500])
caxis([0 4e3])
colormap jet

figure(7)
plot(zKTP*1000,Total_Energy_ff*1e9,'b',zKTP*1000,Total_Energy_shg*1e9,'r',zKTP*1000,Total_Energy_ff*1e9+Total_Energy_shg*1e9,'k','LineWidth',2)
%set(gca,'YTick',0:100:500,'FontSize',15,'fontname','Arial')   % x axis range
set(gca,'XTick',0:2:20,'FontSize',15,'fontname','Arial')   % y axis range
ylabel('Energy(nJ)')
xlabel('Crystal length(mm)')
legend('FF engergy','SHG engergy','Total engergy')
axis([0 20 0 3])
caxis([0 1e4])
colormap jet