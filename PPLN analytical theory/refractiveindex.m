%% no
% a1=5.653;
% a2=0.1185;
% a3=0.2091;
% a4=89.61;
% a5=10.85;
% a6=1.97e-2;
% b1=7.941e-7;
% b2=3.134e-8;
% b3=-4.641e-9;
% b4=-2.188e-6;


%% ne
a1=5.756;
a2=0.0983;
a3=0.2020;
a4=189.32;
a5=12.52;
a6=1.32e-2;
b1=2.860e-6;
b2=4.7e-8;
b3=6.113e-8;
b4=1.516e-4;

T=50;
f=(T-24.5)*(T+570.82);
wavelength=0.5E-6:0.01E-6:4E-6;
neff=sqrt(a1+b1*f+(a2+b2*f)./((wavelength.*1E6).^2-(a3+b3*f)^2)+(a4+b4*f)./((wavelength.*1E6).^2-a5^2)-a6.*((wavelength.*1E6).^2))
figure
plot(wavelength,neff,'r');
grid on;
xlabel('wavelength(m)');
ylabel('Refractive Index');
title('5% MgO doped CLN Refractive Index');
hold on;

a1=5.653;
a2=0.1185;
a3=0.2091;
a4=89.61;
a5=10.85;
a6=1.97e-2;
b1=7.941e-7;
b2=3.134e-8;
b3=-4.641e-9;
b4=-2.188e-6;

T=50;
f=(T-24.5)*(T+570.82);
wavelength=0.5E-6:0.01E-6:4E-6;
neff=sqrt(a1+b1*f+(a2+b2*f)./((wavelength.*1E6).^2-(a3+b3*f)^2)+(a4+b4*f)./((wavelength.*1E6).^2-a5^2)-a6.*((wavelength.*1E6).^2))
plot(wavelength,neff,'b');
legend('ne 5% MgO doped CLN','no 5% MgO doped CLN');