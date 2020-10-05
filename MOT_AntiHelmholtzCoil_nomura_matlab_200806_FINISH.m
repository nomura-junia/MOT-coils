%**********************************************************************
%  MOT_AntiHelmholtzCoil_nomura_calculation ver.1.00
%**********************************************************************
%
%***********************************************************************

clear all
clf reset

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cvac = 2.99792458e8;            %speed of light in vacuum free space 
mu0  = 1.257*1.0e-6;            %permeability of free space
ep0  = 1.0/(cvac*cvac*mu0);     %permittivity of free space（C/Vm）
e=1.602176634*1.0e-19;          %electron charge
kb=1.38*1.0e-23;                %boltzman
h=6.63*1.0e-23;                 %plank
h_bar=h/(2.0*pi);               %plank/2Pi

%***********************************************************************
%     Coil parameters
%***********************************************************************
dia_in=1.2*1.0e-3;
dia_out=1.2*1.0e-3+0.11*1.0e-3;
T_damega=200+273.15;
T_est=60+273.15;

%***********************************************************************
%     Geometric parameters
%***********************************************************************
n_yoko=16;
n_tate=28;
I=1.0;
rad_0=65*10^-3;     %コイル半径（スタート位置）（m）45*10^-3;
dis_0=45*10^-3;     %コイル中心-磁場中心間距離（m）65*10^-3;
wire_diam=1.2*10^-3;
wire_outer=1.31*10^-3;

%***********************************************************************
%       Calculation
%***********************************************************************

% rad=65*10^-3;     %コイル半径（スタート位置）（m）
% dis=45*10^-3;     %コイル中心-磁場中心間距離（m）
Length=0;
z_MM=zeros(100);
Bz_MM=zeros(100);
rou_MM=zeros(100);
Brou_MM=zeros(100);
% syms z f(z) f1(z) f2(z)
% rou=0;
% k1=4*rad*rou/((rad+rou)^2+((z)-dis)^2);            %楕円積分の変数
% [K,E] = ellipke(k1);  %Kは第一種完全楕円積分、Eは第二種完全楕円積分を返す.
% f1(z)=mu0*I/(2*pi*sqrt((rad+rou)^2+((z)-dis)^2))*[K+E*(rad^2-rou^2-((z)-dis)^2)/(rad^2-rou^2+((z)-dis)^2)];
% f2(z)=mu0*(-I)/(2*pi*sqrt((rad+rou)^2+((z)-(-dis))^2))*[K+E*(rad^2-rou^2-((z)-(-dis))^2)/(rad^2-rou^2+((z)-(-dis))^2)];
% f(z)=f1(z)+f2(z);
% fplot(f1*10^4, [-0.1 0.1])
% fplot(f2*10^4, [-0.1 0.1])
% fplot(f*10^4, [-0.01 0.01])
% 
% syms rou g(rou) g1(rou) g2(rou) 
% z=0;
% k1=4*rad*rou/((rad+rou)^2+((z)-dis)^2);            %楕円積分の変数
% [K,E] = ellipke(k1);  %Kは第一種完全楕円積分、Eは第二種完全楕円積分を返す.
% g1(rou)=mu0*I/(2*pi*(rou))*(z-dis)/sqrt((rad+(rou))^2+((z)-dis)^2)*[-K+E*(rad^2+(rou)^2+((z)-dis)^2)/((rad-(rou))^2+((z)-dis)^2)];
% g2(rou)=mu0*(-I)/(2*pi*(rou))*(z-(-dis))/sqrt((rad+(rou))^2+((z)-(-dis))^2)*[-K+E*(rad^2+(rou)^2+((z)-(-dis))^2)/((rad-(rou))^2+((z)-(-dis))^2)];
% %後ろのE*(rad^2+(rou)^2+…のところが、E*(rad^2+(rou)^2-…かも？（…前の符号）
% g(rou)=g1(rou)+g2(rou);
% fplot(g1*10^4, [-0.1 0.1])
% fplot(g2*10^4, [-0.1 0.1])
% fplot(g*10^4, [-0.1 0.1])

%%B_軸方向%%
z=-0.1;
f=0;
dB_dz=0.0;

for i=1:1:201
for yoko=1:1:n_yoko
    for tate=1:1:n_tate
        rad=rad_0+(yoko-1)*wire_outer;
        dis=dis_0+(tate-1)*wire_outer;
        rou=0;
        %Upper Coil
        k1=4*rad*rou/((rad+rou)^2+((z)-dis)^2);            %楕円積分の変数
        [K,E] = ellipke(k1);  %Kは第一種完全楕円積分、Eは第二種完全楕円積分を返す.
        f1=mu0*I/(2*pi*sqrt((rad+rou)^2+((z)-dis)^2))*[K+E*(rad^2-rou^2-((z)-dis)^2)/(rad^2-rou^2+((z)-dis)^2)];
        %Lower Coil       
        k1=4*rad*rou/((rad+rou)^2+((z)-(-dis))^2);            %楕円積分の変数
        [K,E] = ellipke(k1);  %Kは第一種完全楕円積分、Eは第二種完全楕円積分を返す.
        f2=mu0*(-I)/(2*pi*sqrt((rad+rou)^2+((z)-(-dis))^2))*[K+E*(rad^2-rou^2-((z)-(-dis))^2)/(rad^2-rou^2+((z)-(-dis))^2)];
        %Coilsの足し合わせ
        f=f1+f2;
        Bz_M(tate,yoko)=f; %各コイル１巻きが作る磁場の総和
    end

end
Bz=sum(Bz_M,'all');
Bz_MM(i,1)=z;
Bz_MM(i,2)=Bz;    
if z==0.010;
    dB_dz=Bz/(z*10^2) %G/cm
end
z=z+1*10^-3; 
end


figure(1)
plot(Bz_MM(:,1), Bz_MM(:,2))
xlim([-0.1 0.1])

figure(2)
plot(Bz_MM(:,1), Bz_MM(:,2))
xlim([-0.01 0.01])

dB_dz_axial=dB_dz;

%%B_動径方向%%
rou=0.0;
for i=1:1:100
for yoko=1:1:n_yoko
    for tate=1:1:n_tate
        rad=rad_0+(yoko-1)*wire_outer;
        dis=dis_0+(tate-1)*wire_outer;
        z=0.0;
        k2=4*rad*rou/((rad+rou)^2+(z-dis)^2);            %楕円積分の変数
        [K,E] = ellipke(k2);  %Kは第一種完全楕円積分、Eは第二種完全楕円積分を返す.
        g1=mu0*I/(2*pi*(rou))*(z-dis)/sqrt((rad+(rou))^2+((z)-dis)^2)*[-K+E*(rad^2+(rou)^2+((z)-dis)^2)/((rad-(rou))^2+((z)-dis)^2)];
        g2=mu0*(-I)/(2*pi*(rou))*(z-(-dis))/sqrt((rad+(rou))^2+((z)-(-dis))^2)*[-K+E*(rad^2+(rou)^2+((z)-(-dis))^2)/((rad-(rou))^2+((z)-(-dis))^2)];
        g=g1+g2;
        Brou_M(tate,yoko)=g;
    end
end
Brou=sum(Brou_M,'all');
rou_MM(i,1)=rou;
Brou_MM(i,2)=Brou;
 
rou=rou+0.001;
end
figure(3)
plot(rou_MM(:,1), Brou_MM(:,2))


%subplot(3,1,1)
%plot(Bz_M)
%xlim([0.0 150.0])
%xlabel('Length (m)')
%ylabel('SN ratio')

%%行列の設定%%

%i_max=100;
%for i=1:1:i_max
%    B_z=:
%    B_rou=;
%    I_current=I_current+0.1;
%end
