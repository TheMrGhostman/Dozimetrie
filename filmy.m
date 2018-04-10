% konverze minut na davku pomoci  davkoveho prikonu za 1 minutu  [Gy/min]

d=50.7/60;  % prikon k datu  14.10.2014 v Gy/min
r=datenum(2014,10,14);
m=datenum(2017,3,20);
l12=log(2)/(5.2714*365.25);

davpri=d*exp(-l12*(m-r))  % davkovy prikon  k 20.3.2017

D3=davpri*3;
D6=davpri*6;
D9=davpri*9;

% neozarene vzorky

A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\0_1.tif');
B01=scan2od(A);
A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\0_2.tif');
B02=scan2od(A);
B0=[B01;B02];
clearvars B01 B02

% 3minuty ozarene vzorky

A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\3_1.tif');
B01=scan2od(A);
A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\3_2.tif');
B02=scan2od(A);
B3=[B01;B02];
clearvars B01 B02


% 6 minut ozarene vzorky


A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\6_1.tif');
B01=scan2od(A);
A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\6_2.tif');
B02=scan2od(A);
B6=[B01;B02];
clearvars B01 B02


% 9 minut ozarene vzorky


A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\9_1.tif');
B01=scan2od(A);
A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\9_2.tif');
B02=scan2od(A);
B9=[B01;B02];
clearvars B01 B02


% neznama davka

A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\x_1.tif');
B01=scan2od(A);
A=imread('C:\Users\PC\Documents\Jaderka\od neznameho dobrodince\6 semestr\ZPRA\5\skupina_2\x_2.tif');
B02=scan2od(A);
Bn=[B01;B02];
clearvars B01 B02

% spojeni filmu do jednoho souboru pro kalibracni funkci C

C=[B0;B3;B6;B9];
[cxmax,cymax]=size(C);


% nakresli graf
plot (C(:,1),C(:,2))
%fit(C(:,1),C(:,2),'poly3')
hold on
poly3=fit(C(:,1),C(:,2),'poly3')
plot(poly3,'k');
hold off
b = gca; 
legend(b,'data','polynom 3. stupnì');
xlabel(b,'ROD_{rb} [-]');
ylabel('OD_r [-]');




%  kalibracni   x^3*p3+X^2*p2+X*p1+p0

p0=  201.7;
p1= -687.9;
p2= 781.6;
p3= -292; %-291.7

%   neozareny 

[xmax,ymax]=size(B0);
for i=1:xmax
   B0(i,2)=B0(i,1)^3*p3+B0(i,1)^2*p2+B0(i,1)*p1+p0;  
   
end

%  3 minuty 

[xmax,ymax]=size(B3);
for i=1:xmax
   B3(i,2)=B3(i,1)^3*p3+B3(i,1)^2*p2+B3(i,1)*p1+p0;  
   
end

%  6 minuty 

[xmax,ymax]=size(B6);
for i=1:xmax
   B6(i,2)=B6(i,1)^3*p3+B6(i,1)^2*p2+B6(i,1)*p1+p0;  
   
end

%  7 minuty 

[xmax,ymax]=size(B9);
for i=1:xmax
   B9(i,2)=B9(i,1)^3*p3+B9(i,1)^2*p2+B9(i,1)*p1+p0;  
   
end


%  neznamy vzorek 

[xmax,ymax]=size(Bn);
for i=1:xmax
   Bn(i,2)=Bn(i,1)^3*p3+Bn(i,1)^2*p2+Bn(i,1)*p1+p0;  
   
end


%  promena na klibracni data  (min, ODr , std)

calib=zeros(4,3);

%  naplni davku 
i=0;
for s=[0,3,6,9]
i=i+1;
  calib(i,1)=s*davpri;
end

% stredni hodnoty a  rozptyl 
calib(1,2)=mean(B0(:,2));
calib(2,2)=mean(B3(:,2));
calib(3,2)=mean(B6(:,2));
calib(4,2)=mean(B9(:,2));
calib(1,3)=std(B0(:,2));
calib(2,3)=std(B3(:,2));
calib(3,3)=std(B6(:,2));
calib(4,3)=std(B9(:,2));

figure;

plot(calib(:,2),calib(:,1),'*', 'MarkerSize',11,'Color', [0,0.35,0.33]) %pridat fir lomeny a exp
%legend('data');
xlim([2.8,3.6]);
ylim([-1,8]);
xlabel('OD_{r,cor}[-]');
ylabel('D [Gy]');
hold on
exp1 = fit(calib(:,2),calib(:,1),'exp1', 'StartPoint', [200,-5])
% plot(exp1)
zkouska = linspace(2.8,4,10000);
plot(zkouska,exp1(zkouska),'m','LineWidth',1.5)
hold on
fitpoly2=fit(calib(:,2),calib(:,1),'poly2')
plot(fitpoly2,'r')
hold on
poly3=fit(calib(:,2),calib(:,1),'poly3')
plot(poly3,'g')
hold on
rat02=fit(calib(:,2),calib(:,1),'rat02')
plot(rat02,'b')
hold off
b = gca; 
legend(b,'data','exponenciální','polynom 2. stupnì','polynom 3. stupnì','lomená funkce');
xlabel(b,'OD_{r,cor}[-]');
ylabel('D [Gy]');

x=mean(Bn(:,2))

% f=@(b,x)(b(1)*exp(b(2)*x+b(3)));
% [b_fit,R,J,CovB, MSE]=nlinfit(calib(:,2),calib(:,1), f,[10,5,10]);
% alpha = 0.32;
% betaci=nlparci(b_fit,R,'covar',CovB,'alpha',alpha);
% betaci2=nlparci(b_fit,R,'jacobian',J,'alpha',alpha);
% plot(zkouska,f(b_fit,zkouska),'y')
% hold off

davka= 1.436e+08*exp(-5.91*x)
 a =   1.436e+08 
 bb= -5.91
 CH1=exp(bb*x)
 CH2=exp(a*exp(bb*x))
 chyba=sqrt(((2.6*10^8)*CH1)^2+((0.6)*CH2)^2)