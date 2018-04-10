%% naètení obrazu
L0a=imread('L-0-1029.tif');
L0b=imread('L-0-2030.tif');
L6a=imread('L-6-1035.tif');
L6b=imread('L-6-2036.tif');
L12a=imread('L-12-1032.tif');
L12b=imread('L-12-2033.tif');
L18a=imread('L-18-1031.tif');
L18b=imread('L-18-2034.tif');
LXa=imread('L-X-1037.tif');
LXb=imread('L-X-2038.tif');
velky=imread('velky039.tif');

%% tøídìní kanálù (red a blue)
L0aR=double(L0a(:,:,1));
L0aB=double(L0a(:,:,3));

L0bR=double(L0b(:,:,1));
L0bB=double(L0b(:,:,3));

L6aR=double(L6a(:,:,1));
L6aB=double(L6a(:,:,3));

L6bR=double(L6b(:,:,1));
L6bB=double(L6b(:,:,3));

L12aR=double(L12a(:,:,1));
L12aB=double(L12a(:,:,3));

L12bR=double(L12b(:,:,1));
L12bB=double(L12b(:,:,3));

L18aR=double(L18a(:,:,1));
L18aB=double(L18a(:,:,3));

L18bR=double(L18b(:,:,1));
L18bB=double(L18b(:,:,3));

LXaR=double(LXa(:,:,1));
LXaB=double(LXa(:,:,3));

LXbR=double(LXb(:,:,1));
LXbB=double(LXb(:,:,3));

velkyR=double(velky(:,:,1));
velkyB=double(velky(:,:,3));

%% pøevod na OD
OD={L0aR, L0aB, L0bR, L0bB, L6aR, L6aB, L6bR, L6bB, L12aR, L12aB, L12bR, L12bB, L18aR, L18aB, L18bR, L18bB, LXaR, LXaB, LXbR, LXbB, velkyR, velkyB};

k=1/65535;

for i=1:22
    OD{i}=-log10(OD{i}*k);
end

%dimenze jednotlivých matic v OD (øádky->prvky OD, sloupce->rozmìry prvkù)
s=zeros(22,2);
for i=1:22
    s(i,1)=size(OD{i},1);
end

for j=1:22
    s(j,2)=size(OD{j},2);
end

%celkový poèet namìøených hodnot ODr
z=0;
for i=1:2:21        
    z=z+s(i,1)*s(i,2);
end

%% výpoèet ROD
ROD={OD{1}, OD{3}, OD{5}, OD{7}, OD{9}, OD{11}, OD{13}, OD{15}, OD{17}, OD{19}, OD{21}};

for i=1:11
    ROD{i}=OD{2*i-1}./OD{2*i};
end

%% závislost ODr na ROD
ODrvec=[reshape(OD{1},[],1);reshape(OD{3},[],1);reshape(OD{5},[],1);reshape(OD{7},[],1);reshape(OD{9},[],1);reshape(OD{11},[],1);reshape(OD{13},[],1);reshape(OD{15},[],1);reshape(OD{17},[],1);reshape(OD{19},[],1);reshape(OD{21},[],1)];
RODvec=[reshape(ROD{1},[],1);reshape(ROD{2},[],1);reshape(ROD{3},[],1);reshape(ROD{4},[],1);reshape(ROD{5},[],1);reshape(ROD{6},[],1);reshape(ROD{7},[],1);reshape(ROD{8},[],1);reshape(ROD{9},[],1);reshape(ROD{10},[],1);reshape(ROD{11},[],1);];

fitkub=fit(RODvec,ODrvec,'poly3');

%% korekce ODr a tøídìní
ODrcor=fitkub(RODvec);

kat=zeros(6);
kat(1)=1;
for i=1:5
    kat(i+1)=kat(i)+s(4*i-3,1)*s(4*i-3,2)+s(4*i-1,1)*s(4*i-1,2);
end
kat(6)=kat(5)+s(21,1)*s(21,2);
%kat(1):(kat(2)-1)
%seskupení ODrcor podle dávek
ODc0=ODrcor(kat(1):(kat(2)-1));
ODc6=ODrcor(kat(2):(kat(3)-1));
ODc12=ODrcor(kat(3):(kat(4)-1));
ODc18=ODrcor(kat(4):(kat(5)-1));
ODcX=ODrcor(kat(5):(kat(6)-1));
%ODcV=ODrcor(kat(6):z);

%% výpoèet ref. dávek a konstrukce kalibraèní køivky
min=[0 6 12 18];
lambda=log(2)/(365*5.27);               %T(60-Co)=365*5.27 dni
Ddot=(33.13/60)*exp(-lambda*21);         %21 dni mezi stanovením pøíkonu a ozáøením
D=Ddot*min;

Dvec=zeros(kat(5)-1,1);
for i=kat(1):(kat(2)-1)     %0 min
    Dvec(i)=D(1);
end
for i=kat(2):(kat(3)-1)     %6 min
    Dvec(i)=D(2);
end
for i=kat(3):(kat(4)-1)     %12 min
    Dvec(i)=D(3);
end
for i=kat(4):(kat(5)-1)     %18 min
    Dvec(i)=D(4);
end

ODrcorfit=ODrcor(1:(kat(5)-1));
kal=fit(ODrcorfit, Dvec,'poly2');

%% støední hodnoty dávek, urèení neznámé dávky
ODmean=[mean(ODc0) mean(ODc6) mean(ODc12) mean(ODc18)];
ODXmean=mean(ODcX);
DX=kal(ODXmean);     %neznámá dávka

%% plot prvni cast
figure
plot(ODXmean,DX,'ob',ODmean,D,'x')
hold on
plot(kal)
legend('X','ref','fit')
xlabel('OD_{r,cor}')
ylabel('Dávka [Gy]')

figure
plot(fitkub)
legend('fit pro OD(ROD)')
xlabel('ROD_{rb}')
ylabel('OD_{r,cor}')

%% vypracovani pro velky film		
cely=imread('velky039cely.tif');
LVR=double(cely(:,:,1));
LVB=double(cely(:,:,3));
ODVR=-log10(LVR*k);
ODVB=-log10(LVB*k);
RODV=ODVR./ODVB;
ODcV=fitkub(RODV);
sV=size(RODV);
ODcVmat=reshape(ODcV,sV(1,1),sV(1,2));
ODVmax=max(ODcV);
ODVnorm=ODcVmat/ODVmax;

%% plot cast 2
figure
plot(ODVnorm(457,:))
xlabel('Pozice')
ylabel('OD_{norm}')

figure
plot(ODVnorm(:,789))
xlabel('Pozice')
ylabel('OD_{norm}')

figure
imagesc(ODVnorm)
colorbar

%% nacteni dat cast B
load('castB.txt');
W0=[castB(:,1) castB(:,2)];
W6=[castB(:,1) castB(:,3)];
W12=[castB(:,1) castB(:,4)];
W18=[castB(:,1) castB(:,5)];
WX=[castB(:,1) castB(:,6)];

%% ploty absorpcni spektra
figure
plot(W0(:,1),W0(:,2),'x')
hold on
plot(W12(:,1),W12(:,2),'x')
legend('D=0 Gy','D=6,58 Gy')
xlabel('Vlnová délka [nm]')
ylabel('Absorbance')

figure
plot(W0(:,1),W0(:,2))
hold on
plot(W6(:,1),W6(:,2))
plot(WX(:,1),WX(:,2))
plot(W12(:,1),W12(:,2))
plot(W18(:,1),W18(:,2))
legend('0 Gy','3,29 Gy', '3,34 Gy','6,58 Gy','9,86 Gy')
xlabel('Vlnová délka [nm]')
ylabel('Absorbance')

%% lokální extrémy
[pks0,locs0]=findpeaks(W0(:,2),W0(:,1));
[pks6,locs6]=findpeaks(W6(:,2),W6(:,1));
[pks12,locs12]=findpeaks(W12(:,2),W12(:,1));
[pks18,locs18]=findpeaks(W18(:,2),W18(:,1));
[pksX,locsX]=findpeaks(WX(:,2),WX(:,1));

figure
plot(W0(:,1),W0(:,2))
hold on
plot(locs0,pks0,'x')

figure
plot(W6(:,1),W6(:,2))
hold on
plot(locs6,pks6,'x')

figure
plot(W12(:,1),W12(:,2))
hold on
plot(locs12,pks12,'x')

figure
plot(W18(:,1),W18(:,2))
hold on
plot(locs18,pks18,'x')

figure
plot(WX(:,1),WX(:,2))
hold on
plot(locsX,pksX,'x')







