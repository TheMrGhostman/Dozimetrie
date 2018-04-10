%%naètení obrázkù

F01 = double(imread('L-0-1029.tif'));
F0_1r = F01(:,:,1);
F0_1b = F01(:,:,3);
F02 = double(imread('L-0-2030.tif'));
F0_2r = F02(:,:,1);
F0_2b = F02(:,:,3);
F61 = double(imread('L-6-1035.tif'));
F6_1r = F61(:,:,1);
F6_1b = F61(:,:,3);
F62 = double(imread('L-6-2036.tif'));
F6_2r = F62(:,:,1);
F6_2b = F62(:,:,3);
F121 = double(imread('L-12-1032.tif'));
F12_1r = F121(:,:,1);
F12_1b = F121(:,:,3);
F122 = double(imread('L-12-2033.tif'));
F12_2r = F122(:,:,1);
F12_2b = F122(:,:,3);
F181 = double(imread('L-18-1031.tif'));
F18_1r = F181(:,:,1);
F18_1b = F181(:,:,3);
F182 = double(imread('L-18-2034.tif'));
F18_2r = F182(:,:,1);
F18_2b = F182(:,:,3);
FX1 = double(imread('L-X-1037.tif'));
FX_1r = FX1(:,:,1);
FX_1b = FX1(:,:,3);
FX2 = double(imread('L-X-2038.tif'));
FX_2r = FX2(:,:,1);
FX_2b = FX2(:,:,3);
FR = double(imread('velky039.tif'));
F_R = FR(:,:,1);
F_B = FR(:,:,3);

%% optická hustota OD
Fr = {F0_1r, F0_2r, F6_1r, F6_2r, F12_1r, F12_2r, F18_1r, F18_2r, FX_1r, FX_2r, F_R};
Fb = {F0_1b, F0_2b, F6_1b, F6_2b, F12_1b, F12_2b, F18_1b, F18_2b, FX_1b, FX_2b, F_B};

c = 65535;

DOr = cell(size(Fr));
DOb = cell(size(Fb));

for i=1:numel(Fr)%ODr
    %ODr{i} = -log(ODr{i}/c);
    %ODb{i} = -log(ODb{i}/c);
    ODr{i} = -log(Fr{i}/c);
    ODb{i} = -log(Fb{i}/c); 
end

%dimenze matic OD
size_r = zeros(11,2);
size_b = zeros(11,2);

for i=1:numel(ODr)
    [size_r(i,1), size_r(i,2)] = size(ODr{i});
    [size_b(i,1), size_b(i,2)] = size(ODb{i});
end

% total_dim_r = sum(prod(size_r,2));
% total_dim_b = sum(prod(size_b,2));

%size_r = [size_r,prod(size_r,2)];
%size_b = [size_b,prod(size_b,2)];

size_rod = prod(size_r,2);


%waitforbuttonpress;
%% výpoèet RODrb závislosti RODrb na ODr
ROD = cell(size(ODr));
for i=1:11
    ROD{i}=ODr{i}./ODb{i};
end

ODr_1D_vector = [];
ROD_1D_vector = [];

for i = 1:11
    ODr_1D_vector = vertcat(ODr_1D_vector, reshape(ODr{i},[],1));
    ROD_1D_vector = vertcat(ROD_1D_vector, reshape(ROD{i},[],1));
    
end
      
ftt = fit(ROD_1D_vector,ODr_1D_vector, 'poly3');
ftt
plot(ftt)
%% ODr,cor
ODrcor = ftt(ROD_1D_vector);

temp = cumsum(vertcat(1,size_rod));
%temp = cumsum([1;size_r(:,3)]);

shluky = [];
for i = 1:2:11
    shluky = [shluky ; temp(i)];
end
shluky(numel(shluky)) = temp(numel(temp))-2*size_rod(10);
%shluky(numel(shluky)) = temp(numel(temp))-2*size_r(10,3);


% kat=zeros(6,1);
% kat(1)=1;
% for i=1:2:11
%     kat(i+1)=kat(i)+size_r(i,1)*size_r(i,2)+size_r(i+1,1)*size_r(i+1,2);
% end
% 
% kat(6)=kat(5)+s(11,1)*s(11,2);
% 
% % kat
% %kat(1):(kat(2)-1)
% %seskupení ODrcor podle dávek
% ODc0=ODrcor(kat(1):(kat(2)-1));
% ODc6=ODrcor(kat(2):(kat(3)-1));
% ODc12=ODrcor(kat(3):(kat(4)-1));
% ODc18=ODrcor(kat(4):(kat(5)-1));
% ODcX=ODrcor(kat(5):(kat(6)-1));
% %ODcV=ODrcor(kat(6):z);

ODc0=ODrcor(shluky(1):(shluky(2)-1));
ODc6=ODrcor(shluky(2):(shluky(3)-1));
ODc12=ODrcor(shluky(3):(shluky(4)-1));
ODc18=ODrcor(shluky(4):(shluky(5)-1));
ODcX=ODrcor(shluky(5):(shluky(6)-1));

% for i=1:5
%     ODindex = ODrcor(shluky(i):(shluky(i+1)-1),1);
% end

%% vypoèet dávky pøíkonu ozaøovaèe
E_0=33.13/60;  %dávkový pøíkon za minutu
t=21*24*60;       %poèet dnù v minutách mezi kalibrací a ozáøením
T=5.27*365*24*60;   %poloèas rozpadu v minutách
E=E_0*exp(-(log(2)/T)*t);

D0=E*0;     %dávka neozáøeného vzorku
D6=E*6;     %dávka 6 min ozáøeného vzorku
D12=E*12;   %dávka 12 min ozáøeného vzorku
D18=E*18;   %dávka 18 min ozáøeného vzorku

strop = shluky(numel(shluky)-1);
Davky_vec = zeros(strop,1);

for i= 1:strop
    if(i <= (shluky(2)-1))
        Davky_vec(i) = D0;
    elseif((i >= shluky(2)) && (i <= (shluky(3)-1)))
        Davky_vec(i) = D6;
    elseif((i >= shluky(3)) && (i <= (shluky(4)-1)))
        Davky_vec(i) = D12;
    elseif((i >= shluky(4)) && (i <= (shluky(5)-1)))
        Davky_vec(i) = D18;
    end
end




%% výpoèet ref. dávek a konstrukce kalibraèní køivky
% min=[0 6 12 18];
% lambda=log(2)/(365*5.27);               %T(60-Co)=365*5.27 dni
% Ddot=(33.13/60)*exp(-lambda*21);         %21 dni mezi stanovením pøíkonu a ozáøením
% D=Ddot*min;
% 
% Dvec=zeros(kat(5)-1,1);
% for i=kat(1):(kat(2)-1)     %0 min
%     Dvec(i)=D(1);
% end
% for i=kat(2):(kat(3)-1)     %6 min
%     Dvec(i)=D(2);
% end
% for i=kat(3):(kat(4)-1)     %12 min
%     Dvec(i)=D(3);
% end
% for i=kat(4):(kat(5)-1)     %18 min
%     Dvec(i)=D(4);
% end

ODrcorfit=ODrcor(1:(shluky(5)));
kal=fit(ODrcorfit, Davky_vec,'poly2'); % prolo¾ení exponenciálou neni nic moc ('exp2'nebo 'exp')
kal

%% støední hodnoty dávek, urèení neznámé dávky
ODmean = [mean(ODc0) mean(ODc6) mean(ODc12) mean(ODc18)];

ODXmean = mean(ODcX);
DX = kal(ODXmean);     %neznámá dávka

%D_vec = [D0 ,D6 ,D12 ,D18];

%% plot prvni cast
figure
plot(ODXmean, DX,'ob' ,ODmean ,[D0 ,D6 ,D12 ,D18] ,'x')
hold on
plot(kal)
legend('X','ref','fit')
xlabel('OD_{r,cor}')
ylabel('Dávka [Gy]')

figure
plot(ftt)
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
ODcV=ftt(RODV);
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