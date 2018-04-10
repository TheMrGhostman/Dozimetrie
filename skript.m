clc, close all
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
Fr = {F0_1r, F0_2r, F6_1r, F6_2r, F12_1r, F12_2r, F18_1r, F18_2r, FX_1r, FX_2r};%, F_R};
Fb = {F0_1b, F0_2b, F6_1b, F6_2b, F12_1b, F12_2b, F18_1b, F18_2b, FX_1b, FX_2b};%, F_B};

c = 65535;

ODr = cell(size(Fr));
ODb = cell(size(Fb));

for i=1:numel(Fr)
    ODr{i} = -log(Fr{i}/c);
    ODb{i} = -log(Fb{i}/c); 

end
%dimenze matic OD
size_r = zeros(numel(ODr),2);
size_b = zeros(numel(ODr),2);

for i=1:numel(ODr)
    [size_r(i,1), size_r(i,2)] = size(ODr{i});
    %[size_b(i,1), size_b(i,2)] = size(ODb{i});
end

size_rod = prod(size_r,2);
%% výpoèet RODrb závislosti RODrb na ODr
ROD = cell(size(ODr));

for i=1:numel(ODr)
    ROD{i}=ODr{i}./ODb{i};

end

ODr_1D_vector = [];
ROD_1D_vector = [];

for i = 1:numel(ROD)
    ODr_1D_vector = vertcat(ODr_1D_vector, reshape(ODr{i},[],1));
    ROD_1D_vector = vertcat(ROD_1D_vector, reshape(ROD{i},[],1));
    
end
      
corel_fit = fit(ROD_1D_vector,ODr_1D_vector, 'poly3');
corel_fit
%% ODr,cor
OD_r_corel = corel_fit(ROD_1D_vector);

temp = cumsum(vertcat(1,size_rod));

indexy_davek = [];
for i = 1:2:numel(ROD)+1
    indexy_davek = vertcat(indexy_davek, temp(i));
    %skáèu od dvì proto¾e potøebuju zapoèítat F0_1 i F0_2 dohromady, atd
end

%ODRrcor je vektor nafitovaných hodnot, následujících 6 výpoètù je k
%pøiøazení nafitovaných hodnot hustot k dávkám

ODc0 = OD_r_corel(indexy_davek(1):(indexy_davek(2)-1));
ODc6 = OD_r_corel(indexy_davek(2):(indexy_davek(3)-1));
ODc12 = OD_r_corel(indexy_davek(3):(indexy_davek(4)-1));
ODc18 = OD_r_corel(indexy_davek(4):(indexy_davek(5)-1));
ODcX = OD_r_corel(indexy_davek(5):(indexy_davek(6)-1));

%% vypoèet dávky pøíkonu ozaøovaèe
E_0=33.13/60;  %dávkový pøíkon za minutu
t=21*24*60;       %poèet dnù v minutách mezi kalibrací a ozáøením
T=5.27*365*24*60;   %poloèas rozpadu v minutách
E=E_0*exp(-(log(2)/T)*t);

D0=E*0;     %dávka neozáøeného vzorku
D6=E*6;     %dávka 6 min ozáøeného vzorku
D12=E*12;   %dávka 12 min ozáøeného vzorku
D18=E*18;   %dávka 18 min ozáøeného vzorku

strop = indexy_davek(numel(indexy_davek)-1);
Davky_vec = zeros(strop,1);

for i= 1:strop
    if(i <= (indexy_davek(2)-1))
        Davky_vec(i) = D0;
    elseif((i >= indexy_davek(2)) && (i <= (indexy_davek(3)-1)))
        Davky_vec(i) = D6;
    elseif((i >= indexy_davek(3)) && (i <= (indexy_davek(4)-1)))
        Davky_vec(i) = D12;
    elseif((i >= indexy_davek(4)) && (i <= (indexy_davek(5)-1)))
        Davky_vec(i) = D18;
    end
end



%%
%OD_r_corelfit = OD_r_corel(1:(indexy_davek(5)));
kalibrace_fit = fit(OD_r_corel(1:(indexy_davek(5))), Davky_vec,'poly2'); % prolo¾ení exponenciálou neni nic moc ('exp2'nebo 'exp')
kalibrace_fit

%% støední hodnoty dávek, urèení neznámé dávky
ODmean = [mean(ODc0) mean(ODc6) mean(ODc12) mean(ODc18)];

DX = kalibrace_fit(mean(ODcX));     %neznámá dávka

%% plot prvni cast
figure
plot(mean(ODcX), DX,'ob' ,ODmean ,[D0 ,D6 ,D12 ,D18] ,'x')
hold on
plot(kalibrace_fit)
legend('X','ref','fit')
xlabel('OD_{r,cor}')
ylabel('Dávka [Gy]')

figure
hold on
plot(ROD_1D_vector,ODr_1D_vector,'x')
plot(corel_fit)
legend('fit pro OD(ROD)')
xlabel('ROD_{rb}')
ylabel('OD_{r,cor}')

%% vypracovani pro velky film	
ODVr = -log(F_R/c);
ODVb = -log(F_B/c);

RODV = ODVr./ODVb;
ODcV = corel_fit(RODV);

OD_norma = reshape(ODcV,size(RODV,1),size(RODV,2))/max(ODcV);
%% plot cast 2
figure
plot(OD_norma(457,:))
xlabel('Pozice')
ylabel('OD_{norm}')

figure
plot(OD_norma(:,789))
xlabel('Pozice')
ylabel('OD_{norm}')

figure
imagesc(OD_norma)
colorbar

%% nacteni dat cast B
data0 = load('A1_D0.txt');
data6 = load('A2_D6.txt');
data12 = load('A3_D12.txt');
data18 = load('A4_18.txt');
datax = load('A6_DX.txt');

%% ploty absorpcni spektra
figure
plot(data0(:,1),data0(:,2),'x')
hold on
plot(data12(:,1),data12(:,2),'x')
legend('D=0 Gy','D=6,58 Gy')
xlabel('Vlnová délka [nm]')
ylabel('Absorbance')

figure
hold on
plot(data0(:,1),data0(:,2))
plot(data6(:,1),data6(:,2))
plot(datax(:,1),datax(:,2))
plot(data12(:,1),data12(:,2))
plot(data18(:,1),data18(:,2))
legend('0 Gy','3,29 Gy', '3,34 Gy','6,58 Gy','9,86 Gy')
xlabel('Vlnová délka [nm]')
ylabel('Absorbance')

%% lokální extrémy
[max0,argmax0] = findpeaks(data0(:,2),data0(:,1));
[max6,argmax6] = findpeaks(data6(:,2),data6(:,1));
[max12,argmax12] = findpeaks(data12(:,2),data12(:,1));
[max18,argmax18] = findpeaks(data18(:,2),data18(:,1));
[maxx,argmaxx] = findpeaks(datax(:,2),datax(:,1));

figure
plot(data0(:,1),data0(:,2))
hold on
plot(argmax0,max0,'x')

figure
plot(data6(:,1),data6(:,2))
hold on
plot(argmax6,max6,'x')

figure
plot(data12(:,1),data12(:,2))
hold on
plot(argmax12,max12,'x')

figure
plot(data18(:,1),data18(:,2))
hold on
plot(argmax18,max18,'x')

figure
plot(datax(:,1),datax(:,2))
hold on
plot(argmaxx,maxx,'x')