file1 = "Freq_Response.txt";
file2 = "Beam_Pattern.txt";

%% Global paramaters
Beam_Pattern = importdata(file2, "\t", 1);
Freq_Response = importdata(file1, "\t", 1);

Degrees = [90 75 60 45 30 15 0 -15 -30 -45 -60 -75 -90];
Degrees_rad = Degrees * (pi / 180);

freq = [8 9 10 11 12 13 14 15 16];

%% Fetch the data
Freq_8KHz = Beam_Pattern.data(1,2:27);
Freq_12KHz = Beam_Pattern.data(2,2:27);
Freq_16KHz = Beam_Pattern.data(3,2:27);

deg_0 = Freq_Response.data(1,1:18);
deg_45 = Freq_Response.data(2,1:18);
deg_90 = Freq_Response.data(3,1:18);

%% Scalar values used to normalise the data set
scalar1 = max([deg_0 deg_45 deg_90]);
scalar2 = max([Freq_16KHz Freq_12KHz Freq_8KHz]);
scalar_0  = max(deg_0);
scalar_45 = max(deg_45);
scalar_90 = max(deg_90);
scalar_8 = max(Freq_8KHz);
scalar_12 = max(Freq_12KHz);
scalar_16 = max(Freq_16KHz);


%% Generate two temporary vectors to store the values in
Freq_8KHz_mean = zeros(1,13);
Freq_12KHz_mean = zeros(1,13);
Freq_16KHz_mean = zeros(1,13);

deg_0_mean = zeros(1,8);
deg_45_mean = zeros(1,8);
deg_90_mean = zeros(1,8);

%% Average the data and remove the zeros
for i=[1:2:size(Freq_8KHz,2)]
    Freq_8KHz_mean(i) = mean([Freq_8KHz(i) Freq_8KHz(i+1)]);
    Freq_12KHz_mean(i) = mean([Freq_12KHz(i) Freq_12KHz(i+1)]);
    Freq_16KHz_mean(i) = mean([Freq_16KHz(i) Freq_16KHz(i+1)]);
end
for i=[1:2:size(deg_0,2)]
    deg_0_mean(i) = mean([deg_0(i) deg_0(i+1)]);
    deg_45_mean(i) = mean([deg_45(i) deg_45(i+1)]);
    deg_90_mean(i) = mean([deg_90(i) deg_90(i+1)]);
end

%% Scale values for plot
Freq_8KHz_mean = Freq_8KHz_mean(Freq_8KHz_mean~=0) ./ scalar2;
Freq_12KHz_mean = Freq_12KHz_mean(Freq_12KHz_mean~=0) ./ scalar2;
Freq_16KHz_mean = Freq_16KHz_mean(Freq_16KHz_mean~=0) ./scalar2;

%% Scale values for D


%% Normalize the two sets by the biggest value for both sets sepertely
Freq_8KHz_mean_db = 20*log10(Freq_8KHz_mean);
Freq_12KHz_mean_db = 20*log10(Freq_12KHz_mean);
Freq_16KHz_mean_db = 20*log10(Freq_16KHz_mean);

deg_0_mean = 20*log10(deg_0_mean(deg_0_mean~=0) ./ scalar_0);
deg_45_mean = 20*log10(deg_45_mean(deg_45_mean~=0) ./ scalar_45);
deg_90_mean = 20*log10(deg_90_mean(deg_90_mean~=0) ./ scalar_90);

%% Plot beam patter
figure(1)
subplot(1,2,1)
First = polardb(Degrees, Freq_12KHz_mean_db, -15, 2,"-g")
hold on
second = polardb(Degrees, Freq_8KHz_mean_db, -15, 2, "-b")
third = polardb(Degrees, Freq_16KHz_mean_db, -15, 2, "-r")
hold off
legend("", "", "", "","","","12KHz", "8KHz", "16KHz", "location", "best");
set(gca,'fontsize',12,'fontweight','bold');
%set(gcf,'units','centimeters','position',[2,1,29.7,11.0])

%% Plot frequency response for the three distinct angles
figure(1)
subplot(1,2,2)
semilogx(freq, deg_0_mean);
hold on
semilogx(freq, deg_45_mean);
semilogx(freq, deg_90_mean);

hold off
grid on
xlabel("Frequency [KHz]");
ylabel("Magnitude [dB]");
title("Frequency response");
legend("0 degrees", "45 degrees", "90 degrees","location","best");
set(gca,'fontsize',12,'fontweight','bold');
set(gcf,'units','centimeters','position',[2,1,29.7,11.0])

exportgraphics(figure(1), ['Lab3.png'],'Resolution',450)


%% Calculate the directivity D
%First_Quadrant 
theta_1 = flip(Degrees(1:7),2);
rho_8khz_1 =  flip(Freq_8KHz_mean(1:7) ,2);
rho_12khz_1 = flip(Freq_12KHz_mean(1:7),2);
rho_16khz_1 = flip(Freq_16KHz_mean(1:7),2);

%% D 
D_8khz_1 = 2 / (sum((rho_8khz_1.^2) .* sind(theta_1)));
D_12khz_1 = 2 / (sum((rho_12khz_1.^2) .* sind(theta_1)));
D_16khz_1 = 2 / (sum((rho_16khz_1.^2) .* sind(theta_1)));
%% DI
DI_8kh_1 = 10*log10(D_8khz_1)
DI_12khz_1 = 10*log10(D_12khz_1);
DI_16khz_1 = 10*log10(D_16khz_1);

%Second_Quadrant
theta_2     = Degrees(7:13);
rho_8khz_2  = Freq_8KHz_mean(7:13);
rho_12khz_2 = Freq_12KHz_mean(7:13);
rho_16khz_2 = Freq_16KHz_mean(7:13);

%% D
D_8khz_2 = 2 / (sum((rho_8khz_2.^2) .* sind(theta_1)));
D_12khz_2 = 2 / (sum((rho_12khz_2.^2) .* sind(theta_1)));
D_16khz_2 = 2 / (sum((rho_16khz_2.^2) .* sind(theta_1)));
%% DI
DI_8kh_2 = 10*log10(D_8khz_2)
DI_12khz_2 = 10*log10(D_12khz_2);
DI_16khz_2 = 10*log10(D_16khz_2);


%% Bessel function
a = 0.054;
c = 1464;
%lambda_8khz = 
ka_8khz = (2 * pi * 8000 / c) * (abs(Degrees_rad));
ka_12khz = (2 * pi * 12000 / c) * (abs(Degrees_rad));
ka_16khz = (2 * pi * 16000 / c) * (abs(Degrees_rad));

Bessel_8 = besselj(1, ka_8khz);
Bessel_12 = besselj(1, ka_12khz);
Bessel_16 = besselj(1, ka_16khz);

B_8 = abs(Bessel_8 ./ ka_8khz);
B_12 = abs(Bessel_12 ./ ka_12khz);
B_16 = abs(Bessel_16 ./ ka_16khz);

B_8(7) = 1;
B_12(7) = 1;
B_16(7) = 1;

%% Scale 

Max_all = max([B_8 B_12 B_16]);

B_8_s = 20*log10(B_8 / max(B_8));
B_12_s = 20*log10(B_12 / max(B_12));
B_16_s = 20*log10(B_16 / max(B_16));

B_8 = 20*log10(B_8 / Max_all);
B_12 = 20*log10(B_12 / Max_all);
B_16 = 20*log10(B_16 / Max_all);

%B_8 = B_8 / Max_all;
%B_12 = B_12 / Max_all;
%B_16 = B_16 / Max_all;
N = 120;
step = 30;

figure(5)
subplot(1,3,1)
First = polardb(Degrees, B_8, -N, step,"-b");
legend("", "", "", "","","","8KHz", "location", "best");
set(gcf,'units','centimeters','position',[3,1,29.7,11.0]);
set(gca,'fontsize',12,'fontweight','bold');
figure(5)
subplot(1,3,2)
second = polardb(Degrees, B_12, -N, step, "-g");
legend("", "", "", "","","","12KHz", "location", "best");
set(gcf,'units','centimeters','position',[3,1,29.7,11.0])
set(gca,'fontsize',12,'fontweight','bold');
figure(5)
subplot(1,3,3)
third = polardb(Degrees, B_16, -N, step, "-r");
legend("", "", "", "","","","16KHz", "location", "best");
set(gcf,'units','centimeters','position',[3,1,29.7,11.0]);

%legend("", "", "", "","","","12KHz", "8KHz", "16KHz", "location", "best");
set(gca,'fontsize',12,'fontweight','bold');
exportgraphics(figure(5), ['Math_Directivity.png'],'Resolution',450)



