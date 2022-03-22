file1 = "Freq_Response.txt";
file2 = "Beam_Pattern.txt";

%% Global paramaters
Beam_Pattern = importdata(file2, "\t", 1);
Freq_Response = importdata(file1, "\t", 1);

Degrees = [90 75 60 45 30 15 0 -15 -30 -45 -60 -75 -90];
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

Freq_8KHz_mean = Freq_8KHz_mean(Freq_8KHz_mean~=0) ./ scalar2;
Freq_12KHz_mean = Freq_12KHz_mean(Freq_12KHz_mean~=0) ./ scalar2;
Freq_16KHz_mean = Freq_16KHz_mean(Freq_16KHz_mean~=0) ./scalar2;

%% Normalize the two sets by the biggest value for both sets sepertely
Freq_8KHz_mean_db = 20*log10(Freq_8KHz_mean);
Freq_12KHz_mean_db = 20*log10(Freq_12KHz_mean);
Freq_16KHz_mean_db = 20*log10(Freq_16KHz_mean);

deg_0_mean = log10(deg_0_mean(deg_0_mean~=0) ./ scalar1);
deg_45_mean = log10(deg_45_mean(deg_45_mean~=0) ./ scalar1);
deg_90_mean = log10(deg_90_mean(deg_90_mean~=0) ./ scalar1);

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


