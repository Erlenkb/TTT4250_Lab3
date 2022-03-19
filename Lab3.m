file1 = "Freq_Response.txt";
file2 = "Beam_Pattern.txt";

Beam_Pattern = importdata(file2, "\t", 1);
Freq_Response = importdata(file1, "\t", 1);
Degrees = [90 75 60 45 30 15 0 -15 -30 -45 -60 -75 -90];
freq = [8000 9000 10000 11000 12000 13000 14000 15000 16000];
Freq_8KHz = Beam_Pattern.data(1,2:27);
Freq_12KHz = Beam_Pattern.data(2,2:27);
Freq_16KHz = Beam_Pattern.data(3,2:27);

deg_0 = Freq_Response.data(1,1:18);
deg_45 = Freq_Response.data(2,1:18);
deg_90 = Freq_Response.data(3,1:18);

scalar1 = 5.44;
scalar2 = 6.36;

Freq_8KHz_mean = zeros(1,13);
Freq_12KHz_mean = zeros(1,13);
Freq_16KHz_mean = zeros(1,13);

deg_0_mean = zeros(1,8);
deg_45_mean = zeros(1,8);
deg_90_mean = zeros(1,8);

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

%% Normalize the two sets by the biggest value for both sets sepertely
Freq_8KHz_mean = log10(Freq_8KHz_mean(Freq_8KHz_mean~=0) ./ scalar2);
Freq_12KHz_mean = log10(Freq_12KHz_mean(Freq_12KHz_mean~=0) ./ scalar2);
Freq_16KHz_mean = log10(Freq_16KHz_mean(Freq_16KHz_mean~=0) ./scalar2);

deg_0_mean = log10(deg_0_mean(deg_0_mean~=0) ./ scalar1);
deg_45_mean = log10(deg_45_mean(deg_45_mean~=0) ./ scalar1);
deg_90_mean = log10(deg_90_mean(deg_90_mean~=0) ./ scalar1);

%% Plot beam patter
figure(1)
subplot(1,2,1)
First = polardb(Degrees, Freq_8KHz_mean, -5, 1,"-g")
hold on
second = polardb(Degrees, Freq_12KHz_mean, -5, 1, "-b")
third = polardb(Degrees, Freq_16KHz_mean, -5, 1, "-r")
hold off
legend("", "", "", "","","","8KHz", "12KHz", "16KHz");

%% Plot frequency response for the three distinct angles
figure(1)
subplot(1,2,2)
plot(freq, deg_0_mean);
hold on
plot(freq, deg_45_mean);
plot(freq, deg_90_mean);
hold off
grid on
xlabel("Frequency [KHz]");
ylabel("Magnitude [dB]");
title("Frequency response");
legend("0 degrees", "45 degrees", "90 degrees");



