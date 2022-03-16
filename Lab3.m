file1 = "DAta_Lab3.csv";
file2 = "Beam_Pattern.txt";

Beam_Pattern = importdata(file2, "\t", 1);
Degrees = [90 75 60 45 30 15 0 -15 -30 -45 -60 -75 -90];
Freq_8KHz = Beam_Pattern.data(1,2:27);
Freq_12KHz = Beam_Pattern.data(2,2:27);
Freq_16KHz = Beam_Pattern.data(3,2:27);
scalar = 6.36;

Freq_8KHz_mean = zeros(1,13);
Freq_12KHz_mean = zeros(1,13);
Freq_16KHz_mean = zeros(1,13);
%cnt = 0;
for i=[1:2:size(degree_0,2)]

    Freq_8KHz_mean(i) = mean([Freq_8KHz(i) Freq_8KHz(i+1)]);
    Freq_12KHz_mean(i) = mean([Freq_12KHz(i) Freq_12KHz(i+1)]);
    Freq_16KHz_mean(i) = mean([Freq_16KHz(i) Freq_16KHz(i+1)]);
    %cnt = int8(cnt+1);
end
Freq_8KHz_mean = log10(Freq_8KHz_mean(Freq_8KHz_mean~=0) ./ scalar);
Freq_12KHz_mean = log10(Freq_12KHz_mean(Freq_12KHz_mean~=0) ./ scalar);
Freq_16KHz_mean = log10(Freq_16KHz_mean(Freq_16KHz_mean~=0) ./scalar);

First = polardb(Degrees, Freq_8KHz_mean, -5, 1,"-g")
hold on
second = polardb(Degrees, Freq_12KHz_mean, -5, 1, "-b")
third = polardb(Degrees, Freq_16KHz_mean, -5, 1, "-r")
hold off
legend("", "", "", "","","","8KHz", "12KHz", "16KHz");

