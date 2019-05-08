% comparison of power response results as determined from Tykla (equation 3) and Vituix CAD
% 5 May 2019, Matthias Brennwald

% Notes:
% - first run the tykla.m program using the test data in the directory "testdata" in order to calculate the power response according to Tykla
% - the power response calculated using Vituix was determined using the same data as those in "testdata"


% load power resonse as calculated in Vituix:
[vf,vPR,vPHASE] = textread ('vituix_results/VituixCAD_Power.txt','%f%f%f','headerlines',1);

% delta = 10*log10(4*pi); % dB offset used to plot the Vituix curve
% delta = 11;


% plot power response curves:
figure(3)

semilogx ( f,PR,'k-' , vf,vPR,'r-' );
xlim([20 20e3])
xlabel ('Frequency (Hz)')
ylabel ('Power response (dB)')
title ('Power Response')
legend ('Tylka','Vituix CAD')
