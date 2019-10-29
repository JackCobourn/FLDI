A2 = load('C2tea-lens.dat'); %Data from a candle near the field lens
A3 = load('C3tea-lens.dat');

B2 = load('C2tea-cen.dat'); %Data from a candle near the center of the FLDI system
B3 = load('C3tea-cen.dat');

C2 = load('C2jet-lens.dat'); %Data from a jet near the field lens
C3 = load('C3jet-lens.dat');

D2 = load('C2jet-cen.dat'); %Data from a jet near the center of the FLDI system
D3 = load('C3jet-cen.dat');

DelPhiA = asin(2.*((A2(:,2)-A3(:,2))./(A2(:,2)+A3(:,2)))-1);
DelPhiB = asin(2.*((B2(:,2)-B3(:,2))./(B2(:,2)+B3(:,2)))-1);
DelPhiC = asin(2.*((C2(:,2)-C3(:,2))./(C2(:,2)+C3(:,2)))-1);
DelPhiD = asin(2.*((D2(:,2)-D3(:,2))./(D2(:,2)+D3(:,2)))-1);

rho = 1.225; %kg/m3 at stp;
lambda = 632.8E-9; %m
K = 2.257E-4; %Gladstone dale const for air
dX = 2E-3; %estimated speration distance

FFTA = fft(DelPhiA);
NA = length(FFTA);

FFTB = fft(DelPhiB);
NB = length(FFTB);

FFTC = fft(DelPhiC);
NC = length(FFTC);

FFTD = fft(DelPhiD);
ND = length(FFTD);


TIA = (1/rho)*sqrt((1/NA^2)*sum(abs(lambda.*(1/(2*pi*K*dX)).*FFTA).^2))*100;
TIB = (1/rho)*sqrt((1/NB^2)*sum(abs(lambda.*(1/(2*pi*K*dX)).*FFTB).^2))*100;
TIC = (1/rho)*sqrt((1/NC^2)*sum(abs(lambda.*(1/(2*pi*K*dX)).*FFTC).^2))*100;
TID = (1/rho)*sqrt((1/ND^2)*sum(abs(lambda.*(1/(2*pi*K*dX)).*FFTD).^2))*100;