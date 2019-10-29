cd 'FLDI-20190521'
A2 = load('C2tea-lens.dat'); %Data from a candle near the field lens
A3 = load('C3tea-lens.dat');
Afs = 1./(mean([mean(diff(A2(:,1))),mean(diff(A3(:,1)))]));

B2 = load('C2tea-cen.dat'); %Data from a candle near the center of the FLDI system
B3 = load('C3tea-cen.dat');
Bfs = 1./(mean([mean(diff(B2(:,1))),mean(diff(B3(:,1)))]));

C2 = load('C2jet-lens.dat'); %Data from a jet near the field lens
C3 = load('C3jet-lens.dat');
Cfs = 1./(mean([mean(diff(C2(:,1))),mean(diff(C3(:,1)))]));

D2 = load('C2jet-cen.dat'); %Data from a jet near the center of the FLDI system
D3 = load('C3jet-cen.dat');
Dfs = 1./(mean([mean(diff(D2(:,1))),mean(diff(D3(:,1)))]));


[A,Af] = mscohere(A2(:,2),A3(:,2),[],[],[],Afs,'onesided');
[B,Bf] = mscohere(B2(:,2),B3(:,2),[],[],[],Bfs,'onesided');
[C,Cf] = mscohere(C2(:,2),C3(:,2),[],[],[],Cfs,'onesided');
[D,Df] = mscohere(D2(:,2),D3(:,2),[],[],[],Dfs,'onesided');

% [Ap,Apf] = cpsd(A2(:,2),A3(:,2),[],[],[],Afs,'onesided');
% [Bp,Bpf] = cpsd(B2(:,2),B3(:,2),[],[],[],Bfs,'onesided');
% [Cp,Cpf] = cpsd(C2(:,2),C3(:,2),[],[],[],Cfs,'onesided');
% [Dp,Dpf] = cpsd(D2(:,2),D3(:,2),[],[],[],Dfs,'onesided');

[Ap,Apf] = pwelch(A2(:,2),[],[],[],Afs,'onesided');
[Bp,Bpf] = pwelch(B2(:,2),[],[],[],Bfs,'onesided');
[Cp,Cpf] = pwelch(C2(:,2),[],[],[],Cfs,'onesided');
[Dp,Dpf] = pwelch(D2(:,2),[],[],[],Dfs,'onesided');

figure(1)
plot(Af,A)


figure(2)
plot(Bf,B)
%ylim([0,2E-7])

figure(3)
plot(Cf,C)
%ylim([0,2E-7])

figure(4)
plot(Df,D)
%ylim([0,2E-7])

figure(5)
plot(Apf,Ap)
ylim([0,2E-7])

figure(6)
plot(Bpf,Bp)
ylim([0,2E-7])

figure(7)
plot(Cpf,Cp)
ylim([0,2E-7])

figure(8)
plot(Cpf,Cp)
ylim([0,2E-7])

figure(9)
plot(A2(:,1),A2(:,2),A3(:,1),A3(:,2))

figure(10)
plot(B2(:,1),B2(:,2),B3(:,1),B3(:,2))

figure(11)
plot(C2(:,1),C2(:,2),C3(:,1),C3(:,2))

figure(12)
plot(D2(:,1),D2(:,2),D3(:,1),D3(:,2))