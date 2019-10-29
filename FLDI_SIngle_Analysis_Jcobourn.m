A2 = load('FLDI-20190529\C21 diode w jet cen00000.dat');
Afs = 1./(mean(diff(A2(:,1))));
LA = length(A2(:,2));
[Ap,Apf] = pwelch(A2(:,2),hann(round(LA/8)),[],[],Afs);

B2 = load('FLDI-20190531\C20531-flowon00000.dat');
Bfs = 1./(mean(diff(B2(:,1))));
LB = length(B2(:,2));
[Bp,Bpf] = pwelch(B2(:,2),hann(round(LB/8)),[],[],Bfs);

C2 = load('FLDI-20190531\C20531-noflow00000.dat');
Cfs = 1./(mean(diff(C2(:,1))));
LC = length(C2(:,2));
[Cp,Cpf] = pwelch(C2(:,2),hann(round(LC/8)),[],[],Cfs);

D2 = A2(A2(:,1)<=-0.047,:);
Dfs = 1./(mean(diff(D2(:,1))));
LD = length(D2(:,2));
[Dp,Dpf] = pwelch(D2(:,2),hann(round(LD/8)),[],[],Dfs);

E2 = A2(A2(:,1)>-0.047,:);
Efs = 1./(mean(diff(E2(:,1))));
LE = length(E2(:,2));
[Ep,Epf] = pwelch(E2(:,2),hann(round(LE/8)),[],[],Efs);

figure(4)
subplot(2,1,1)
plot(D2(:,1),D2(:,2),'r',E2(:,1),E2(:,2),'g')
xlabel('time (s)')
ylabel('Voltage (V)')
title('Raw data')
subplot(2,1,2)
plot(Dpf,Dp,'r',Epf,Ep,'g',Apf,Ap,'b')
ylim([0,1E-7])
xlabel('Hz')
ylabel('PSD')
title('Pwelch()')
sgtitle('M = 1.5 Jet in center of FLDI, data from 05/29')

figure(1)
subplot(2,1,1)
plot(A2(:,1),A2(:,2))
xlabel('time (s)')
ylabel('Voltage (V)')
title('Raw data')
subplot(2,1,2)
plot(Apf,Ap)
ylim([0,2E-7])
xlabel('Hz')
ylabel('PSD')
title('Pwelch()')
sgtitle('M = 1.5 Jet in center of FLDI, data from 05/29')

figure(2)
subplot(2,1,1)
plot(B2(:,1),B2(:,2))
xlabel('time (s)')
ylabel('Voltage (V)')
title('Raw data')
subplot(2,1,2)
plot(Bpf,Bp)
ylim([0,2E-8])
xlabel('Hz')
ylabel('PSD')
title('Pwelch()')
sgtitle('M = 1.5 Jet in center of FLDI, data from 05/31')

figure(3)
subplot(2,1,1)
plot(C2(:,1),C2(:,2))
xlabel('time (s)')
ylabel('Voltage (V)')
title('Raw data')
subplot(2,1,2)
plot(Cpf,Cp)
ylim([0,2E-8])
xlabel('Hz')
ylabel('PSD')
title('Pwelch()')
sgtitle('No flow, data from 05/31')