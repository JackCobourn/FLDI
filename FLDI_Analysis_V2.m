%% FLDI Analysis V2
cd FLDI-20190710
ID = 0.176;
Xloc = [(5+5/16),(5+5/16),(5+5/16),7,7,7,8.75,8.75,8.75];
Xloc_d = Xloc./ID;
    
A1 = importdata('20190710_JetPosition1_1.txt',',',5);
A2 = importdata('20190710_JetPosition1_2.txt',',',5);
A3 = importdata('20190710_JetPosition1_3.txt',',',5);
B1 = importdata('20190710_JetPosition2_1.txt',',',5);
B2 = importdata('20190710_JetPosition2_2.txt',',',5);
B3 = importdata('20190710_JetPosition2_3.txt',',',5);
C1 = importdata('20190710_JetPosition3_1.txt',',',5);
C2 = importdata('20190710_JetPosition3_2.txt',',',5);
C3 = importdata('20190710_JetPosition3_3.txt',',',5);
N1 = importdata('20190710_Nothing.txt',',',5);
F1 = importdata("20190710 Flame.txt",',',5);

cd ..

A1MAX = max(A1.data(:,2));
A1MIN =  min(A1.data(:,2));
A1AVE = mean(A1.data(:,2));

A2MAX = max(A2.data(:,2));
A2MIN =  min(A2.data(:,2));
A2AVE = mean(A2.data(:,2));

A3MAX = max(A3.data(:,2));
A3MIN =  min(A3.data(:,2));
A3AVE = mean(A3.data(:,2));

B1MAX = max(B1.data(:,2));
B1MIN =  min(B1.data(:,2));
B1AVE = mean(B1.data(:,2));

B2MAX = max(B2.data(:,2));
B2MIN =  min(B2.data(:,2));
B2AVE = mean(B2.data(:,2));

B3MAX = max(B3.data(:,2));
B3MIN =  min(B3.data(:,2));
B3AVE = mean(B3.data(:,2));

C1MAX = max(C1.data(:,2));
C1MIN =  min(C1.data(:,2));
C1AVE = mean(C1.data(:,2));

C2MAX = max(C2.data(:,2));
C2MIN =  min(C2.data(:,2));
C2AVE = mean(C2.data(:,2));

C3MAX = max(C3.data(:,2));
C3MIN =  min(C3.data(:,2));
C3AVE = mean(C3.data(:,2));

N1MAX = max(N1.data(:,2));
N1MIN =  min(N1.data(:,2));
N1AVE = mean(N1.data(:,2));

F1MAX = max(F1.data(:,2));
F1MIN =  min(F1.data(:,2));
F1AVE = mean(F1.data(:,2));

MaxVolt = max([A1MAX,A2MAX,A3MAX,B1MAX,B2MAX,B3MAX,C1MAX,C2MAX,C3MAX,N1MAX,F1MAX]);
MinVolt = min([A1MIN,A2MIN,A3MIN,B1MIN,B2MIN,B3MIN,C1MIN,C2MIN,C3MIN,N1MIN,F1MIN]);
AveVolt = mean([A1AVE,A2AVE,A3AVE,B1AVE,B2AVE,B3AVE,C1AVE,C2AVE,C3AVE]);

cA = (MaxVolt-MinVolt)./2;
cD = (MaxVolt+MinVolt)./2;

Pint = acos(((AveVolt-cD)./cA));

A1P = acos((A1.data(:,2)-cD)./cA)-Pint;
A2P = acos((A2.data(:,2)-cD)./cA)-Pint;
A3P = acos((A3.data(:,2)-cD)./cA)-Pint;
B1P = acos((B1.data(:,2)-cD)./cA)-Pint;
B2P = acos((B2.data(:,2)-cD)./cA)-Pint;
B3P = acos((B3.data(:,2)-cD)./cA)-Pint;
C1P = acos((C1.data(:,2)-cD)./cA)-Pint;
C2P = acos((C2.data(:,2)-cD)./cA)-Pint;
C3P = acos((C3.data(:,2)-cD)./cA)-Pint;
N1P = acos((N1.data(:,2)-cD)./cA)-Pint;
F1P = acos((F1.data(:,2)-cD)./cA)-Pint;

A1_var = var(A1.data(:,2));
A2_var = var(A2.data(:,2));
A3_var = var(A3.data(:,2));
B1_var = var(B1.data(:,2));
B2_var = var(B2.data(:,2));
B3_var = var(B3.data(:,2));
C1_var = var(C1.data(:,2));
C2_var = var(C2.data(:,2));
C3_var = var(C3.data(:,2));

A1P_AVE = mean(A1P);
A2P_AVE = mean(A2P);
A3P_AVE = mean(A3P);
B1P_AVE = mean(B1P);
B2P_AVE = mean(B2P);
B3P_AVE = mean(B3P);
C1P_AVE = mean(C1P);
C2P_AVE = mean(C2P);
C3P_AVE = mean(C3P);
N1P_AVE = mean(N1P);
F1P_AVE = mean(F1P);

A1P_VAR = var(A1P);
A2P_VAR = var(A2P);
A3P_VAR = var(A3P);
B1P_VAR = var(B1P);
B2P_VAR = var(B2P);
B3P_VAR = var(B3P);
C1P_VAR = var(C1P);
C2P_VAR = var(C2P);
C3P_VAR = var(C3P);
N1P_VAR = var(N1P);
F1P_VAR = var(F1P);

A1P_rms = rms(A1P-mean(A1P));
A2P_rms = rms(A2P-mean(A2P));
A3P_rms = rms(A3P-mean(A3P));
B1P_rms = rms(B1P-mean(B1P));
B2P_rms = rms(B2P-mean(B2P));
B3P_rms = rms(B3P-mean(B3P));
C1P_rms = rms(C1P-mean(C1P));
C2P_rms = rms(C2P-mean(C2P));
C3P_rms = rms(C3P-mean(C3P));
N1P_rms = rms(N1P-mean(N1P));
F1P_rms = rms(F1P-mean(F1P));

figure(1)
plot(smoothdata([A1.data(:,2),A2.data(:,2),A3.data(:,2),B1.data(:,2),B2.data(:,2),B3.data(:,2),C1.data(:,2),C2.data(:,2),C3.data(:,2),N1.data(:,2),F1.data(:,2)]));
legend('Location 1 data A','Location 1 data B','Location 1 data C','Location 2 data A','Location 2 data B','Location 2 data C','Location 3 data A','Location 3 data B','Location 3 data C', 'No Flow', 'Flame')

figure(2)
%plot(Xloc(1:3),[A1P_AVE,A2P_AVE,A3P_AVE],'b*')
plot(Xloc_d([1,4,7]),[A1P_AVE,B1P_AVE,C1P_AVE],'b*',Xloc_d([1,4,7]),[A2P_AVE,B2P_AVE,C2P_AVE],'r*',Xloc_d([1,4,7]),[A3P_AVE,B3P_AVE,C3P_AVE],'g*')
legend('First sample','Second Sample','Third Sample')  
xlabel('x/d')
ylabel('\Delta P, Ave')

figure(8)
plot(Xloc_d([1,4,7]),[A1P_VAR,B1P_VAR,C1P_VAR],'b*',Xloc_d([1,4,7]),[A2P_VAR,B2P_VAR,C2P_VAR],'r*',Xloc_d([1,4,7]),[A3P_VAR,B3P_VAR,C3P_VAR],'g*')
legend('First sample','Second Sample','Third Sample')  
xlabel('x/d')
ylabel('\Delta P, var')


figure(9)
plot(Xloc_d([1,4,7]),[A1P_rms,B1P_rms,C1P_rms],'b*',Xloc_d([1,4,7]),[A2P_rms,B2P_rms,C2P_rms],'r*',Xloc_d([1,4,7]),[A3P_rms,B3P_rms,C3P_rms],'g*')
legend('First sample','Second Sample','Third Sample')  
xlabel('x/d')
ylabel('\Delta P, mean subtracted rms')


figure(3)
plot([A1P_AVE,A2P_AVE,A3P_AVE,B1P_AVE,B2P_AVE,B3P_AVE,C1P_AVE,C2P_AVE,C3P_AVE],'b*')

figure(4)
plot(Xloc,[A1P_AVE,A2P_AVE,A3P_AVE,B1P_AVE,B2P_AVE,B3P_AVE,C1P_AVE,C2P_AVE,C3P_AVE],'b*')

%figure(5)
%plot(Xloc,[A1P_AVE,B1P_AVE,C1P_AVE,A2P_AVE,B2P_AVE,C2P_AVE,A3P_AVE,B3P_AVE,C3P_AVE],'b*')

[Pw,f] = pwelch([A1P,A2P,A3P,B1P,B2P,B3P,C1P,C2P,C3P],[],[],[],1./abs((A1.data(1,1)-A1.data(2,1))));
figure(6)
hold on
plot(f(7:end),Pw(7:end,[1,2,3]),'b')
plot(f(7:end),Pw(7:end,[4:6]),'r')
plot(f(7:end),Pw(7:end,[7:9]),'g')
legend('Location 1 data A','Location 1 data B','Location 1 data C','Location 2 data A','Location 2 data B','Location 2 data C','Location 3 data A','Location 3 data B','Location 3 data C')
xlabel('Hz')
ylabel('')