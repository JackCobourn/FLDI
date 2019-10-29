JetDia=(3/16)*25.4;

Max = importdata('maximum.txt',',',5);
Min = importdata('minimum.txt',',',5);

figure(1)
histogram(Max.data(:,2),25)
figure(2)
histogram(Min.data(:,2),25)

AVE_Max=mean(Max.data(:,2));
AVE_Min=mean(Min.data(:,2));
setting_ideal = mean([AVE_Max,AVE_Min]);

set_med = importdata('Setting.txt',',',5);
figure(3)
histogram(set_med.data(:,2),25)
act_setting=mean(set_med.data(:,2))

%%Data List
File_list = ls();
File_list_C = mat2cell(File_list,ones(1,length(File_list)));

Data_list = regexp(File_list_C,'step(?<num>\d*).\w*','match');
Data_num = regexp(File_list_C,'step(?<num>\d*).\w*','tokens');
Data_list = Data_list(~cellfun('isempty',Data_list));
Data_num = Data_num(~cellfun('isempty',Data_num));

for i=1:length(Data_num) %%this whole nighnamre fixes the nested cells from regexp(cells)
    Data_num_M{i,:}=Data_num{i}{1}{1}
    
end
Data=struct('Step',[Data_list{1:end}]','Number',Data_num_M);
for i=1:length(Data_num)
    Data(i).Number=str2num(Data(i).Number);
end

[B,I] = sort([Data.Number]);

Data(:)=Data(I);
%
%% Data Import
for i=1:length(Data_num)
    temp=importdata(sprintf('%s',Data(i).Step),',',5);
    Data(i).time=temp.data(:,1);
    Data(i).Volt=temp.data(:,2);
end
for i=1:length(Data_num)
    Data(i).Loc = Data(i).Number.*2;
    Data(i).RMS = rms(Data(i).Volt-mean(Data(i).Volt));
    Data(i).STD = std(Data(i).Volt-mean(Data(i).Volt));
    Data(i).RMS_T = rms(Data(i).Volt);
end
for i=1:length(Data_num)
    RMS(i)=Data(i).RMS;
    Loc(i)=Data(i).Loc./JetDia;
end

figure(1)
p1 = plot(Loc,RMS,'k*--');
grid on
hold on
p2 = plot(25./JetDia,rms(Max.data(:,2)-mean(Max.data(:,2))),'r*',25./JetDia,rms(Min.data(:,2)-mean(Min.data(:,2))),'g*',25./JetDia,rms(set_med.data(:,2)-mean(set_med.data(:,2))),'c*');
xlabel('Jet position,y/D')
ylabel('Mean-Subtracted RMS Voltage')
legend(sprintf('y/d=%.2f',14.6*10./JetDia),'Wind Off Signal Max','Wind Off Signal Min','Wind Off Signal Med','Location','SouthWest')
xlim([-1,11])

for i=1:length(Data_num)
    figure(i)
    histogram(Data(i).Volt)
    title(sprintf('Histogram of data point %u',i))
end
