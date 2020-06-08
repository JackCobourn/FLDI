function RGB = UTKcolors(NameorNumber)
ListofNames = {'Valley','Torch','Globe','Limestone','River',...
    'Leconte','Regalia','Sunsphere','Rock','Legacy','Summitt',...
    'Buckskin','Energy','Switchgrass','Fountain','Eureka!','Orange','White','Smokey'};

ListofRGB = [0 116 111; 230 89 51; 0 108 147; 240 237 227; 81 124 150;...
    141 32 72; 117 74 126; 254 213 53; 167 169 172; 87 149 132; 185 225 226;...
    112 85 80; 238 62 128; 171 193 120; 33 151 169; 235 234 100; 255 130 0; 255 255 255; 88 89 91]./255;
if ischar(NameorNumber)
    if sum(nonzeros(strcmp(ListofNames,NameorNumber)))>0
        Num = find(strcmp(ListofNames,NameorNumber));
    else
        error('typo in name')
    end
elseif ceil(NameorNumber)<=length(ListofRGB)
    Num = NameorNumber; 
else
    error('unknown imput')
end

RGB = ListofRGB(Num,:);
    
end