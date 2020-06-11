function AVIwrite(filename,ImageData)

v = VideoWriter([filename,'.AVI'],'Grayscale AVI');
open(v)
writeVideo(v,ImageData)
close(v)

end 
