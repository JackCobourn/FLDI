function AVIwrite(filename,ImageData,Framerate)

v = VideoWriter([filename,'.avi'],'Grayscale AVI');
v.FrameRate=Framerate;
open(v)
writeVideo(v,ImageData)
close(v)

end 
