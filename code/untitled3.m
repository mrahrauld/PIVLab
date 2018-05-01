file = ''
videoFileReader = vision.VideoFileReader(file + '.MOV');
videoPlayer = vision.VideoPlayer('Position',[100,100,680,520]);
fixedFrame = step(videoFileReader);
imwrite(fixedFrame,file+'.png');