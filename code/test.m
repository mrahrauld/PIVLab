I = imread('fisheye.JPG');
% [imagePoints,boardSize] = detectCheckerboardPoints(images.Files);
% squareSize = 29; % millimeters
% boardSize
% worldPoints = generateCheckerboardPoints(boardSize,squareSize);
% I = readimage(images,1);
% imageSize = [size(I,1) size(I,2)];
% params = estimateFisheyeParameters(imagePoints,worldPoints,imageSize);
[J,newOrigin] = undistortImage(I,cameraParams);
figure
imshow(I);
figure
imshow(J);
title('Full Output View')

tgzraehifoùAFZGHRBFJZADJF=GZQJFkjazryGSVQ