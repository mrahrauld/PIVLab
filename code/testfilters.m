close all;
RGB = imread('../../000370.jpg');
imshow(RGB);
gray = rgb2gray(RGB);
imshow(gray);
threshold = 90 % custom threshold value
gray(gray < threshold)=0;
imshow(gray);
R = RGB;
R(:,:,[1 3])=0;
%imshow(R);