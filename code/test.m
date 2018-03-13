
videoFileReader = vision.VideoFileReader('../../DJI_0005.MOV');
videoPlayer = vision.VideoPlayer('Position',[100,100,680,520]);
fixedFrame = step(videoFileReader);
fixed_ref = imref2d(size(fixedFrame));

imgh = imshow(fixedFrame);
objectRegion=round(getPosition(imrect));
Trackpoints = detectMinEigenFeatures(rgb2gray(fixedFrame),'ROI',objectRegion);
lengths = [0 ;length(Trackpoints)];
fixedPoints = [mean(Trackpoints.Location)];
for i = 1:3
    objectRegion = round(getPosition(imrect));
    pts = detectMinEigenFeatures(rgb2gray(fixedFrame),'ROI',objectRegion);
    Trackpoints = [Trackpoints; pts];
    lengths = [lengths ; length(Trackpoints)];
    fixedPoints = [fixedPoints ; mean(pts.Location)];
end
close(imgh.Parent.Parent);
tic
pointImage = insertMarker(fixedFrame,Trackpoints.Location,'+','Color','white');
tracker = vision.PointTracker('MaxBidirectionalError',1);
initialize(tracker,Trackpoints.Location,fixedFrame);

v = VideoWriter('myFile','Motion JPEG AVI');
open(v)

k = 0;
movingPoints = zeros(4,2);
while ~isDone(videoFileReader)
%while k<100
    k = k+1;
    frame = step(videoFileReader);
    [Trackpoints, validity] = step(tracker,frame);
    for i= 1:4
        movingPoints(i,:) = mean(Trackpoints(lengths(i)+1:lengths(i+1),:));
    end
    t_concord = fitgeotrans(movingPoints,fixedPoints,'projective');
    moving_registered = imwarp(frame,t_concord,'OutputView',fixed_ref);
    imshowpair(moving_registered,fixedFrame,'blend');
    writeVideo(v,moving_registered);
    
    out = insertMarker(frame,movingPoints,'+');
    step(videoPlayer,out);
end
release(videoPlayer);
release(videoFileReader);
close(v)
toc
% moving = imread('000300.jpg');
% fixed = imread('000001.jpg');
% [movingPoints,fixedPoints] = ...
%        cpselect(moving,fixed,...
%                 'Wait',true);
% t_concord = fitgeotrans(movingPoints,fixedPoints,'projective');
% 
% fixed_ref = imref2d(size(fixed)); %relate intrinsic and world coordinates
% moving_registered = imwarp(moving,t_concord,'OutputView',fixed_ref);
% figure, imshowpair(moving_registered,fixed,'blend')