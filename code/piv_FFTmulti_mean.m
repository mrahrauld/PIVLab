function [xtable ytable utable vtable typevector] = piv_FFTmulti_mean (image1,image2,interrogationarea, step, subpixfinder, mask_inpt, roi_inpt,result_conv,passes,imdeform)
%profile on
%this funtion performs the  PIV analysis.
warning off %#ok<*WNOFF> %MATLAB:log:logOfZero
size(result_conv)
if numel(roi_inpt)>0
    xroi=roi_inpt(1);
    yroi=roi_inpt(2);
    widthroi=roi_inpt(3);
    heightroi=roi_inpt(4);
    image1_roi=double(image1(yroi:yroi+heightroi,xroi:xroi+widthroi));
    image2_roi=double(image2(yroi:yroi+heightroi,xroi:xroi+widthroi));
else
    xroi=0;
    yroi=0;
    image1_roi=double(image1);
    image2_roi=double(image2);
end
gen_image1_roi = image1_roi;
gen_image2_roi = image2_roi;

if numel(mask_inpt)>0
    cellmask=mask_inpt;
    mask=zeros(size(image1_roi));
    for i=1:size(cellmask,1)
        masklayerx=cellmask{i,1};
        masklayery=cellmask{i,2};
        mask = mask + poly2mask(masklayerx-xroi,masklayery-yroi,size(image1_roi,1),size(image1_roi,2)); %kleineres eingangsbild und maske geshiftet
    end
else
    mask=zeros(size(image1_roi));
end
mask(mask>1)=1;
gen_mask = mask;

miniy=1+(ceil(interrogationarea/2));
minix=1+(ceil(interrogationarea/2));
maxiy=step*(floor(size(image1_roi,1)/step))-(interrogationarea-1)+(ceil(interrogationarea/2)); %statt size deltax von ROI nehmen
maxix=step*(floor(size(image1_roi,2)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));

numelementsy=floor((maxiy-miniy)/step+1);
numelementsx=floor((maxix-minix)/step+1);

LAy=miniy;
LAx=minix;
LUy=size(image1_roi,1)-maxiy;
LUx=size(image1_roi,2)-maxix;
shift4centery=round((LUy-LAy)/2);
shift4centerx=round((LUx-LAx)/2);
if shift4centery<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
    shift4centery=0;
end
if shift4centerx<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
    shift4centerx=0;
end
miniy=miniy+shift4centery;
minix=minix+shift4centerx;
maxix=maxix+shift4centerx;
maxiy=maxiy+shift4centery;

image1_roi=padarray(image1_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)], min(min(image1_roi)));
image2_roi=padarray(image2_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)], min(min(image1_roi)));
mask=padarray(mask,[ceil(interrogationarea/2) ceil(interrogationarea/2)],0);

if (rem(interrogationarea,2) == 0) %for the subpixel displacement measurement
    SubPixOffset=1;
else
    SubPixOffset=0.5;
end
xtable=zeros(numelementsy,numelementsx);
ytable=xtable;
utable=xtable;
vtable=xtable;
typevector=ones(numelementsy,numelementsx);

%% MAINLOOP
try %check if used from GUI
    handles=guihandles(getappdata(0,'hgui'));
    GUI_avail=1;
catch %#ok<CTCH>
    GUI_avail=0;
end

% divide images by small pictures
% new index for image1_roi and image2_roi
s0 = (repmat((miniy:step:maxiy)'-1, 1,numelementsx) + repmat(((minix:step:maxix)-1)*size(image1_roi, 1), numelementsy,1))'; 
s0 = permute(s0(:), [2 3 1]);
s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
ss1 = repmat(s1, [1, 1, size(s0,3)])+repmat(s0, [interrogationarea, interrogationarea, 1]);

image1_cut = image1_roi(ss1);
image2_cut = image2_roi(ss1);

%do fft2
% result_conv = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
% minres = permute(repmat(squeeze(min(min(result_conv))), [1, size(result_conv, 1), size(result_conv, 2)]), [2 3 1]);
% deltares = permute(repmat(squeeze(max(max(result_conv))-min(min(result_conv))),[ 1, size(result_conv, 1), size(result_conv, 2)]), [2 3 1]);
% result_conv = ((result_conv-minres)./deltares)*255;


%apply mask
ii = find(mask(ss1(round(interrogationarea/2+1), round(interrogationarea/2+1), :)));
jj = find(mask((miniy:step:maxiy)+round(interrogationarea/2), (minix:step:maxix)+round(interrogationarea/2)));
typevector(jj) = 0;
result_conv(:,:, ii) = 0;

[y, x, z] = ind2sub(size(result_conv), find(result_conv==255));

 % we need only one peak from each couple pictures
[z1, zi] = sort(z);
dz1 = [z1(1); diff(z1)];
i0 = find(dz1~=0);
x1 = x(zi(i0));
y1 = y(zi(i0));
z1 = z(zi(i0));

xtable = repmat((minix:step:maxix)+interrogationarea/2, length(miniy:step:maxiy), 1);
ytable = repmat(((miniy:step:maxiy)+interrogationarea/2)', 1, length(minix:step:maxix));

if subpixfinder==1
    [vector] = SUBPIXGAUSS (result_conv,interrogationarea, x1, y1, z1, SubPixOffset);
elseif subpixfinder==2
    [vector] = SUBPIX2DGAUSS (result_conv,interrogationarea, x1, y1, z1, SubPixOffset);
end
vector = permute(reshape(vector, [size(xtable') 2]), [2 1 3]);

utable = vector(:,:,1);
vtable = vector(:,:,2);


xtable=xtable-ceil(interrogationarea/2);
ytable=ytable-ceil(interrogationarea/2);

xtable=xtable+xroi;
ytable=ytable+yroi;
%profile viewer
%p = profile('info');
%profsave(p,'profile_results')

function [vector] = SUBPIXGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2)));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2);
    if(numel(x)~=0)
        ip = sub2ind(size(result_conv), y, x, z);
        %the following 8 lines are copyright (c) 1998, Uri Shavit, Roi Gurka, Alex Liberzon, Technion � Israel Institute of Technology
        %http://urapiv.wordpress.com
        f0 = log(result_conv(ip));
        f1 = log(result_conv(ip-1));
        f2 = log(result_conv(ip+1));
        peaky = y + (f1-f2)./(2*f1-4*f0+2*f2);
        f0 = log(result_conv(ip));
        f1 = log(result_conv(ip-xmax));
        f2 = log(result_conv(ip+xmax));
        peakx = x + (f1-f2)./(2*f1-4*f0+2*f2);

        SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
        SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
        vector(z, :) = [SubpixelX, SubpixelY];  
    end
    
function [vector] = SUBPIX2DGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2)));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2);
    if(numel(x)~=0)
        c10 = zeros(3,3, length(z));
        c01 = c10;
        c11 = c10;
        c20 = c10;
        c02 = c10;
        ip = sub2ind(size(result_conv), y, x, z);

        for i = -1:1
            for j = -1:1
                %following 15 lines based on
                %H. Nobach � M. Honkanen (2005)
                %Two-dimensional Gaussian regression for sub-pixel displacement
                %estimation in particle image velocimetry or particle position
                %estimation in particle tracking velocimetry
                %Experiments in Fluids (2005) 38: 511�515
                c10(j+2,i+2, :) = i*log(result_conv(ip+xmax*i+j));
                c01(j+2,i+2, :) = j*log(result_conv(ip+xmax*i+j));
                c11(j+2,i+2, :) = i*j*log(result_conv(ip+xmax*i+j));
                c20(j+2,i+2, :) = (3*i^2-2)*log(result_conv(ip+xmax*i+j));
                c02(j+2,i+2, :) = (3*j^2-2)*log(result_conv(ip+xmax*i+j));
                %c00(j+2,i+2)=(5-3*i^2-3*j^2)*log(result_conv_norm(maxY+j, maxX+i));
            end
        end
        c10 = (1/6)*sum(sum(c10));
        c01 = (1/6)*sum(sum(c01));
        c11 = (1/4)*sum(sum(c11));
        c20 = (1/6)*sum(sum(c20));
        c02 = (1/6)*sum(sum(c02));
        %c00=(1/9)*sum(sum(c00));

        deltax = squeeze((c11.*c01-2*c10.*c02)./(4*c20.*c02-c11.^2));
        deltay = squeeze((c11.*c10-2*c01.*c20)./(4*c20.*c02-c11.^2));
        peakx = x+deltax;
        peaky = y+deltay;

        SubpixelX = peakx-(interrogationarea/2)-SubPixOffset;
        SubpixelY = peaky-(interrogationarea/2)-SubPixOffset;

        vector(z, :) = [SubpixelX, SubpixelY];
    end