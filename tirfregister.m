
function ISOR=tirfregister(nd2file,imregalgotf)

%Register nd2 movie (or tif timeseries stack) from TIRF experiment
% containing 3 channels (DNA, ch1, ch2)by Fourier transform drift
% correction or translation affine transform
% INPUT                 -nd2file: multichannel nd2 timelapse dataset
%                       -imregalgotf: registering options by drift
%                       correction (0), or affine transform (1)
%
% OUTPUT                ISOR: x-by-y-by-Timeframes-by-1-by-numchannels
%                       matrix containing the registered image data
%                       The image data is further saved as .ome.tiff


close all;
delete *.ome.tiff % registered files are save as .ome.tiff

IS = bfopen(nd2file);
B = cat(3,IS{1,1}{1:3:end,1}); % DNA channel (#1)
R = cat(3,IS{1,1}{2:3:end,1}); % channel#2
G = cat(3,IS{1,1}{3:3:end,1}); % channel#3

[x, y, t] = size(B);

Breg = zeros(size(B),'like',B);
Breg(:,:,1) = B(:,:,1);
Rreg = zeros(size(R),'like',R);
Rreg(:,:,1) = R(:,:,1);
Greg = zeros(size(G),'like',G);
Greg(:,:,1) = G(:,:,1);

if imregalgotf == 1

    tformall = imreg2(B); %getting the geometric transformations

    for i = 2:t
        Breg(:,:,i) = imwarp(B(:,:,i),tformall{i},'OutputView',imref2d([x y]));
        Rreg(:,:,i) = imwarp(R(:,:,i),tformall{i},'OutputView',imref2d([x y]));
        Greg(:,:,i) = imwarp(G(:,:,i),tformall{i},'OutputView',imref2d([x y]));
    end

elseif imregalgotf == 0

    shiftcoord = imdrift(B,0); %getting the xy shift for drift correction

    for i = 2:t
        if any(shiftcoord{i})
            S = shiftcoord{i};
            SY = S(1); SX = S(2);
            Br = circshift(B(:,:,i),[SY SX]);
            Rr = circshift(R(:,:,i),[SY SX]);
            Gr = circshift(G(:,:,i),[SY SX]);
            if SY > 0
                Br(1:SY,:) = 0;
                Rr(1:SY,:) = 0;
                Gr(1:SY,:) = 0;
            elseif SY < 0
                Br(end+SY+1:end,:) = 0;
                Rr(end+SY+1:end,:) = 0;
                Gr(end+SY+1:end,:) = 0;
            end
            if SX > 0
                Br(:,1:SX) = 0;
                Rr(:,1:SX) = 0;
                Gr(:,1:SX) = 0;
            elseif SX < 0
                Br(:,end+SX+1:end) = 0;
                Rr(:,end+SX+1:end) = 0;
                Gr(:,end+SX+1:end) = 0;
            end
            Breg(:,:,i) = Br;
            Rreg(:,:,i) = Rr;
            Greg(:,:,i) = Gr;
        else
            Breg(:,:,i) = B(:,:,i);
            Rreg(:,:,i) = R(:,:,i);
            Greg(:,:,i) = G(:,:,i);
        end

    end

end
% Cropping the resulting image
BFcrop = min(Breg,[],3);
xminp = sum(BFcrop,1);
xtf = xminp == 0;
yminp = sum(BFcrop,2);
ytf = yminp == 0;

Breg(:,xtf,:) = [];
Breg(ytf,:,:) = [];
Rreg(:,xtf,:) = [];
Rreg(ytf,:,:) = [];
Greg(:,xtf,:) = [];
Greg(ytf,:,:) = [];

[sx, sy, st] = size(Breg);

ISOR = zeros(sx, sy, st, 1, 3,'like',Breg);
ISOR(:,:,:,1,1) = Breg;
ISOR(:,:,:,1,2) = Rreg;
ISOR(:,:,:,1,3) = Greg;

%saving the file in .ome.tiff
fname = nd2file;
sname = [fname(1:end-4) '.ome.tiff'];
disp('BEWARE: saving the file may take several minutes!')
bfsave(ISOR, sname, 'dimensionOrder', 'XYTZC');

end

%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shiftcoord = imdrift(imstack,modetf)

% Demo of drift correction
% 2018/11/15  by Renfong


% INPUT:    -imstack 3D x-by-y-by-t grayscale image stack to correct
%           -modetf drift correction mode 0: standard (Ref is t0) 1:
%           recursive (Ref is t-1)
%
% OUTPUT:   -shiftcoord numframe-by-1 cell arrray vector containing 2-by-1
% vector of pixel drift correction in Y and X direction


[sy, sx, ni] = size(imstack);

shiftcoord = cell(ni,1);
progressText(0,'Drift Correction');

switch modetf
    case 0
        imfixed = imstack(:,:,1);
        
        for i= 2:ni
            immoving = imstack(:,:,i);
            % Using fft cross correlations to detect the moving distance[1]
            fftim1 = fft2(imfixed);
            fftim2 = fft2(immoving);
            cc = fftshift(ifft2(fftim1 .* conj(fftim2)));
            [shiftY,shiftX] = find(cc == max(cc(:)));
            shiftY = shiftY - fix(sy / 2) - 1;
            shiftX = shiftX - fix(sx / 2) - 1;
            
            shiftcoord{i} = [shiftY shiftX];
            progressText(i / ni);
        end
        
    case 1
        imfixed = imstack(:,:,1);
        
        for i = 2:ni
            immoving = imstack(:,:,i);
            % Using fft cross correlations to detect the moving distance[1]
            fftim1 = fft2(imfixed);
            fftim2 = fft2(immoving);
            cc = fftshift(ifft2(fftim1 .* conj(fftim2)));
            [shiftY,shiftX] = find(cc == max(cc(:)));
            shiftY = shiftY - fix(sy / 2) - 1;
            shiftX = shiftX - fix(sx / 2) - 1;
            
            shiftcoord{i} = [shiftY shiftX];
            imfixed = circshift(immoving,[shiftY,shiftX]);
            progressText(i / ni);
        end
        
end

end


function tformall=imreg2(imstack)

% INPUT:    -imstack  3D x-by-y-by-t grayscale image stack to correct
%
% OUTPUT:   -tformall  numframe-by-1 cell arrray vector containing
% translation transformation matrix for each time frame


ni = size(imstack,3);

tformall = cell(ni,1);
imfixed = imstack(:,:,1);
tforminit = affine2d;
tformall{1,2} = tforminit;

[o, m] = imregconfig('monomodal');
progressText(0,'Image Registration');
for j = 2:ni
    immoving = imstack(:,:,j);
    tform = imregtform(immoving,imfixed,'translation',o,m);
    tformall{j,1} = tform;

    progressText(j / ni);
end

end