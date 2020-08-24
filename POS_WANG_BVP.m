function [BVP, check] = POS_WANG_BVP(VideoFile, FS, roiDetAlg)
% POS_WANG The Plane Orthogonal to Skin-Tone (POS) Method from: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491. DOI: 10.1109/TBME.2016.2609282
%
%   Inputs:
%       VideoFile               = Video file path.
%       FS                      = Video framerate (fps).
%
%   Outputs:
%       BVP                     = Processed Blood Volume Pulse (BVP).
%
%   Requires - Signal Processing Toolbox
%
% Adjusted from: Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

%% Parameters
SkinSegmentTF = true;

WinSec=1.6;%(based on refrence's 32 frame window with a 20 fps camera)

%% Load Video:
VidObj = VideoReader(VideoFile);
% VidObj.CurrentTime = StartTime;

FramesToRead = floor(VidObj.Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Read Video and Spatially Average:
T = zeros(FramesToRead,1);%initialize time vector
RGB = zeros(FramesToRead,3);%initialize color signal
FN = 0;

% Create a cascade detector object.
roiDetector = vision.CascadeObjectDetector(roiDetAlg);

bbox_last = [0 0 VidObj.Height VidObj.Width];

save_every_n_frame = 100;
check = struct();

while hasFrame(VidObj) %&& (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in
    %reference as a CSK detector from Henriques et al., 2012
    % Read a video frame and run the detector.
    if FN == 1 || mod(FN, 20) == 0
        bbox = step(roiDetector, VidFrame);
        if size(bbox, 1) > 0
            bbox_last = bbox(end, :);
        else
            bbox = bbox_last;
        end
    end
    
    VidROI = VidFrame(bbox(end,2):bbox(end,2)+bbox(end,4), bbox(end,1):bbox(end,1)+bbox(end,3), :);
    
    
    
    if(SkinSegmentTF)%skin segmentation - originally specified in reference as an OC-SVM from Wang et al. 2015
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
    else
        RGB(FN,:) = sum(sum(VidROI,2)) ./ (size(VidROI,1)*size(VidROI,2));
    end
    
    if mod(FN, save_every_n_frame) == 0
        check(end+1).frameROI = VidROI;
        if SkinSegmentTF
            check(end).skinROI = ROISkin;
        end
        check(end).bbox = bbox;
        check(end).frameIndex = FN;
    end
end
%% POS:
% Transform from: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017, May). Color-distortion filtering for remote photoplethysmography. In Automatic Face & Gesture Recognition (FG 2017), 2017 12th IEEE International Conference on (pp. 71-78). IEEE.
useFGTransform=false;
if useFGTransform
    RGBBase = mean(RGB);
    RGBNorm = bsxfun(@times,RGB,1./RGBBase)-1;
    FF = fft(RGBNorm);
    F = (0:size(RGBNorm,1)-1)*FS/size(RGBNorm,1);
    H = FF*[-1/sqrt(6);2/sqrt(6);-1/sqrt(6)];
    W = (H.*conj(H))./sum(FF.*conj(FF),2);
    FMask = (F >= LPF)&(F <= HPF);
    % FMask(length(FMask)/2+1:end)=FMask(length(FMask)/2:-1:1);
    FMask = FMask + fliplr(FMask);
    W=W.*FMask';%rectangular filter in frequency domain - not specified in original paper
    FF = FF.*repmat(W,[1,3]);
    RGBNorm=real(ifft(FF));
    RGBNorm = bsxfun(@times,RGBNorm+1,RGBBase);
    
    RGB=RGBNorm;
end

%lines and comments correspond to pseudo code algorithm on reference page 7       
N = size(RGB,1);%line 0 - RGB is of length N frames
H = zeros(1,N);%line 1 - initialize to zeros of length of video sequence
l = ceil(WinSec*FS);%line 1 - window length equivalent to reference: 32 samples of 20 fps camera (default 1.6s)
C = zeros(length(l),3);
for n = 1:N-1%line 2 - loop from first to last frame in video sequence
    %line 3 - spatial averaging was performed when video was read
    m = n - l + 1;%line 4 condition
    if(m > 0)%line 4
        Cn = ( RGB(m:n,:) ./ mean(RGB(m:n,:)) )';%line 5 - temporal normalization
        S = [0, 1, -1; -2, 1, 1] * Cn;%line 6 - projection
        h = S(1,:) + ((std(S(1,:)) / std(S(2,:))) * S(2,:));%line 7 - tuning
        H(m:n) = H(m:n) + (h - mean(h));%line 8 - overlap-adding
    end%line 9 - end if
end%line 10 - end for

BVP=H;
end%end function
