function [BVP, check] = CHROM_DEHAAN_BVP(VideoFile, FS, roiDetAlg, useFGT)
% useFGT is used for consistency with POS_WANG_BVP function
% CHROM_DEHAAN The Chrominance Method from: De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886. DOI: 10.1109/TBME.2013.2266196
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
% Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

%% Parameters
SkinSegmentTF = true;

LPF = 0.7; %low cutoff frequency (Hz) - specified as 40 bpm (~0.667 Hz) in reference
HPF = 4; %high cutoff frequency (Hz) - specified as 240 bpm (~4.0 Hz) in reference

WinSec=1.6; %(was a 32 frame window with 20 fps camera)

%% Load Video:
VidObj = VideoReader(VideoFile);
% VidObj.CurrentTime = StartTime;

FramesToRead=floor(VidObj.Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Read Video and Spatially Average:
T = zeros(FramesToRead,1);%initialize time vector
RGB = zeros(FramesToRead,3);%initialize color signal
FN = 0;

% Create a cascade detector object.
if ~strcmp(roiDetAlg, 'NoDet')
    roiDetector = vision.CascadeObjectDetector(roiDetAlg);
end

bbox_last = [1 1 VidObj.Width-1 VidObj.Height-1];
bbox = bbox_last;

save_every_n_frame = 200;
check = struct();
cnt = 0;

while hasFrame(VidObj) %&& (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    
    %position for optional face detection/tracking - originally specified in
    %reference as a CSK detector from Henriques et al., 2012
    % Read a video frame and run the detector.
    if FN == 1 || mod(FN, 20) == 0
        if ~strcmp(roiDetAlg, 'NoDet')
            bbox = step(roiDetector, VidFrame);
            if size(bbox, 1) > 0
                bbox_last = bbox(end, :);
            else
                bbox = bbox_last;
            end
        end
    end
    
    VidROI = VidFrame(bbox(end,2):bbox(end,2)+bbox(end,4), bbox(end,1):bbox(end,1)+bbox(end,3), :);
    
    
    if(SkinSegmentTF)%skin segmentation - not specified in reference
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
        cnt = cnt + 1;
        check(cnt).frameROI = VidROI;
        if SkinSegmentTF
            check(cnt).skinROI = ROISkin;
        end
        check(cnt).bbox = bbox;
        check(cnt).frameIndex = FN;
    end
end
%% CHROM:
NyquistF = 1/2*FS;
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter - originally specified as an a FIR band-pass filter with cutoff frequencies 40-240 BPM

%Window parameters - overlap, add with 50% overlap
WinL = ceil(WinSec*FS);
if(mod(WinL,2))%force even window size for overlap, add of hanning windowed signals
    WinL=WinL+1;
end
NWin = floor((FN-WinL/2)/(WinL/2));
S = zeros(NWin,1);
WinS = 1;%Window Start Index
WinM = WinS+WinL/2;%Window Middle Index
WinE = WinS+WinL-1;%Window End Index

for i = 1:NWin
    TWin = T(WinS:WinE,:);
    
    RGBBase = mean(RGB(WinS:WinE,:));
    RGBNorm = bsxfun(@times,RGB(WinS:WinE,:),1./RGBBase)-1;
    
    % CHROM
    Xs = squeeze(3*RGBNorm(:,1)-2*RGBNorm(:,2));%3Rn-2Gn
    Ys = squeeze(1.5*RGBNorm(:,1)+RGBNorm(:,2)-1.5*RGBNorm(:,3));%1.5Rn+Gn-1.5Bn
    
    Xf = filtfilt(B,A,double(Xs));
    Yf = filtfilt(B,A,double(Ys));
    
    Alpha = std(Xf)./std(Yf);
    
    SWin = Xf - Alpha.*Yf;
    
    SWin = hann(WinL).*SWin;
    %overlap, add Hanning windowed signals
    if(i==1)
        S = SWin;
        TX = TWin;
    else
        S(WinS:WinM-1) = S(WinS:WinM-1)+SWin(1:WinL/2);%1st half overlap
        S(WinM:WinE) = SWin(WinL/2+1:end);%2nd half
        TX(WinM:WinE) = TWin(WinL/2+1:end);
    end
    
    WinS = WinM;
    WinM = WinS+WinL/2;
    WinE = WinS+WinL-1;
end

BVP=S;
end%end function
