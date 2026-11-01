function [result, SI, startTime] = readScanImageTiffLinLab(fileName)
% readScanImageTiffLinLab:
% Uses ScanImageTiffReader to read in a ScanImage file with the name
% fileName
% We use this because scanimage.util.opentif is unbearably slow
% 
% Usage:
% [result, SI] = readScanImageTiffLinLab(fileName);
% - result is a a 5D matrix in the order: x, y, z, channels, time
% If any of these are singletons, eg only one channel, they are kept as 
% singletons (the result is not squeezed)
% - SI holds the metadata as original created by ScanImage (possibly with
% some extra fields if the file was modified by our Fiji macros)
%
% Usage notes:
% If the metadata isn't there (eg because the .tif was created in Fiji),
% the function looks for an accompanying metadata text file
% If fileName is example.tif, the metadata file should be called
% example-metadata.txt
% If you lack the metadata file from a file created in Fiji, use the macro 
% 'saveMetadataAsText' in Fiji

% Written by Andrew Lin, andrew.lin@sheffield.ac.uk

tic

% un-comment the line below if your ScanImageTiffReader is from 2021 or later
import ScanImageTiffReader.ScanImageTiffReader;

reader = ScanImageTiffReader(fileName);

metadataText = reader.metadata;
desc = reader.descriptions();

if isempty(metadataText)
    % look for the accompanying text file
    suffixIndex = strfind(fileName,'.');
    suffixIndex = suffixIndex(end);
    metadataFileName = strcat(fileName(1:(suffixIndex-1)),'-metadata.txt');
    [SI,startTime] = parseFijiScanImageMetadata(metadataFileName);
else
    [SI,startTime] = parseScanImageMetadata(metadataText, desc);
end

movie = reader.data();
dimensions = size(movie);

% linesPerFrame = SI.hRoiManager.linesPerFrame;
% pixelsPerLine = SI.hRoiManager.pixelsPerLine;
numChannels = length(SI.hChannels.channelSave);

% in raw ScanImage data, the discarded flyback slices are retained in the
% movie file
% thus, if you want 10 slices and have 2 flyback slices, the parameters
% will be:
% SI.hFastZ.discardFlybackFrames = true
% SI.hStackManager.numSlices = 10
% SI.hFastZ.numFramesPerVolume = 12
% SI.hFastZ.numDiscardFlybackFrames = 2
% however, the Lin lab Z-registration Fiji plugin deletes the blank slices
% so in the above example, there would only be 10 slices, and no blank
% slices
% The code below deletes the blank slices but only if they are there.

% first check if ScanImage blanked the flyback frames, then check if the
% flybak frames have been deleted. If not, delete them

% in SI.TIFF_FORMAT_VERSION 3, numFramesPerVolume is in
% SI.hFastZ.numFramesPerVolume and it includes the flyback frames
% for backward compatbility, in in SI.TIFF_FORMAT_VERSION 4 we will use 
% SI.hStackManager.numFramesPerVolumeWithFlyback, which also includes the
% fly back frames

if SI.TIFF_FORMAT_VERSION==4
    numFramesPerVolume = SI.hStackManager.numFramesPerVolumeWithFlyback;
    numSlices = SI.hStackManager.actualNumSlices;
else
    numFramesPerVolume = SI.hFastZ.numFramesPerVolume;
    numSlices = SI.hStackManager.numSlices;
end

if (SI.hFastZ.discardFlybackFrames && ~SI.LinLab.flybackFramesDeleted)
    if ((numFramesPerVolume - SI.hFastZ.numDiscardFlybackFrames) ~= SI.hStackManager.numSlices)
        numFramesPerVolume
        SI.hFastZ.numDiscardFlybackFrames
        numSlices
        error('Something is wrong with the ScanImage z parameters\n');
    end
    % reshape into dimensions: x y c z t
    result = reshape(reader.data(), dimensions(1), dimensions(2), ...
       numChannels,  numFramesPerVolume, ...
        ceil(dimensions(3)/(numFramesPerVolume*numChannels)));
    % now delete the last SI.hFastZ.numDiscardFlybackFrames z-slices
    result = result(:,:,:,1:(end-SI.hFastZ.numDiscardFlybackFrames),:);
else
    if mod(dimensions(3), numSlices*numChannels)
        error('numFrames does not divide evenly by numSlices and numChannels');
    end
    % reshape into dimensions: x y c z t
    result = reshape(reader.data(), dimensions(1), dimensions(2), ...
        numChannels, numSlices, ceil(dimensions(3)/(numSlices*numChannels)));
end


% reorder into x y z c t
% the reason reshape used x y c z t instead of x y z c t is that
% channels are acquired simultaneously whereas z is acquired
% sequentially using the piezo. Thus channels are closer together in
% the interleaved file than z slices.
result = permute(result, [1 2 4 3 5]);

if SI.TIFF_FORMAT_VERSION==4
    % Ariel's measurements 251028
    % at 8x zoom, 512 x 512 pixels, ~0.225 microns per pixel
    SI.LinLab.rowCal = 0.225 * SI.hRoiManager.scanAngleMultiplierFast / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.pixelsPerLine/512);
    SI.LinLab.columnCal = 0.225 * SI.hRoiManager.scanAngleMultiplierSlow / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.linesPerFrame/256);
else
    % x and y calibration based on measurements of 170929
    % at 8x zoom, 256 x 256 pixels, 0.32 microns per pixel
    % this is the first dimension which is x or fast scan (but
    % Matlab treats it as rows)
    SI.LinLab.rowCal = 0.32 * SI.hRoiManager.scanAngleMultiplierFast / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.pixelsPerLine/256);
    % this is the second dimension which is y or slow scan (but
    % Matlab treats it as columns)
    SI.LinLab.columnCal = 0.32 * SI.hRoiManager.scanAngleMultiplierSlow / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.linesPerFrame/256);
end

disp('Total time to read in file was: ');
toc
end