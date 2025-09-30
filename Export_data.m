
ex=Raw_1;
% for i=1:size(ex,3)
%     ex(:,:,i) = rescale(ex(:,:,i),0,1);
% end


outputFolder = 'exported_images';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


for idx = 1:size(ex,3)
    % ??????
    currentImageData = ex(:,:,idx);
    currentImage = (currentImageData - min(currentImageData(:))) / (max(currentImageData(:)) - min(currentImageData(:)));
     currentImage = ind2rgb(im2uint8(currentImage), jet(256));
%     thresholdValue = graythresh(currentImage);
%     binaryImage = imbinarize(currentImage,thresholdValue);
    % ??????????????image_001.png, image_002.png??
    fileName = sprintf('image_%03d.png', idx);
    
    % ?????????
    outputPath = fullfile(outputFolder, fileName);
    
    % ???????
%     imwrite(binaryImage, outputPath);
    imwrite(currentImage, outputPath);
end

disp('over');














