function [tform,match_num] = matlabSIFTregist(moving,fixed,options)
arguments
    moving
    fixed
    options.selectStrongest = 500;
    options.ContrastThreshold = 0.0133;
    options.EdgeThreshold = 10;
    options.NumLayersInOctave = 3;
    options.sigma = 1.6;
    options.MatchThreshold = [];
    options.MaxDistance = 1.5; 
    options.Confidence = 99; 
    options.MaxNumTrials = 1000;
    options.kaze = false;
    options.verbose = 0;
end

%%
%MATLABSIFTREGIST Summary of this function goes here
%   Detailed explanation goes here
%     tic;
    if ~options.kaze
        points1 = detectSIFTFeatures(fixed, ContrastThreshold = options.ContrastThreshold, ...
            EdgeThreshold = options.EdgeThreshold, ...
            NumLayersInOctave = options.NumLayersInOctave, sigma = options.sigma);
    else
        points1 = detectKAZEFeatures(fixed,Diffusion='region',Threshold=0.00001,NumOctaves=3,NumScaleLevels=4);
    end

%     toc

    if isempty(options.selectStrongest) || isinf(options.selectStrongest)
        [features1,validPoints1] = extractFeatures(fixed,points1);
    else
        [features1,validPoints1] = extractFeatures(fixed,points1.selectStrongest(options.selectStrongest));
    end
    %% debug
%     figure;imshow(fixed);
%     hold on;
%     plot(points1)
    %%
    % detect SIFT feature for moving
    if ~options.kaze
        points2 = detectSIFTFeatures(moving, ContrastThreshold = options.ContrastThreshold, ...
            EdgeThreshold = options.EdgeThreshold, ...
            NumLayersInOctave = options.NumLayersInOctave, sigma = options.sigma);
    else
        points2 = detectKAZEFeatures(moving,Diffusion='region',Threshold=0.00001,NumOctaves=3,NumScaleLevels=4);
    end

    if isempty(options.selectStrongest) || isinf(options.selectStrongest)
        [features2,validPoints2] = extractFeatures(moving,points2);
    else
        [features2,validPoints2] = extractFeatures(moving,points2.selectStrongest(options.selectStrongest));
    end  
    %% debug
%     figure;imshow(moving);
%     hold on;
%     plot(points2)
    %%
    % match features
    if isempty(options.MatchThreshold)
        indexPairs = matchFeatures(features1,features2,Unique=true);
    else
        indexPairs = matchFeatures(features1,features2,MatchThreshold=options.MatchThreshold,Unique=true);
    end
    matchedPoints1 = validPoints1(indexPairs(:,1));
    matchedPoints2 = validPoints2(indexPairs(:,2));
    %% debug
    if options.verbose
        figure; 
        showMatchedFeatures(fixed,moving,matchedPoints1,matchedPoints2);
    end
    %%
    % tranformation
    [tform,inlierIndex] = estgeotform2d(matchedPoints2,matchedPoints1,'rigid',...
        MaxDistance=options.MaxDistance,Confidence=options.Confidence,MaxNumTrials=options.MaxNumTrials); % estimate transformation
    match_num = sum(inlierIndex);
    disp(['Number of matched features:' num2str(match_num)])

end

