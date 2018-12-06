function run_example_MA()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Luca Giancardo MIT-IIT gianca@mit.edu
% and by Alessandro Crimi IIT alessandro.crimi@iit.it

% All source code unless otherwise specified in the source code are copyrighted to Luca Giancardo

% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    
    %--- Imports
    addpath('lib');
    load; 
    %--

    display('run trainTestSLDA, writeTracksSel  or compareSLDAtTest:');
    display(' ');
    display('    [sldaFeat, sldaWeight, sldaWeightStd] = trainTestSLDA( featMat, classVec ); % To identify the features using SLDA');
    display(' ');
    display('    writeTracksSel( ''BTBR-control'', sldaFeat, trkInfo, dictVec ); % To save the indified features for Trackvis');
    keyboard;
end



function featSelTtest = selWithTtest( featMat, classVec, pVal )
    if( nargin < 3 )
        pVal = 0.0001;
    end
    [~,p] = ttest2(featMat(:,classVec==0)', featMat(:,classVec==1)');
    featSelTtest = find(p < pVal);
end

function [featSelOut, weightFeatSelOutAvg, weightFeatSelOutStd] = trainTestSLDA( featMat, classVec, stopVal, minRoundNum , trkInfo, dictVec)
    if( nargin < 4 )
        stopVal = -150;
        minRoundNum = size(featMat,2)/2;
    end
 
    % set parameters
    lambda = 1e-6; % l2-norm
    % stopVal; % l1-norm. negative: number of vars in LARS-EN
    maxiter = 250; % max iter in SDCA alg.§
    
    % Invert featMat to samples X feat
    featMat = featMat';
    
    %---- Perform leave-2-out
    featSelLst = [];
    featSelWeightLst = [];
    
    classIdxPos = find(classVec==0);
    classIdxNeg = find(classVec==1);
    %if( length(classIdxPos)~=length(classIdxNeg))
    %    error('num of classes not even');
    %end
    for i=1:length(classIdxPos)
        % convert class ids to cells
        classIdxPosCell = mat2cell( classIdxPos, 1,ones(1,length(classIdxPos)) );
        classIdxNegCell = mat2cell( classIdxNeg, 1,ones(1,length(classIdxNeg)) );
        % get testing ids
        testClassIdx = [classIdxPos(i), classIdxNeg(i) ];
        % get training ids
        classIdxPosCell(i) = [];
        classIdxNegCell(i) = [];
        trainClassIdx = [cell2mat( classIdxPosCell ), cell2mat( classIdxNegCell )];
        
        % split featMat
        featMatTrain = featMat(trainClassIdx, :);
        featMatTest = featMat(testClassIdx, :);
        
        % Create class labels for SLDA
        classY = zeros(size(featMatTrain,1), 2);
        classY(classVec(trainClassIdx)==0,1) = 1;
        classY(classVec(trainClassIdx)==1,2) = 1;

        % normalise
        [X mu d] = normalize( featMatTrain );

        % Train SLDA
        [sl theta ridgeCost] = slda( X, classY, lambda, stopVal, 1, maxiter, 1e-6,0 );
  
        % test SLDA
        featMatTestNorm = normalizeWithOldData( featMatTest, mu, d );
        
        % Project data onto the sparse directions (dim=1)
        DC = featMatTestNorm*sl;
           
        featSel = find(sl~=0);
        featSelWeight = sl(featSel);
        
        featSelLst = [featSelLst; featSel(:)];
        featSelWeightLst = [featSelWeight; featSelWeightLst(:)]; % store weight corresponding to featSelLst
    end
    
    % select features
    uFeat = unique(featSelLst(:));
    hist( featSelLst(:), uFeat );
    xlabel('Features ID');
    ylabel('Occurrence of features');
    [occ, bins] = hist( featSelLst(:), uFeat );
    [occSorted,idxSorted] = sort(occ, 'descend');
    featSelOut = bins( idxSorted( find(occSorted>=minRoundNum) ) );
    
    % calculate average feature weight
    weightFeatSelOutAvg = zeros(size(featSelOut));
    weightFeatSelOutStd = zeros(size(featSelOut));
    for i=1:length(featSelOut)
        currWeights = featSelWeightLst(featSelOut(i) == featSelLst);
        weightFeatSelOutAvg(i) = mean(currWeights);
        weightFeatSelOutStd(i) = std(currWeights);
    end
    
    % order features selected by weight
    [~,idxWeightSorted] = sort(abs(weightFeatSelOutAvg), 'descend');
    weightFeatSelOutAvg = weightFeatSelOutAvg(idxWeightSorted);
    weightFeatSelOutStd = weightFeatSelOutStd(idxWeightSorted);
    featSelOut = featSelOut(idxWeightSorted);

    disp('Output features computed, the histogram shows the occurences of the features');
 
end


function Xout = normalizeWithOldData( X, mu, d )
% normalise with the trained vectors
% X -> samples x features
% mu -> mean
% d -> std
    n = size(X,1);
    Xout = X - repmat(mu, n, 1);
    Xout = Xout./(ones(n,1)*d);
end

% write tracks for trackvis
function writeTracksSel( lbl, featSel, trkInfo, dictVec )
    dirOutVol = {   'output_data/1' ... %1
                    'output_data/2' ... %2
                    'output_data/3' ... %3
                    'output_data/4' ... %4
                    'output_data/5' ... %5
                    'output_data/6' ... %6 
                    'output_data/7' ... %7
                    'output_data/8' ... %8
                    'output_data/9' ... %9
                    'output_data/10' ... %10
                    'output_data/11' ... %11
                    'output_data/12' ... %12 
                    'output_data/13' ... %13
                    'output_data/14' ... %14
                    'output_data/15' ... %15
                    'output_data/16' ... %16
        };  
    for volIdx=1:length(dirOutVol)
        writeTracks( [ dirOutVol{volIdx} '/' lbl '_vol' num2str(volIdx) '.trk'], featSel, volIdx, trkInfo, dictVec  );
    end

end
 

function compareSLDAtTest( featMat, classVec )
    % SLDA param
    stopVal = -150;
    minRoundNum = 8;
    %
    alphaVec = [1:5:200];
    intersecting = zeros(1,length(alphaVec));
    intersecting2 = zeros(1,length(alphaVec));
 
    for kk = 1:length(alphaVec)
        kk
        stopVal=alphaVec(kk);
        % run SLDA
        [sldaFeat, sldaWeight, sldaWeightStd] = trainTestSLDA( featMat, classVec, -stopVal, minRoundNum );
        % and ttest
        [~,p] = ttest2(featMat(:,classVec==0)', featMat(:,classVec==1)');
        [~,ttestFeat] = sort(p, 'ascend'); % sort features based on their significance

        intersecting(kk)  = abs( length(intersect(sldaFeat, ttestFeat(1:length(sldaFeat)))) -length(sldaFeat)) ; % Detected by SLDA but rejected by T-test
        intersecting2(kk) = abs( length(intersect(sldaFeat, ttestFeat(1:length(sldaFeat)))) -(sum(p<0.001))) ;  % Detected by T-test but rejected by SLDA
    end
 
    % plot intersection
    figure;
    hold on
         plot(alphaVec, intersecting , 'r' );
         plot(alphaVec, intersecting2, 'b' );
    hold off
    ylabel('Features detected by SLDA but rejected by T-test');
    xlabel('L-1 regularization parameter');
 
end
