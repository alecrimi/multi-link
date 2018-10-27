function  [featMat, classVec, dictVec, volInfo, spaceInfo, trkInfo] = ...
    buildFeatures_wm( minLength, cubeSize )
%BUILDFEATURES build tracts feature vector
%featMat -> featId x classes
%classVec -> classes
%dictMat -> maps the featId to the index of the connections
%volInfo -> volInfo.volNames names of the volumes, volInfo.regMat registration matrices
%spaceInfo-> spaceInfo.globMaxPnt,spaceInfo.globMinPnt, spaceInfo.partNum, spaceInfo.cubeSize , spaceInfo.minLength
%trkInfo-> trkInfo.tkvHeader, trkInfo.tkvTracks, trkInfo.connToTrkIdLst

% This script requires some files of the Along-tract-stats of John Colby https://github.com/johncolby/along-tract-stats

    %global partNum globMaxPnt globMinPnt cubeSize %#ok<REDEF>
    global partNum
    
    %-- Param
    if( nargin < 2)
        minLength = 15; % minimum length of tract
        cubeSize = 10; % size of a cube 
    end
    
    tractDir = 'dti_res_orig/';

    
    % single affine matrix (ridig)
    matDir = 'mat_wm_orig_sm_rig/'; 
    tractsFileSuffix = 'regRigid_dipy_affOnly.trk';
%     % double affine matrix (non-ridig)
%     matDir = 'mat_wm_orig/'; 
%     tractsFileSuffix = 'regNonRigid_dipy.trk';
%     % single affine matrix (non-ridig)
%     matDir = 'mat_wm_orig_sm/'; 
%     tractsFileSuffix = 'regNonRigid_dipy_sm.trk';
    
%     % % Domain for cubeSize -> 5
%     globMaxPnt = [134, 115, 158]';
%     globMinPnt = [29, -5, 3]';

    % % Domain for cubeSize -> 10 (adapted to template)
    globMaxPnt = [100, 85, 105]';
    globMinPnt = [20, 35, 5]';
    
    partNum = (globMaxPnt - globMinPnt) / cubeSize; % partitions number
    
    if( any(fix(partNum) - partNum) )
        error('parts must be integers');
    end
    %--
    
    
    %---Subject names and classes
    % All volumes 
    volNames = {'ag120213a'  'ag120307d'  'ag120329c'  'ag120403c'  'ag120410a'  'ag120416b' ...
                 'ag120214c'  'ag120309b'  'ag120330b'  'ag120404d'  'ag120411b'  'ag120417c' ...
                 'ag120215c'  'ag120313c'  'ag120402c'  'ag120405b'  'ag120412c'};
     
    classVec = [0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1]; %0 -> c57, 1 -> BTB 
  
    %store volNames
    volInfo.volNames = volNames;
    %---
    
    
    %-- Store space info
    spaceInfo.globMaxPnt = globMaxPnt;
    spaceInfo.globMinPnt = globMinPnt;
    spaceInfo.partNum = partNum;
    spaceInfo.cubeSize = cubeSize;
    spaceInfo.minLength = minLength;
    %--
    
    % Histogram of all possible connections
    totConn = (partNum(1) * partNum(2) * partNum(3))^2;
    connHistGlob = sparse(totConn,1);
    % Holders 
    connHistLst = cell(length(volNames), 1);
    for volIdx=1:length(volNames)
        %for connections in each volume
        connHistLst{volIdx} = sparse(totConn,1);
        % for track references
        trkInfo.connToTrkIdLst{volIdx} = [];
        % for tracks
        trkInfo.tkvTracksLst{volIdx} = [];
        % for track Headers
        trkInfo.tkvHeaderLst{volIdx} = [];
    end
    
    for volIdx=1:length(volNames)
        display(['vol ' num2str(volIdx) ' of ' num2str(length(volNames))]);
        subjName = volNames{volIdx};
        fileMat = [matDir subjName '_l' num2str(cubeSize) '_' num2str(minLength) '.mat'];
        
        %-- Get stEndMat
        stEndMat = [];
        stEndMatReg = [];
        if( ~exist(fileMat, 'file') )
            tractsFile = [tractDir subjName '/' tractsFileSuffix];
            [stEndMat, tkvHeader, tkvTracks] = generateStEndMat ( tractsFile, minLength );  
            trkInfo.tkvHeaderLst{volIdx} = tkvHeader;
            trkInfo.tkvTracksLst{volIdx} = tkvTracks;
            % save
            save(fileMat, 'stEndMat', 'tkvHeader', 'tkvTracks');
        else
            % load
            loaded = load(fileMat);
            stEndMat = loaded.stEndMat;
            trkInfo.tkvHeaderLst{volIdx} = loaded.tkvHeader;
            trkInfo.tkvTracksLst{volIdx} = loaded.tkvTracks;
        end
        %--

        %-- store without registration
        stEndMatReg(1:3, :) = stEndMat(1:3,:);
        stEndMatReg(4:6, :) = stEndMat(4:6,:);
        %--
        
%         %-- Register to common space with rigid transformation (WM)        
%         reg1FileName = findFile( subjName, reg1Dir );
%         reg1Mat = importdata(reg1FileName);
%         reg2FileName = findFile( subjName, reg2Dir );
%         reg2Mat = importdata(reg2FileName);
%         
%         regMat = reg2Mat * reg1Mat;
%         resSt = regMat * [stEndMat(1:3,:); ones(1,size(stEndMat,2))];
%         resEnd = regMat * [stEndMat(4:6,:); ones(1,size(stEndMat,2))];
%         stEndMatReg(1:3, :) = resSt(1:3,:);
%         stEndMatReg(4:6, :) = resEnd(1:3,:);
%         %store regMat
%         volInfo.regMats{volIdx} = regMat;
%         %--
        
%         %-- Register to common space (PRNI)
%         regFileName = 'data_registered.txt';
%         regMat = importdata([tractDir subjName '_iso/' regFileName]);
% 
%         resSt = regMat * [stEndMat(1:3,:); ones(1,size(stEndMat,2))];
%         resEnd = regMat * [stEndMat(4:6,:); ones(1,size(stEndMat,2))];
%         stEndMatReg(1:3, :) = resSt(1:3,:);
%         stEndMatReg(4:6, :) = resEnd(1:3,:);
%         %store regMat
%         volInfo.regMats{volIdx} = regMat;
%         %--
        
        % set domain
        maxPnt = globMaxPnt;
        minPnt = globMinPnt;
        
        %-- Fill in histogram all connections
        trkInfo.connToTrkIdLst{volIdx} = zeros( size(stEndMatReg,2), 1 ); % map connection id to trkId
        for i=1:size(stEndMatReg,2)
            % find sampled positions
            stPos = ceil(((stEndMatReg(1:3,i) - minPnt) ./ (maxPnt-minPnt)) .* partNum);
            endPos = ceil(((stEndMatReg(4:6,i) - minPnt) ./ (maxPnt-minPnt)) .* partNum);

            % avoid 0s
            stPos(stPos==0) = 1;
            endPos(endPos==0) = 1;
            
            % remove fibres out of range
            if( any(stEndMatReg(:,i) > [globMaxPnt;globMaxPnt]) || any(stEndMatReg(:,i) < [globMinPnt;globMinPnt])  )
                continue;
            end
            
            % Add to matrices
            connId = getIdx( stPos(1),stPos(2),stPos(3),endPos(1),endPos(2),endPos(3) );
            connHistGlob( connId ) = connHistGlob( connId ) + 1;
            connHistLst{volIdx}( connId ) = connHistLst{volIdx}( connId ) + 1;
                        
            % map connectionId to trackid
            trkInfo.connToTrkIdLst{volIdx}(i) = connId;
        end
        %--
    end
    
    % dictionary vector
    dictVec = find(connHistGlob>0);
    
    % fill in featMat
    featMat = zeros(length(dictVec), length(volNames));
    for volIdx=1:length(volNames)
        featMat(:,volIdx) = connHistLst{volIdx}(dictVec);
    end
    
    % transform connectionId with dictionary
%     connToTrkId = dictVec(connToTrkId);
    % store track information
%     trkInfo.tkvHeader = tkvHeader;
%     trkInfo.tkvTracks = tkvTracks;
%     trkInfo.connToTrkIdLst = connToTrkIdLst;
    
end

function linIdx = getIdx( fX, fY, fZ, tX, tY, tZ )
    global partNum 
    
    linIdx = sub2ind([partNum; partNum]', fX, fY, fZ, tX, tY, tZ);
end

function [fX, fY, fZ, tX, tY, tZ] = revIdx( linIdx )
    global partNum 
    
    [fX, fY, fZ, tX, tY, tZ] = ind2sub([partNum; partNum]', linIdx);
end

function [stPnt, endPnt] = revIdxDeReg( linIdx, dictVec, regMat, spaceInfo )

    % init
    globMaxPnt = spaceInfo.globMaxPnt;
    globMinPnt = spaceInfo.globMinPnt;
    partNum = spaceInfo.partNum;
    cubeSize = spaceInfo.cubeSize;
    
    % get indexes
    [fX_, fY_, fZ_, tX_, tY_, tZ_] = revIdx( dictVec(linIdx) );
    
    % convert to real coordinates
    totSize = (globMaxPnt-globMinPnt);
    stIdx = [fX_, fY_, fZ_]';
    endIdx = [tX_, tY_, tZ_]';
%     stPntReg = globMinPnt + ((stIdx./partNum) .* totSize) - cubeSize/2;
%     endPntReg = globMinPnt + ((endIdx./partNum) .* totSize) - cubeSize/2;
    stPntReg = globMinPnt + ((stIdx./partNum) .* totSize);
    endPntReg = globMinPnt + ((endIdx./partNum) .* totSize);


    % Unregister
    regMatI = inv(regMat);
    stPnt = regMatI * [stPntReg;1];
    endPnt = regMatI * [endPntReg;1];
    
    % Convert to coordinates (valid in trackvis)
    stPnt =  diag([1.29,1.29,1.29,1]) \ stPnt; % inv(diag([1.29,1.29,1.29,1])) * stPnt;
    endPnt = diag([1.29,1.29,1.29,1]) \ endPnt; %inv(diag([1.29,1.29,1.29,1])) * endPnt;
    
end



function [stEndMat, tkvHeader, tkvTracks] = generateStEndMat ( tractsFile, minLength )
        % % Fsl (to be fixed)

        % read tracts
%         [header, tracks] = trk_read_raw(tractsFile);
        [tkvHeaderOrig, tkvTracksOrig] = trk_read(tractsFile);
        voxSize = tkvHeaderOrig.voxel_size;

        %- Filter by lengths
        tractLengths = trk_length(tkvTracksOrig);
        tractLengthsLbl = tractLengths>minLength;
        goodTractsIDs = find( tractLengthsLbl );
        %-

        %- find tracts domain, create startEnd matrix and create new tkv structures
        tkvHeader = tkvHeaderOrig;
        tkvHeader.n_count = length(goodTractsIDs);
        tkvTracks(tkvHeader.n_count).nPoints = 0; % init struture
        stEndMat = zeros(6, tkvHeader.n_count);
        parfor gtIdx=1:tkvHeader.n_count
            idxTkvTracksOrig = goodTractsIDs(gtIdx);
            currTrk = tkvTracksOrig( idxTkvTracksOrig );
            % rescale points
            currTrk.matrix = currTrk.matrix ./ repmat(voxSize,currTrk.nPoints,1 );
            % add start and end point
%             stEndMat(1:3,gtIdx) = currTrk.matrix(1,:)';
%             stEndMat(4:6,gtIdx) = currTrk.matrix(end,:)';
            stEndMat(:,gtIdx) = [currTrk.matrix(1,:)';currTrk.matrix(end,:)'];
            % add whole track to structure
            tkvTracks(gtIdx).nPoints = currTrk.nPoints;
            tkvTracks(gtIdx).matrix = currTrk.matrix;
        end
        %-
end


function fName = findFile( subjName, regDir )
    fName = ls( [regDir '/*' subjName '*'] );
    fName = fName(1:end-1);
end

