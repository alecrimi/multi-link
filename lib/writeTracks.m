function  writeTracks( fileName, dictIdLst, volId, trkInfo, dictVec  )
%writeTracks Summary of this function goes here
%   Detailed explanation goes here

    %addpath('/mnt/data/matlabLibs/along-tract-stats');

    nAreas = length(dictIdLst);
    trkIds = [];
    for idArea=1:nAreas
        % translate area to connection
        connId = dictVec(dictIdLst(idArea));
        % find track ids
        trkIdsForConn = find (trkInfo.connToTrkIdLst{volId} == connId);
        trkIds = [trkIds; trkIdsForConn]; %#ok<AGROW>
        % find track ids Rev (if they exist)
        if( isfield(trkInfo, 'connToTrkIdLstRev') )
            trkIdsForConn = find (trkInfo.connToTrkIdLstRev{volId} == connId);
            trkIds = [trkIds; trkIdsForConn]; %#ok<AGROW>
        end
    end
    
    if( ~isempty(trkIds) )
        % make trkIds unique
        trkIds = unique(trkIds);
    
        %-- add only tracks identified by trkIds
        tkvHeaderNew = trkInfo.tkvHeaderLst{volId};
        tkvHeaderNew.n_count = length(trkIds);
        tkvTracksNew(tkvHeaderNew.n_count).nPoints = 0; % init structure
        voxSize = tkvHeaderNew.voxel_size;
        for tIdx=1:length(trkIds)
            currTrk = trkInfo.tkvTracksLst{volId}(trkIds(tIdx));
            % rescale points
            currTrk.matrix = currTrk.matrix .* repmat(voxSize,currTrk.nPoints,1 );
            % add whole track to structure
            tkvTracksNew(tIdx).nPoints = currTrk.nPoints;
            tkvTracksNew(tIdx).matrix = currTrk.matrix;
        end
        % writing volume
        display( ['Writing '  fileName] );
        trk_write_new(tkvHeaderNew, tkvTracksNew, fileName );
    else
        display( [ fileName ' has no tracts in the areas selected'] );
    end
end

