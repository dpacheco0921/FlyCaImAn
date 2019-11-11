function [mC, mCi] = update_mC(mC, mCi, Cinit, C, imCi)
% [mC, mCi] = update_mC(mC, mCi, Cinit, C, imCi)
% This function updates the cells that stores merged C (mC), and their indexes to the C they generate
% inputs: 
% mC = cell with all the C that orignated a new component mCi (index of the current C)
% mCi = index of the new component in Cinit.
% Cinit = previous to merging.
% C = Cinit after merging.
% imCi = index of new merged components
if ~isempty(imCi)
    for j = 1:numel(imCi)
        % see if its the result of a merged component
        igate = numel(intersect(imCi{j}, mCi));
        if igate > 0
            % collect previous merges
            [~, ~, ib] = intersect(imCi{j}, mCi);
            tmC{j, 1} = cell2mat(mC(ib));
            % add new components
            [~, ia] = setdiff(imCi{j}, mCi);
            tmC{j, 1} = [tmC{j, 1}; Cinit(imCi{j}(ia), :)];
            % delete old ones and their indexes
            mCi(ib) = []; mC(ib) = [];
        else
            tmC{j, 1} = Cinit(imCi{j}, :);
        end
    end
    mC = [mC; tmC];
    mCi = [mCi, (size(C, 1) - numel(imCi) + 1):size(C, 1)];
end
end