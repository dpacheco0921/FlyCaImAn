function [mCo, mCio] = getCfrommC(mergedC, Cpre, Cpost)
% [mCo, mCio]= getCfrommC(mergedC, Cpre, Cpost)
% related to storing merged components
mCo = [];
if ~isempty(mergedC) && iscell(mergedC)
    for j = 1:size(mergedC, 1)
        mCo{j, 1} = Cpre(mergedC{j}, :);
    end
end
if ~isempty(mCo) && ~isempty(Cpost)
        mCio = (size(Cpost, 1) - numel(mCo) + 1):size(Cpost, 1);
else
    mCio = [];
end
end