function [mCi, mCi2del] = update_mCi(mCi, i2del)
% This function updates the mCi depending on how many components below this
% value are eliminated
mCi2del = [];
if ~isempty(i2del)
    [~, mCi2del, ~] = intersect(mCi, i2del);
    for i = 1:numel(mCi)
        mCi(i) = mCi(i) - sum(i2del < mCi(i));
    end
end
end