function Y = imperlabel_to_labeledim(iY, siz_iY)
% imperlabel_to_labeledim: it does the reverse of labeledim_to_imperlabel
%
% Usage:
%   Y = imperlabel_to_labeledim(iY, siz_iY)
%
% Args:
%   iY: input 2D matrix
%   siz_iY: original dimensions

Y = zeros(siz_iY);

for i = 1:size(iY, 2)
    Y(iY(:, i) ~= 0) = i;
end

end
