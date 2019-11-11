function [A_or, C_or, srt, srt_val] = order_ROIs_int(A, C, type)
%[A_or, C_or, srt, srt_val] = order_ROIs_int(A, C, type)
% ordering of the found components based on their maximum temporal
% activation and their size (through their l_inf norm)
% you can also pre-specify the ordering sequence
if ~exist('type', 'var'); type = 0; end
%% get l_inf norm
nA = full(sqrt(sum(A.^2)));
nr = length(nA);
A = A/spdiags(nA(:),0,nr,nr);
mA = sum(A.^4).^(1/4);
%% normalize C
if type == 0
    mC = max(C,[],2);
else
    C = spdiags(nA(:),0,nr,nr)*C;
    mC = max(C,[],2);
end
%% get and apply rank
if ~exist('srt', 'var')||isempty(srt)
    [srt_val, srt] = sort(mC.*mA', 'descend');
end
A_or = A(:,srt);
C_or = C(srt,:);
srt_val = full(srt_val);
end