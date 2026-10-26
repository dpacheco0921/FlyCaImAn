function p = parse_opto_fictrac_txtfile(filename)
% parse_opto_fictrac_txtfile: parse control file with stimulus information returns a structure p with
%   fieldnames given by first word of headers in control file and contents
%   given by subsequent rows
%
% USAGE
%     p = parse_opto_fictrac_txtfile(filename)
%
% PARAMETERS
%     filename - name of the tab-delimited control (.txt) file, with the
%                       path if it is not in the current folder
%
% RETURNS
%     p - structure with one field per column of the control file
%         (named after the first word of the header, e.g. stimFileName,
%         intensity, rate, silencePre, silencePost), keeping only the
%         rows of the opto stimulus (second stimulus, after ';').
%         p.intensity is replaced by its second value (opto intensity).

% parse inputs
controlFilePath = filename;

% scan file
fid = fopen(controlFilePath);

if fid < 0
    error('parse_opto_fictrac_txtfile:nofile', ...
        'Could not open control file %s', controlFilePath)
end

% parse the first line as header
Chead = textscan(fid, '%s%s%s%s%s%s%s%*[^\n]', ...
    1, 'CommentStyle','#', 'Delimiter','\t');

% parse all following lines as content
Cbody = textscan(fid, '%s%f%f%f%f%f%s%*[^\n]', ...
    'CommentStyle','#', 'Delimiter','\t');
fclose(fid);

for c = 1:length(Chead)    
    
    % use first "word" of the header name only to allow for units etc in the header names
    fieldName = strsplit(Chead{c}{1},' ');                                                
    
    % populate param struct, use header as field name and body as field contents
    p.(fieldName{1}) = Cbody{c};
    
end

% opto stimuli is always the second stimuli (after ';')
% (cell2mat would drop rows without ';' and misalign the rows, so test
%   every row separately)
stimtype = cellfun(@(c) ~isempty(c) && c(1) == 1, ...
    strfind(p.stimFileName, ';'));

% remove non-opto stimulus
for c = 1:length(Chead)    
    
    fieldName = strsplit(Chead{c}{1}, ' ');                                                
    p.(fieldName{1}) = p.(fieldName{1})(stimtype);
    
end

if isempty(p.intensity)
    error('parse_opto_fictrac_txtfile:noopto', ...
        'No opto stimulus rows (stimFileName starting with '';'') in %s', ...
        controlFilePath)
end

for i = 1:numel(p.intensity)
    intensity_num(i, :) = eval(p.intensity{i});
end

p.intensity = intensity_num(:, 2);

end
