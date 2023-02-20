function p = parse_opto_fictrac_txtfile(filename)
% parse_opto_fictrac_txtfile: parse control file with stimulus information returns a structure p with
%   fieldnames given by first word of headers in control file and contents
%   given by subsequent rows
%
% USAGE
%     p = parse_opto_fictrac_txtfile(filename)
%
% PARAMETERS
%     filename - file name, needs to contain the full path to the
%                       control file if no p is passed, otherwise will get
%                       path from p

% parse inputs
controlFilePath = filename;

% scan file
% control file needs to be in ctrlDir
fid = fopen(controlFilePath);

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
stimtype = cell2mat(strfind(p.stimFileName, ';'));
stimtype = stimtype == 1;

% remove non-opto stimulus
for c = 1:length(Chead)    
    
    fieldName = strsplit(Chead{c}{1}, ' ');                                                
    p.(fieldName{1}) = p.(fieldName{1})(stimtype);
    
end

for i = 1:numel(p.intensity)
    intensity_num(i, :) = eval(p.intensity{i});
end

p.intensity = intensity_num(:, 2);

end
