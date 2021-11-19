function result = pull_gsheet(DOCID, GID)
% result = pull_gsheet(DOCID)
% Download a google spreadsheet as csv and import into a Matlab cell array.
%
% [DOCID] see the value after 'key=' in your spreadsheet's url
%           e.g. '0AmQ013fj5234gSXFAWLK1REgwRW02hsd3c'
% [GID] value after 'gid=' in url - optional
%
% [result] cell array of the the values in the spreadsheet
%
% IMPORTANT: The spreadsheet must be shared with the "anyone with the link" option
%
% This has no error handling and has not been extensively tested.
% Please report issues on Matlab FX.
%
% DM, Jan 2013

loginURL = 'https://www.google.com';

if exist('GID', 'var') && ~isempty(GID)
    csvURL = ['https://docs.google.com/spreadsheet/ccc?key=', ...
        DOCID '&gid=' GID '&output=csv&pref=2'];
else
    csvURL = ['https://docs.google.com/spreadsheet/ccc?key=', ...
        DOCID '&output=csv&pref=2'];
end

%Step 1: go to google.com to collect some cookies
cookieManager = java.net.CookieManager([], java.net.CookiePolicy.ACCEPT_ALL);
java.net.CookieHandler.setDefault(cookieManager);
handler = sun.net.www.protocol.https.Handler;
connection = java.net.URL([],loginURL,handler).openConnection();
connection.getInputStream();

%Step 2: go to the spreadsheet export url and download the csv
connection2 = java.net.URL([],csvURL,handler).openConnection();
result = connection2.getInputStream();
result = char(readstream(result));

%Step 3: convert the csv to a cell array
result = parseCsv(result);

% if you are only dealing with numerical data
% result = cellfun(@str2num, result);

end

function data = parseCsv(data)
% splits data into individual lines
data = textscan(data,'%s','whitespace','\n');
% data = split(data, sprintf('\r\n')) % to test

data = data{1};

% trim any empty last line
% if isempty(data{end})
%     data = data(1:end-1);
% end

for ii=1:length(data)
   %for each line, split the string into its comma-delimited units
   %the '%q' format deals with the "quoting" convention appropriately.
   tmp = textscan(data{ii},'%q','delimiter',',');
   data(ii,1:length(tmp{1})) = tmp{1};
end

end

function out = readstream(inStream)
%READSTREAM Read all bytes from stream to uint8
%From: http://stackoverflow.com/a/1323535

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
byteStream = java.io.ByteArrayOutputStream();
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
isc.copyStream(inStream, byteStream);
inStream.close();
byteStream.close();
out = typecast(byteStream.toByteArray', 'uint8'); 

end
