function [xClip,yClip] = clipdata2
% clipdata.m
% 20060924- update to to bugs with new matlab edition
% 13-Oct-2005 JS: script written
% pulls xy data between 2 datatips from the current figure
% a figure must already be open with 2 datatips on a line
% see 'doc datacursormode' for more info on matlab datatips
%
% the line to be clipped should either contain the datatips or be the 
% current object
%
% TODO: interactive datatip placement, automated default datatips
%
% USAGE:
% 1) select the line that you would like to extract data from
% 2) if only a portion of that line is desired, use data tips
% 3) syntax [xclip,yclip] = clipdata;
fig = gcf;
cursorInfo = getCursorInfo(datacursormode(fig));
nDataTip = size(cursorInfo,2);
xDataTip = [];
for iDataTip = 1:nDataTip
    xDataTip(iDataTip) = cursorInfo(iDataTip).Position(1);
end
if strcmp(get(gco,'Type'),'line')
    disp('clipdata: clipping selected line');
    xData = get(gco,'XData');
    yData = get(gco,'YData');
else
    disp('clipdata: clipping line with datatip #1');
    xData = get(cursorInfo(1).Target,'XData');
    yData = get(cursorInfo(1).Target,'YData');
end
if xDataTip
    iClipMin = find(xData==min(xDataTip));
    iClipMax = find(xData==max(xDataTip));
    xClip = xData(iClipMin:iClipMax);
    yClip = yData(iClipMin:iClipMax);
else
    xClip = xData;
    yClip = yData;
end