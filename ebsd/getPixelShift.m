function [pixelShift,cornerNum] = getPixelShift(imageSize,defectList,corner)
% GETPIXELSHIFT   Calculate the distance from each defect to image corners
%       and select the defect that results in the minimum total shift. This
%       trusts the user to have performed the EBSD with essentially similar
%       scans, as all subsequent ANG files will use the same corner
% PIXELSHIFT - an [i,j] matrix indicating the pixel shifts to apply to this
%       particular image in order to place one defect on one corner

nDefects = size(defectList,1);
if exist('corner','var') % corner to use already supplied
    distList = zeros(nDefects,1);
    shiftList = zeros(nDefects,2,1);
    for i = 1:nDefects
        switch corner
            case 1 % origin / top left
                shiftList(i,1,corner) = defectList(i,1);
                shiftList(i,2,corner) = defectList(i,2);
            case 2 % top right
                shiftList(i,1,corner) = (defectList(i,1)-imageSize(1));
                shiftList(i,2,corner) = defectList(i,2);
            case 3 % bot right
                shiftList(i,1,corner) = (defectList(i,1)-imageSize(1));
                shiftList(i,2,corner) = (defectList(i,2)-imageSize(2));
            case 4 % bot left
                shiftList(i,1,corner) = defectList(i,1);
                shiftList(i,2,corner) = (defectList(i,2)-imageSize(2));
        end
        distList(i)=sqrt(shiftList(i,1,corner)^2+shiftList(i,2,corner)^2);
    end
    cornerNum = corner;
    minDistRow = find(distList == min(distList));
else
    % Build general distance list
    distList = zeros(nDefects,4); % 4 corners = 4 distances
    shiftList = zeros(nDefects,2,4);
    for i = 1:nDefects
        for j = 1:4
            switch j
                case 1 % origin / top left
                    shiftList(i,1,j) = defectList(i,1);
                    shiftList(i,2,j) = defectList(i,2);
                case 2 % top right
                    shiftList(i,1,j) = (defectList(i,1)-imageSize(1));
                    shiftList(i,2,j) = defectList(i,2);
                case 3 % bot right
                    shiftList(i,1,j) = (defectList(i,1)-imageSize(1));
                    shiftList(i,2,j) = (defectList(i,2)-imageSize(2));
                case 4 % bot left
                    shiftList(i,1,j) = defectList(i,1);
                    shiftList(i,2,j) = (defectList(i,2)-imageSize(2));
            end
            distList(i,j)=sqrt(shiftList(i,1,j)^2+shiftList(i,2,j)^2);  
        end
    end
    % Find minimum row & column
    [minDistRow,cornerNum] = find(distList == min(min(distList)));
end



pixelShift = [shiftList(minDistRow,1,cornerNum), shiftList(minDistRow,2,cornerNum)];
