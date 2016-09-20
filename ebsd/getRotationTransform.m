function tform = getRotationTransform(rotationAngle, rcX, rcY)
% GETROTATIONTRANSFORM   return the transformation matrix for a given
% rotation about a specified point, passed through affine2d in order to generate a tform
% struct
% ROTATIONANGLE - 2D rotation angle in radians
% RCX           - x-value for the point to rotate about (world coords)
% RCY           - y-value for the point to rotate about (world coords)

shiftTo = [1 0 0; 0 1 0; -rcX -rcY 1;];
shiftFrom = [1 0 0; 0 1 0; rcX rcY 1;];
transMatrix = [cos(rotationAngle) -sin(rotationAngle) 0; ...
               sin(rotationAngle) cos(rotationAngle) 0; ...
               0 0 1];
tform = affine2d(shiftTo*transMatrix*shiftFrom);
