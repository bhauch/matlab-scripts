function rot = getRotation(angfile)
%  GETROTATION   Uses the precomputed defect centers to determine a
%  best-fit connecting line, and then computes the angle with the x-axis. 
%  Angle given in radians

trnsX = angfile.Defect_Centers(angfile.Defect_Centers(:,1)>0,1);
trnsY = size(angfile.IQimage,2)-angfile.Defect_Centers(angfile.Defect_Centers(:,1)>0,2); % image y is "inverted"
defectsFit = polyfit(trnsX,trnsY,1); % yields m, b from y=mx+b
rot = acos(1/sqrt(1+defectsFit(1)^2)); % from dot product
