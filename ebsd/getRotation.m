function rot = getRotation(angfile)
%  GETROTATION   Uses the precomputed defect centers to determine a
%  best-fit connecting line, and then computes the angle with the x-axis. 
%  Angle given in radians

trnsX = angfile.Defect_Centers(:,1)-angfile.XY_Translation(1);
trnsY = angfile.Defect_Centers(:,2)-angfile.XY_Translation(2);
defectsFit = polyfit(trnsX,trnsY,1); % yields m, b from y=mx+b
rot = acos(1/sqrt(1+defectsFit(1)));