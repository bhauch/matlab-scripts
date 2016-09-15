function xydata = angfile2xydata(angfile,selector)
%  ANGFILE2XYDATA     convert the by-row representation of a TSL EBSD scan
%  into an image representation for the Image Quality and the Confidence
%  Index datasets. If CI, the raw data is thresholded using the current
%  value in the angfile structure
% ANGFILE - STRUCT: the .ANG file exported from TSL OIM Analysis
% SELECTOR - STRING: either 'IQ' or 'CI' 
  xx = angfile.FileData(:,4); % Column vector of X values
  yy = angfile.FileData(:,5); % Column vector of Y values
    
% X and Y values in the ANG file are not guaranteed to be integers, but 
% image representation requires row,column assignment, e.g. integer 
% subscripts. Divide by the ANG file's self-reported step size, convert
% to integer, and add 1 (because 0,0 isn't a valid matrix index position)
  xxx = 1+int64(xx ./ str2double(strrep(angfile.FileHeader{27},'# XSTEP: ','')));
  yyy = 1+int64(yy ./ str2double(strrep(angfile.FileHeader{28},'# YSTEP: ','')));
% Check for duplicate indices in the [yyy xxx] matrix by creating a test
% matrix that should have 1 at every index If a duplication occurs, there
% will be a 2 and a 0 somewhere (or 3 and two 0's... etc). However, 0's can
% occur from padding of the matrix (e.g. noncontiguous or nonrectangular)
  if (sum(sum(accumarray([yyy xxx],1)>1)) ~= 0)
    disp(['WARNING - ' angfile.FileName ' has duplicated XY indices']);
  end
  
  if (strcmpi(selector,'iq'))  
    selectedColumn = 6;
  elseif (strcmpi(selector,'ci'))
    selectedColumn = 7;
  else
    error('myFn:WrongInput','must enter string ''IQ'' or ''CI'' as the second input to angfile2xydata')
  end
  xydata = xydata2arr(xxx, yyy, angfile.FileData(:,selectedColumn));
end