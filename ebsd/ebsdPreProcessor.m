%% Reconstruction Preprocessor
% Originally written by Benjamin Hauch [bhauch at gmail com]

%{
This script operates on a set of input ANG files to perform an alignment
using defects in the scan, which are assumed to be fiducial markers of
sorts, either microhardness indents or FIB-milled circles.

The general prescription: 
* For each slice: - Read in ANG file to a data
  structure 
    - Identify (X,Y) pairs associated with threshold CI < 0.01 
    - Search identified pairs for circle geometries, and sub-select only 
      those x,y that fall into circular defect clusters 
    - Compute the centers of these clusters 
    - Compute x,y translation to move the scan origin to the center of the
      closest defect cluster & store 
    - Perform the x,y origin translation on all points in that ANG file

* Choose an ANG file to serve as the reference orientation 
* For other ANG files 
    - Compute the rotational transformation needed to make the line between
      first & last defect cluster centers parallel/co-linear with that of 
      the reference slice (could look at other inter-defect lines too) 
    - Perform transformation on the ANG file

* For each slice: 
    - Perform inverse x,y origin translation to remove negative x,y (unless 
      DREAM3D can handle them) 
    - Write new ANG file with prefix to prevent overwriting data inputs

* Write txt file indicating transformations applied to each slice
%}
clear;clc;
MASK_THRESHOLD = 0.01;
FIDUCIAL_DIAMETER = 2;

%% Read in ANG files

bFiles_Selected = false;
% Ask user to select files
while (bFiles_Selected == false)
    angFolder = uigetdir('','Select the directory containing the .ANG files');
    if (angFolder == 0);
       disp('Folder selection aborted. Exiting Preprocessor');
       return;
    else
        fileList = dir(angFolder);
        fileList=fileList(~cellfun('isempty',{fileList.date})); % Strip invalid entries

%{ 
angFiles: struct of info/properties
Column List:
    FileName - name of the file
    FileData - the file's data
    CI_Threshold - custom threshold used for isolating defects
    Defect_Diam - Diameter used for defect-finding
    Defect_Centers - centers of found defects
    XY_Translation - Coordinate shift used to make one defect (0,0)
    isTranslated - whether the file data is currently translated'
    RotTransform - the rotational transformation used to align the line
                    between the first and last defects' centers
    isRotated - whether the file data is currently rotated
%}
        j = numel(fileList);
        hasANGfiles = false;
        for i = numel(fileList):-1:1;
            if regexp(fileList(i).name,'\.[Aa][Nn][Gg]')
                hasANGfiles=true;
                angFiles(j).FileName=fileList(i).name;
                angFiles(j).FileHeader=[];
                angFiles(j).FileData=[];
                angFiles(j).CI_Threshold = MASK_THRESHOLD;
                angFiles(j).Defect_Diam = FIDUCIAL_DIAMETER;
                angFiles(j).Defect_Centers=[];
                angFiles(j).XY_Translation = [];
                angFiles(j).isTranslated = false;
                angFiles(j).RotTransform = [];
                angFiles(j).isRotated = false;
                j = j - 1;
            end
        end
        if (hasANGfiles)
            angFiles=angFiles(~arrayfun(@(s) all(structfun(@isempty,s)), angFiles));
            % Display files selected for preprocessing
            disp(['.ANG files present in ' angFolder]);
            disp({angFiles.FileName}');
            user_response = '';
            while (strcmpi(user_response,''))
                user_response=questdlg('Use these files?');
                if (strcmpi(user_response,'cancel'))
                    disp('Exiting preprocessor...');
                    clear
                    return;
                end
            end
            if (strcmpi(user_response,'yes'))
                bFiles_Selected = true;
            else
                clear fileList angFolder angFiles i j
            end
        else
            disp(['No ANG files in ' angFolder])
        end
    end
end
clear i j user_response hasANGfiles fileList bFiles_Selected
%% File Loop

% Load data from the ANG files now that the user has confirmed the
% selection of files to use
% j1 | F | j2 | x | y | IQ | CI | PhaseID | DetectI | Fit
ss=get(0,'screensize');
for i = 1:numel(angFiles)
    rawfile=importdata(strcat(angFolder,'\\',angFiles(i).FileName),' ',38);
    angFiles(i).FileHeader = rawfile.textdata;
    angFiles(i).FileData=rawfile.data;
    clear rawfile;
    
    %{ 
    Convert X, Y, IQ vectors into rectangular matrix for image generation
          --> ANG files loop x for a given y, then increment y <--
    X and Y values in the ANG file are not guaranteed to be integers, but 
    image representation requires row,column assignment, e.g. integer 
    subscripts. Divide by the ANG file's self-reported step size, convert
    to integer, and add 1 (because 0,0 isn't a valid matrix index position)
    %} 
    xx = angFiles(i).FileData(:,4);
    yy = angFiles(i).FileData(:,5);
    angFiles(i).xxx = 1+int64(xx ./ str2double(strrep(angFiles(i).FileHeader{27},'# XSTEP: ','')));
    angFiles(i).yyy = 1+int64(yy ./ str2double(strrep(angFiles(i).FileHeader{28},'# YSTEP: ','')));
    % Check for duplicate indices in the [yyy xxx] matrix
    if (sum(sum(accumarray([angFiles(i).yyy angFiles(i).xxx],1)>1)) ~= 0)
        % All indices should be 1 if no duplicates, so the sum should
        % equal the number. If a duplication occurs, there will be a 2 and
        % a 0 somewhere (or 3 and two 0's... etc). However, 0's can occur
        % from padding of the matrix (e.g. noncontiguous data or
        % nonrectangular data
        disp(['WARNING - ' angFiles(i).FileName ' has duplicated XY indices']);
        % only a warning is issued at the moment
    end
    % Assemble the array of Image Quality data for display as an image
    IQArray=accumarray([angFiles(i).yyy angFiles(i).xxx],angFiles(i).FileData(:,6));
    angFiles(i).IQimage=mat2gray(IQArray);

    % Render side-by-side plot of IQ and the thresholded pixels
    fheight=size(IQArray,2)+10; 
    fwidth=2*size(IQArray,1)+40;
    fh=figure('name',['Thresholded ' angFiles(i).FileName], ...
        'outerposition',[10 ... Distance to left screen border
                         ss(4)-10-fheight ... Distance from bottom
                         fwidth fheight], ...
        'menubar','none', ...
        'numbertitle','off', ...
        'toolbar','none', ...
        'visible','off');
    colormap('gray');
    iqplot=subplot(1,2,1);
    imshow(angFiles(i).IQimage);title('IQ Map');
    ciplot=subplot(1,2,2,'DataAspectRatio',[1 1 1]);
        
    iterate_threshold=true;
    while iterate_threshold
        % Create ANG-specific mask for CI < $MASK_THRESHOLD
        cipixels=1-(angFiles(i).FileData(:,7)<angFiles(i).CI_Threshold);
        scatter(ciplot,xx,yy,1,cipixels,'filled');
        title(strcat('CI Threshold = ',num2str(angFiles(i).CI_Threshold)));
        xlabel('µm');ylabel('µm');
        ciplot.YDir='reverse';
        figure(fh) % Make visible
        
        % Ask user to confirm threshold value
        user_response = '';
        while (strcmpi(user_response,''))
            user_response=questdlg(strcat('Use threshold value =',num2str(angFiles(i).CI_Threshold)));
        end
        if (strcmpi(user_response,'yes'))
            iterate_threshold = false;
        else
            % ask for new CI_Threshold value
            need_new_threshold=true;
            while need_new_threshold
                user_response=input(['Old threshold was ' num2str(angFiles(i).CI_Threshold) '\nEnter new CI threshold value: '],'s');
                user_response = str2double(user_response);
                if (user_response>=0 && user_response < 1)
                    need_new_threshold=false;
                    disp(['Repeating with threshold value = ' num2str(user_response)]);
                    angFiles(i).CI_Threshold=user_response;
                    clear user_response;
                else
                    disp('Invalid response. Threshold must be numeric, between 0 and 1')
                end
            end
            clear need_new_threshold
        end
    end
    
    % User has iteratively chosen the appropriate CI to use. Convert the
    % cipixels vector into an image matrix as was done for IQimage
    CIArray=accumarray([angFiles(i).yyy angFiles(i).xxx],cipixels);
    angFiles(i).CIimage=mat2gray(CIArray);
    subplot(iqplot) % overwrite IQmap with the CI image based on the threshold
    imshow(angFiles(i).CIimage);title('Original CI image');
    clear cipixels xx yy iterate_threshold
    
    % Use pixel erosion and growth to remove isolated low-CI pixels 
    refine_selection = true; continue_refine = false;
    refine_operator_radius = 1;
    while refine_selection
        if continue_refine == false;
            angFiles(i).dilatorHistory=num2str(refine_operator_radius);
            modCIimage=angFiles(i).CIimage; % load from original thresholded CI image
        else
            angFiles(i).dilatorHistory=strcat(angFiles(i).dilatorHistory,'>',num2str(refine_operator_radius));
        end
        disp(['Using disk of radius ' num2str(refine_operator_radius)])
        dilator=strel('disk',refine_operator_radius);
        modCIimage=imopen(modCIimage,dilator); 
        modCIimage=imclose(modCIimage,dilator);
        
        % Show modified image (dilate = increase white, erode = increase black)
        figure(fh); subplot(ciplot);
        imshow(modCIimage);title(['Cleaned w/radius ' angFiles(i).dilatorHistory]);
        
        % Ask user for input
        user_response = '';
        while (strcmpi(user_response,''))
            user_response=questdlg(['End cleaning step for ' angFiles(i).FileName '?'], ...
                            angFiles(i).FileName, ... prompt title
                            'Yes', 'Do More', 'Start Over', ...
                            'Yes'); % Default to Yes
        end
        if (strcmpi(user_response,'yes')) % user is satisfied with this result
            refine_selection = false; continue_refine = false;
            angFiles(i).modCIimage=modCIimage;
            clear modCIimage dilator refine_operator_radius user_response
        else
            % ask for new disk radius value
            refine_operator_radius=getNewRadius(refine_operator_radius);
            if (strcmpi(user_response,'start over'))
                % User wants to restart the refinement process from the
                % original image... nothing to do here
                continue_refine = false;
                disp(['Processing thresholded ' angFiles(i).FileName ...
                    ' CI image starting with disk radius = ' ...
                    num2str(refine_operator_radius)]);
            else
                % User wants to continue refining the current modified
                % image
                continue_refine = true;
            end
        end
    end
    clear continue_refine refine_selection
    
    % Analyze modified CI image for contiguous black regions of size
    % .Defect_Diam and report centers in .Defect_Centers[x y]
    % Analyze coordinates for circular regions of diameter .Defect_Diam
    angFiles(i).defects=bwboundaries(1-angFiles(i).modCIimage,'noholes');
    centers=zeros(size(angFiles(i).defects,1),2);
    for j=1:size(angFiles(i).defects,1)
        xyr=CircleFitByPratt(angFiles(i).defects{j});
        if xyr(1,3) >= angFiles(i).Defect_Diam
            centers(j,:) = [xyr(1,1) xyr(1,2) sqrt(xyr(1,1)^2+xyr(1,2)^2)];
        end
    end
    disp(['Found ' num2str(sum(centers(:,1)>1)) ' markers in ' ...
            angFiles(i).FileName]);
    angFiles(i).Defect_Centers=centers; % x | y | distance from origin
    minDistRow=find(centers==min(centers));
    angFiles(i).XY_Translation=[centers(minDistRow,1) centers(minDistRow,2)];
    clear xyr centers
    
    % Perform the X-Y transpose by shifting all points in accordance with
    % the XY center of the defect with the shortest distance
    
    
    break % remove to iterate all angFiles
end
disp('EBSD preprocessing has completed');
clear ss;