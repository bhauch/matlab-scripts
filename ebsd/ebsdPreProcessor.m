%% Reconstruction Preprocessor
% Originally written by Benjamin Hauch [bhauch at gmail com]

%{
This script operates on a set of input ANG files to perform an alignment
using defects in the scan, which are assumed to be a linear set of fiducial
markers of sorts, either microhardness indents or FIB-milled circles

The general prescription: 
* For each slice: - Read in ANG file to a data
  structure 
    - Identify (X,Y) pairs associated with threshold CI < 0.01 
    - Search identified pairs for circle geometries, and sub-select only 
      those x,y that fall into circular defect clusters 
    - Compute the centers of these clusters 
    - Compute x,y translation to move the scan origin to the center of the
      closest defect cluster (e.g. world (x,y) --> 0,0 at the defect)
    - Compute the angle of the best-fit line through the centers, relative
      to the x-axis
    - Perform transforms and rotations on the image reference object
    - Extract shifted coordinates from image reference object
    - Write new ANG file
        - Perform additional X,Y shift if DREAM3D cannot handle negative XY
        in its ANG importer
        - Prefix the files and place in a subdirectory
        - Include "history" txt file for rapid reprocessing
%}
clear;clc;
%% Defaults/Repeat Check
% Check for an existing preprocessor history file in the active directory
% that contains instructions for parameters and/or sequences
% Instruction file extension is .preproc
if exist('*.preproc','file') > 0
    INTERACTIVE = false;
else
    INTERACTIVE = true;
    INIT_MASK_THRESHOLD = 0.01; % CI < threshold indicates possible marker
    MIN_FIDUCIAL_DIAMETER = 4;  % number of pixels wide the markers must be
    INIT_CI_DILATOR = 2;
    OUTPUT_DIR = 'aligned';     % will be added to current directory
end
[status,message] = mkdir(OUTPUT_DIR);
if (status == 0)
    error(['Unable to make output directory ' OUTPUT_DIR]);
elseif strcmpi(message,'') == false
    disp(['Output dir ' OUTPUT_DIR ' already exists, files may be overwritten if there are name conflicts.']);
    clear status message
end
%% Read in ANG files
if (~INTERACTIVE)
    % Grab the listed files from the instruction file, the per-file
    % settings (if more than one line, else apply given line to all)
    
else % Ask user to select files
    bFiles_Selected = false;
    while (bFiles_Selected == false)
        angFolder = uigetdir('','Select the directory containing the .ANG files');
        if (angFolder == 0);
           disp('Folder selection aborted. Exiting Preprocessor');
           return;
        else
            fileList = dir(angFolder);
            fileList = fileList(~cellfun('isempty',{fileList.date})); % Strip invalid entries

    %{ 
    angFiles: struct of info/properties
    Column List:
        FileName - name of the file
        FileData - the file's data
        CI_Threshold - custom threshold used for isolating defects
        Defect_Diam - Threshold Diameter used for defect-finding
        Defect_Centers - centers of found defects
        XY_Translation - Coordinate shift used to make one defect center =(0,0)
        RotAngle - the rotation angle for the best-fit line through the defects
    %}
            j = numel(fileList);
            hasANGfiles = false;
            for i = numel(fileList):-1:1;
                if regexp(fileList(i).name,'\.[Aa][Nn][Gg]')
                    hasANGfiles = true;
                    aF(j).FileName = fileList(i).name;
                    aF(j).FileHeader = [];
                    aF(j).FileData = [];
                    aF(j).CI_Threshold = INIT_MASK_THRESHOLD;
                    aF(j).Defect_Diam = MIN_FIDUCIAL_DIAMETER;
                    aF(j).Defect_Centers = [];
                    aF(j).XY_PixelShift = [];
                    aF(j).RotAngle = 0;
                    j = j - 1;
                end
            end
            if (hasANGfiles)
                aF = aF(~arrayfun(@(s) all(structfun(@isempty,s)), aF));
                % Display files selected for preprocessing
                disp(['.ANG files present in ' angFolder]);
                disp({aF.FileName}');
                user_response = '';
                while (strcmpi(user_response,''))
                    user_response = questdlg('Use these files?');
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
end

%% Slice Processing

% Load data from the ANG files now that the user has confirmed the
% selection of files to use
% j1 | F | j2 | x | y | IQ | CI | PhaseID | DetectI | Fit
ss=get(0,'screensize');
for i = 1:numel(aF)
    rawfile = importdata(strcat(angFolder,'\\',aF(i).FileName),' ',38);
    aF(i).FileHeader = rawfile.textdata;
    aF(i).FileData = rawfile.data;
    clear rawfile
    
    % Obtain X-by-Y matrix of IQ values for image rendering
    IQArray = angfile2xydata(aF(i),'iq');
    aF(i).IQimage = mat2gray(IQArray);

    % Generate spatial reference object 
    aF(i).imgRef = imref2d(size(IQArray), ...
             str2double(strrep(aF(i).FileHeader{27},'# XSTEP: ','')), ...
             str2double(strrep(aF(i).FileHeader{28},'# YSTEP: ','')));
    
    % Render side-by-side plot of IQ and the thresholded pixels
    fheight = size(IQArray,2)+10; 
    fwidth = 2*size(IQArray,1)+40;
    fh = figure('name',['Thresholded ' aF(i).FileName], ...
        'outerposition',[10 ... Distance to left screen border
                         ss(4)-10-fheight ... Distance from bottom
                         fwidth fheight], ...
        'menubar','none', ...
        'numbertitle','off', ...
        'toolbar','none', ...
        'visible','off');
    colormap('gray');
    iqplot = subplot(1,2,1);
    imshow(aF(i).IQimage);title('IQ Map');
    ciplot = subplot(1,2,2,'DataAspectRatio',[1 1 1]);
        
    iterate_threshold = true;
    while iterate_threshold
        % Create ANG-specific mask for CI < $MASK_THRESHOLD
        cipixels = 1-(aF(i).FileData(:,7)<aF(i).CI_Threshold);
        scatter(ciplot,aF(i).FileData(:,4),aF(i).FileData(:,5),1,cipixels,'filled');
        title(strcat('CI Threshold = ',num2str(aF(i).CI_Threshold)));
        xlabel('µm'); ylabel('µm');
        ciplot.YDir = 'reverse';
        figure(fh) % Make visible
        
        % Ask user to confirm threshold value
        user_response = '';
        while (strcmpi(user_response,''))
            user_response=questdlg(strcat('Use threshold value = ',num2str(aF(i).CI_Threshold)));
        end
        if (strcmpi(user_response,'yes'))
            iterate_threshold = false;
        else
            % ask for new CI_Threshold value
            need_new_threshold = true;
            while need_new_threshold
                user_response = input(['Old threshold was ' num2str(aF(i).CI_Threshold) '\nEnter new CI threshold value: '],'s');
                user_response = str2double(user_response);
                if (user_response >= 0 && user_response < 1)
                    need_new_threshold = false;
                    disp(['Repeating with threshold value = ' num2str(user_response)]);
                    aF(i).CI_Threshold = user_response;
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
    CIArray = angfile2xydata(aF(i),'ci');
    aF(i).RawCIimage = mat2gray(CIArray);
    aF(i).CIimage = mat2gray(1-(CIArray<aF(i).CI_Threshold));
    subplot(iqplot)         % overwrite IQmap with the thresholded CI image
    imshow(aF(i).CIimage); title('Original Thresholded CI image');
    clear cipixels iterate_threshold CIArray fheight fwidth IQArray
    
    % Use pixel erosion and growth to remove isolated low-CI pixels 
    refine_selection = true; continue_refine = false;
    refine_operator_radius = INIT_CI_DILATOR;
    while refine_selection
        if continue_refine == false;
            aF(i).dilatorHistory = num2str(refine_operator_radius);
            modCIimage = aF(i).CIimage;     % load original thresholded CI image
        else
            aF(i).dilatorHistory = strcat(aF(i).dilatorHistory,'->',num2str(refine_operator_radius));
        end
        disp(['Using disk of radius ' num2str(refine_operator_radius)])
        dilator = strel('disk',refine_operator_radius);
        modCIimage = imopen(modCIimage,dilator); 
        modCIimage = imclose(modCIimage,dilator);
        
        % Show modified image (dilate = increase white, erode = increase black)
        figure(fh); subplot(ciplot);
        imshow(modCIimage); title(['Cleaned w/radius ' aF(i).dilatorHistory]);
        
        % Ask user for input
        user_response = '';
        while (strcmpi(user_response,''))
            user_response=questdlg(['End cleaning step for ' aF(i).FileName '?'], ...
                            aF(i).FileName, ... prompt title
                            'Yes', 'Do More', 'Start Over', ...
                            'Yes'); % Default to Yes
        end
        if (strcmpi(user_response,'yes')) % user is satisfied with this result
            refine_selection = false; continue_refine = false;
            aF(i).modCIimage = modCIimage;
            clear modCIimage dilator refine_operator_radius user_response
        else
            % ask for new disk radius value
            refine_operator_radius = getNewRadius(refine_operator_radius);
            if (strcmpi(user_response,'start over'))
                % User wants to restart the refinement process from the
                % original image... nothing to do here
                continue_refine = false;
                disp(['Processing thresholded ' aF(i).FileName ...
                    ' CI image starting with disk radius = ' ...
                    num2str(refine_operator_radius)]);
            else
                % User wants to continue refining the current modified
                % image
                continue_refine = true;
            end
        end
    end
    clear continue_refine refine_selection ciplot iqplot fh
    
    % Analyze modified CI image for contiguous black regions of size
    % .Defect_Diam and report centers, radius in .Defect_Centers[x y]
    % Analyze coordinates for circular regions of min diameter .Defect_Diam
    aF(i).defects = bwboundaries(1-aF(i).modCIimage,'noholes');
    centers = zeros(size(aF(i).defects,1),3);
    for j = 1:size(aF(i).defects,1)
        xyr = CircleFitByPratt(aF(i).defects{j});
        if ((xyr(1,3)*2) >= aF(i).Defect_Diam)
            centers(j,:) = [xyr(1,2) xyr(1,1) sqrt(xyr(1,1)^2+xyr(1,2)^2)];
        end % images have inverted x,y referencing, hence the 1,2 then 1,1 above
    end
    disp(['Found ' num2str(sum(centers(:,1)>1)) ' markers in ' aF(i).FileName]);
    aF(i).Defect_Centers = centers; % x | y | distance from origin
    clear xyr j centers
    
    % Calculate the XY translation based on the nearest corner to a defect
    % center. imtranslate assumes world-coord shift, so account for pixel
    % pixel size in creating the transform
    if (~exist('SHIFT_CORNER','var'))
        [aF(i).XY_PixelShift, SHIFT_CORNER] = ... 
            getPixelShift(size(aF(i).IQimage),aF(i).Defect_Centers);
    else
        % Shift corner is already picked (e.g. not the first EBSD slice)
        aF(i).XY_PixelShift = ...
            getPixelShift(size(aF(i).IQimage),aF(i).Defect_Centers,SHIFT_CORNER);
    end
    aF(i).XY_Translation = [aF(i).XY_PixelShift(1)*aF(i).imgRef.PixelExtentInWorldX, ...
                            aF(i).XY_PixelShift(2)*aF(i).imgRef.PixelExtentInWorldY];
    
    % Determine the rotational relationship between the datafile and the
    % x-axis and store the rotation angle
    aF(i).RotAngle = getRotation(aF(i));
    
    % Perform image translation and rotation about the corner defect
    aF(i).RT_Transform = getRotationTransform(aF(i).RotAngle,aF(i).XY_Translation(1),aF(i).XY_Translation(2));
    [aF(i).RT_IQ, aF(i).imgRefRT] = imwarp(aF(i).IQimage, aF(i).imgRef, aF(i).RT_Transform);
    figure('numbertitle','off','name',['R+T ' aF(i).FileName]); 
    imshow(aF(i).RT_IQ);        % display result   
    
    % Generate new X & Y columns for output file by passing through the
    % transform (needs work yet, missing initial translate)
    [aF(i).outputX, aF(i).outputY] = aF(i).RT_Transform.transformPointsForward(aF(i).FileData(:,4),aF(i).FileData(:,5));
    aF(i).outputName = ['procd_' aF(i).FileName];
    
    if i>1 ;break;end; % remove to iterate all angFiles
end
%% ANG output
% Assemble & write the ANG files with new XY coordinates. IQ/CI are
% unchanged

for i = 1:numel(aF)
    if size(aF(i).outputX,1) > 0
        outputHeader = aF(i).FileHeader;
        outputData = aF(i).FileData;
        outputData(:,4)=aF(i).outputX;
        outputData(:,5)=aF(i).outputY;
        outputFileName = ['.\\' OUTPUT_DIR '\\' aF(i).outputName];
        fd = fopen(outputFileName,'w');                         % Open file
        fprintf(fd, '%s\n',outputHeader{1:size(outputHeader)}); % Write header
        fclose(fd);                                             % Close file
        dlmwrite(outputFileName,outputData,'-append','delimiter','\t'); % Data
    end
end
clear outputHeader outputData outputFileName fd
%% History Output-to-file
% Write history file into the output dir
% Since output of ang files results in direct overwrite given naming
% conflict, should attempt to not overwrite existing history files so that
% past processing can be recovered

% save the input file listing, ci threshold per file, minimum defect size,
% defect centers (pixel coords!), the corner for shifting, the rotation
% angle per file, and the dilator history

%historyOutput = getHistoryOutput(angFiles);
outputName = ['.\\' OUTPUT_DIR '\\history_' date '.pphist'];
fileExists = exist(outputName,'file');
while fileExists ~= 0
    % mutate history file's name until a unique is available
    outputName = ['.\\' OUTPUT_DIR '\\history_' date '_' ...
        datestr(now,'HH.MM.SS') '.pphist'];
    fileExists = exist(outputName,'file');
end

% Write history file
%fd = fopen(outputName,'w');
%fprintf(fd,'%s\n',['EBSD preprocessing completed: ' datestr(now)]);
%fprintf(fd,'%s\n',historyOutput{1:size(historyOutput)});
%fclose(fd);

disp('EBSD preprocessing has completed');
clear ss fileExists outputName INTERACTIVE;
