tic 
Array1 = []; % Empty matrix defintion for future time points
Array2 = []; % Empty matrix defintion for future time points
happy_test = 0; % this while loop will run until the user is happy with the correlation frame test (until [happy_test] = 1) Checkpoints thorught code.
first_time_point = 0;
Driftcorrevals = [];
Xcorrevals = []; % Opens empty matrix for writing in futur results.
Magmeanvals = [];
% Loop start for individual 3D matrices representing each individual time point.
for k = 0:3
    
    
    if k==0
        [filename, pathname]=uigetfile('*.tif', 'Image sequence'); % Only opens GUI on first loop and then sequentially opens.
        
        % Finding location of first number in filename.
        for w = 1 : numel(filename)
            if  (str2double(filename(w))==1)
                numpos = w;
                break
            end
        end
    end
    
    filewords = filename(1 : numpos-1); % Extracts filename up to number
    
    
    
    filename = strcat(filewords,num2str(k),'.tif'); % Allows filename to be used and k number to be added sequentially with each loop.
    info=imfinfo(strcat(pathname,filename));
    numframes=numel(info);
    clear frame; % clears frame from possible last time point
    frame =zeros(info(1).Height,info(1).Width,numframes); % All the information required in each z stack
    for icount = 1:numframes % string all z stacks into 1 3d matrix.
        frame(:,:,icount) = imread(strcat(pathname,filename),icount,'Info',info);
    end
    
    
    Array2 = frame(:,:,:); % Defines Array2
   % [xsize,ysize,zsize]=size(Array2); % Sizes in 3Dimensions
    
    
    if k == 0
        %first_time_point == 0
        Array1=Array2;
        disp("first time round")
   
    else
        [newxsize,newysize,newzsize]=size(Array2);
       [xsize,ysize,zsize]=size(Array1); 
        % Create Function that can be recalled back in one line
        % Function for Drift correction
        disp("xsize")
        disp(xsize)
        disp("ysize")
        disp(ysize)
        disp("zsize")
        disp(zsize)
        disp("size Array 1")
        disp(size(Array1))
        disp("size Array 2")
        disp(size(Array2))
        
      
        %     frame=0;
        %bring size of Array2 down to size of Array1
        updatex= newxsize-xsize;
        updatey = newysize-ysize;
        updatez = newzsize-zsize;
        
        Array2([1:updatex],:,:) = [];
        Array2(:,[1:updatey],:) = [];
        Array2(:,:,[1:updatez]) = [];
%         [xsize,ysize,zsize]=size(Array2);
%         [newxsize,newysize,newzsize]=size(Array1);
        %run drift correction function
        [True_drift_x,True_drift_y,True_drift_z] = run_drift_correction (xsize, ysize, zsize, Array1, Array2);
        
        True_drift_x = ceil(True_drift_x);
        True_drift_y =  ceil(True_drift_y);
        True_drift_z = ceil(True_drift_z);
        
        %A(:,:,[1:100,301:400]) = [];
        % Drift cutting to original ARRAY2
    
        if True_drift_x>=0
            Array1([1:True_drift_x],:,:) = [];
            Array2 ([(xsize-(True_drift_x-1)):xsize],:,:)= [];
        end
        if True_drift_x<0
            Array1([(xsize + (True_drift_x+1)):xsize],:,:) = [];
            Array2 ([1:((True_drift_x)*-1)],:,:)= [];
        end
        
        if True_drift_y>=0
            Array1(:,[1:True_drift_y],:) = [];
            Array2 (:,[(ysize-(True_drift_y-1)):ysize],:)= [];
        end
        if True_drift_y<0
            Array1(:,[(ysize + (True_drift_y+1)):ysize],:) = [];
            Array2 (:,[1:((True_drift_y)*-1)],:)= [];
        end
        
        if True_drift_z>=0
            Array1(:,:,[1:True_drift_z]) = [];
            Array2 (:,:,[(zsize-(True_drift_z-1)):zsize])= [];
        end
        if True_drift_z<0
            Array1(:,:,[(zsize+(True_drift_z+1)):zsize]) = [];
            Array2 (:,:,[1:((True_drift_z)*-1)])= [];
        end
    
    
        disp("Cut of Arrays")
        disp(size(Array1))
        disp(size(Array2))
        
        
        
        while happy_test == 0 % while happy test is 0 will allow you to modify parameters befare launching analysis.
            
            
            
            % Opens GUI to set parameters for source, search and grid size.
            prompt = {'Source Size XY [pxl]',... % The xy size in pxl to be anlysed
                'Source Size Z [pxl]'... % the z size in slice number to be analysed
                'Search Size XY [pxl]'... % the xy size in which source size will move to find correlation (must be equal or larger than source size)
                'Search Size Z [pxl]'... % the z size in slice number which source size will move to find correlation (equal or larger than source size)
                'Grid Size XY [pxl]'... % Grid xy size in pxl moving through Array
                'Grid Size Z [pxl]'... % Grid z size in slice moving through Array (might be able to replace with search or source)
                'XY pxl->nm'...
                'Z step [nm]'}...
                %             'Time difference';
            % GUI for inputs of parameter
            title = 'Parameters';
            dims = [1 35];
            definput = {'30', '5', '60', '10', '20', '2', '104', '261'}; % standard inputs, can be modified.
            user_answer = inputdlg(prompt,title,dims,definput);
            %Defines parameters from above GUI box
            source_size_xy = str2double(user_answer{1,1});
            source_size_z = str2double(user_answer{2,1});
            search_size_xy = str2double(user_answer{3,1});
            search_size_z = str2double(user_answer{4,1});
            grid_size_xy = str2double(user_answer{5,1});
            grid_size_z = str2double(user_answer{6,1});
            xy_nm = str2double(user_answer{7,1});
            z_step = str2double(user_answer{8,1});
            %         time = str2double(user_answer{9,1});
            % Creating 3D grid to move through in xyz for source and search size.
            maxi = ceil (size(Array2,1)/grid_size_xy);
            maxj = ceil (size(Array2,2)/grid_size_xy);
            maxl = ceil (size(Array2,3)/grid_size_z);
            % Creating empty grid of correct Grid size locations
            Xcorrx = zeros(maxi, maxj, maxl); % For X shift in pxl
            Xcorry = zeros(maxi, maxj, maxl); % For Y shift in pxl
            Xcorrz = zeros(maxi, maxj, maxl); % For Z shift in pxl
            % Creating empty grid of correct Grid size locations
            Xcorrxnm = zeros(maxi, maxj, maxl); % For X shift in nm
            Xcorrynm = zeros(maxi, maxj, maxl); % For Y shift in nm
            Xcorrznm = zeros(maxi, maxj, maxl); % For Z shift in nm
            % Creating empty grid for Magnitude in 3D
            Xcorrmag = zeros (maxi, maxj, maxl);
            Xcorrmagnm = zeros (maxi, maxj, maxl);
            % Asks if you want to visualise parameters on Maximum Intensity
            % Projection (MIP)
            run_test_q = questdlg('Do you want to visualise parameters?', ...
                'Run test',...
                'Yes', 'No', 'No');
            
            if strcmp(run_test_q, 'Yes') % if yes displays MIP
                run_test = 1
                
                if run_test == 1 % if answered yes to display of MIP
                    mip = max(Array2, [], 3); % Creates new matrix of MIP
                    imagesc(mip) % Displays MIP
                    axis on;
                    [rows,columns]=size(mip); % Accurate size of MIP
                    
                    hold on; % Holds image to overlay search and soucre size grid
                    for row = 1 : grid_size_xy : rows
                        line ([1, columns],[row,row], 'Color', 'g'); % Grid rows in Green for Source size
                    end
                    %                 for row = 1 : search_size_xy : rows
                    %                     line ([1, columns],[row,row], 'Color', 'r'); % Grid rows in Red for search size
                    %                 end
                    for col = 1 : grid_size_xy : columns
                        line([col,col], [1,rows], 'Color', 'g'); % Grid columns in Green for Source size
                    end
                    %                 for col = 1 : search_size_xy : columns
                    %                     line([col,col], [1,rows], 'Color', 'r');% Grid columns in red for Search size
                    %                 end
                end
            end
            
            % If you don't want to VIsualise parameters and select NO
            if strcmp(run_test_q,'No')
                run_test = 0;
                
                if run_test == 0
                    happy_test_q = questdlg('Are you happy with the parameters?', ...  % Asks if happy with parameters
                        'Happy?',...
                        'Yes', 'No', 'Yes');
                    
                    if strcmp(happy_test_q, 'Yes') % If YES launches analysis
                        happy_test = 1;
                    end
                    if strcmp(happy_test_q, 'No')% if NO returns to Parameter inputs
                        happy_test = 0;
                    end
                end
            end
        end
        
        %     disp ('thisisk')
        %     disp (k)
        % If k==0 (still on first time point) next phase (the analysis) is skipped as there is no second Array with which to analyse it against.
        if happy_test == 1 % If parameters have been set and accepted
            for i = 1:maxi % Loop moving across in X Grid Distance
                for j = 1:maxj % Loop moving across in Y Grid Distance
                    for l = 1:maxl % Loop moving across in Z Grid Distance
                        
                        %                         n = n+1
                        % Building an index based off Grid Size to apply to
                        % source and search volumes.
                        xindex = 1+(i-1)*grid_size_xy;
                        yindex = 1+(j-1)*grid_size_xy;
                        zindex = 1+(l-1)*grid_size_z;
                        
                        try
                            % 3D source volume moving through Grid location
                            sub_roi3D = Array1...
                                (xindex - ceil(source_size_xy):...
                                xindex + ceil(source_size_xy),...
                                yindex - ceil(source_size_xy):...
                                yindex + ceil(source_size_xy),...
                                zindex - ceil(source_size_z):...
                                zindex + ceil(source_size_z));
                        catch
                            continue
                        end
                        
                        try
                            % 3D Search volume for corresponding Source
                            % volume moving in equal Gid loaction
                            sub_area3D = Array2...
                                (xindex - ceil(search_size_xy):...
                                xindex + ceil(search_size_xy),...
                                yindex - ceil(search_size_xy):...
                                yindex + ceil(search_size_xy),...
                                zindex - ceil(search_size_z):...
                                zindex + ceil(search_size_z));
                            %                             disp ('catch2')
                            
                        catch
                            continue
                        end
                        
                        % Convolution Function in 3D from sub-ROI, SubArea
                        try
                            cor_3d=convn(sub_roi3D,sub_area3D(end:-1:1,end:-1:1,end:-1:1));
                            %                             f= f + 1
                            
                        catch
                            continue
                        end
                        
                        % Convolution offset for final pixel shift +/-
                        xsize_convoffset = size(sub_area3D,1) + size(sub_roi3D,1);
                        ysize_convoffset = size(sub_area3D,2) + size(sub_roi3D,2);
                        zsize_convoffset = size(sub_area3D,3) + size(sub_roi3D,3);
                        
                        % Finding Convolution maximum across xyz.
                        [valz,iz]=max(cor_3d,[],3);
                        [valy,iy]=max(valz,[],2);
                        [valx,ix]=max(valy,[],1);
                        Truex = ix - (xsize_convoffset/2);
                        
                        [valz,iz]=max(cor_3d,[],3);
                        [valx,ix]=max(valz,[],1);
                        [valy,iy]=max(valx,[],2);
                        Truey = iy - (ysize_convoffset/2);
                        
                        [valx,ix]=max(cor_3d,[],1);
                        [valy,iy]=max(valx,[],2);
                        [valz,iz]=max(valy,[],3);
                        Truez = iz - (zsize_convoffset/2);
                        
                        % if valz 0 then no shift in xyz.
                        if valz==0
                            Truex=0;
                            Truey=0;
                            Truez=0;
                        end
                        
                        % Table of Data for k=time point, i=x position, j=y
                        % position, l=z position, Truex= x pixel shift, Truey=y
                        % pixel shift, Truez=z pixel shift, valz= Max number in 3D
                        newvals = [k,i,j,l,Truex,Truey,Truez,valz];
                        Xcorrevals = vertcat(Xcorrevals,newvals);
                        
                        % Shift per grid point in pxl
                        Xcorrx(i,j,l) = Truex;
                        Xcorry(i,j,l)= Truey;
                        Xcorrz(i,j,l)= Truez;
                        % Shift per grid point in nmy
                        
                        Xcorrxnm(i,j,l) = Truex*xy_nm;
                        Xcorrynm(i,j,l) = Truey*xy_nm;
                        Xcorrznm(i,j,l) = Truez*z_step;
                        
                        % Magnitude for each grid point
                        Xcorrmag(i,j,l) = sqrt(Truex^2 + Truey^2 + Truez^2); % all in pixels so not to scale!
                        
                        Xcorrmagnm(i,j,l) = sqrt((Truex*xy_nm)^2 + (Truey*xy_nm)^2 + (Truez*z_step)^2);
                        
                        
                        
                        % h = h + 1
                    end
                end
            end
            
            %[Truex,Truey,Truez,valz] = run_PIV (maxi,maxj,maxl,Array1,Array2,grid_size_xy,grid_size_z,source_size_xy,source_size_z,search_size_xy,search_size_z);
            
            
            
            Driftvals = [k,xsize,ysize,zsize,True_drift_x,True_drift_y,True_drift_z];
            Driftcorrevals = vertcat(Driftcorrevals,Driftvals);
            
            
            %Average mag per min in microns
            %save in empty cell and append!
            Magmean = mean(Xcorrmagnm(:));
            RealMagMean = Magmean/(10*1000);
            
            Magmeanval = [Magmean,RealMagMean];
            Magmeanvals = vertcat(Magmeanvals,Magmeanval);
        end
        
        % Important to enable the upload of following 3D time point.
        Array1 = Array2;
        ctc=300
    end
end

%% try to save Xcorrevals

% create output folders if not already there
if ~exist(fullfile(pathname, 'Results'))
    mkdir(fullfile(pathname, 'Results'))
end

% save all data to .mat file
save(fullfile(pathname, 'Results', ...
    ['Data_', filewords,'.mat']));

% save Mean in micron/min to .mat file
save(fullfile([pathname '/Results'], ...
    ['Mag_', filewords, '.mat']), ...
    'Magmeanvals');
save(fullfile([pathname '/Results'], ...
    ['Drift_', filewords, '.mat']), ...
    'Driftcorrevals');
toc 

%         %end
%         %first_time_point = first_time_point + 1;
%         % Saving Magnitude of shifts in pxl at Index position in 3D for each time point.
%         if ~exist(fullfile(pathname, 'Shift'))
%             mkdir(fullfile(pathname, 'Shift'))
%         end
%
%         save(fullfile([pathname '/Shift'], ...
%             ['Magnm', filename, '.mat']), ...
%             'Xcorrmagnm');
%         save(fullfile([pathname '/Shift'], ...
%             ['Xcorrx', filename, '.mat']), ...
%             'Xcorrx');
%         save(fullfile([pathname '/Shift'], ...
%             ['Xcorry', filename, '.mat']), ...
%             'Xcorry');
%         save(fullfile([pathname '/Shift'], ...
%             ['Xcorrz', filename, '.mat']), ...
%             'Xcorrz');
%
%     save(fullfile([pathname '/Shift'], ...
%         ['Mag', filename, '.mat']), ...
%         'Xcorrmag');
% Include ones for xcorrx, xcorry, xcorrz.