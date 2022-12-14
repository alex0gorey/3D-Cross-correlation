function [Truex,Truey,Truez,valz] = run_PIV (maxi,maxj,maxl,Array1,Array2,grid_size_xy,grid_size_z,source_size_xy,source_size_z,search_size_xy,search_size_z)

for i = 1:maxi % Loop moving across in X Grid Distance
    for j = 1:maxj % Loop moving across in Y Grid Distance
        for l = 1:maxl % Loop moving across in Z Grid Distance
            xindex = 1+(i-1)*grid_size_xy;
            yindex = 1+(j-1)*grid_size_xy;
            zindex = 1+(l-1)*grid_size_z;
            
            % try
            % 3D source volume moving through Grid location
            sub_roi3D = Array1...
               (xindex - ceil(source_size_xy):...
                xindex + ceil(source_size_xy),...
                yindex - ceil(source_size_xy):...
                yindex + ceil(source_size_xy),...
                zindex - ceil(source_size_z):...
                zindex + ceil(source_size_z));
            % catch
            %     continue
            % end
            %
            % try
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
            
            % catch
            %     continue
            % end
            
            % Convolution Function in 3D from sub-ROI, SubArea
            % try
            cor_3d=convn(sub_roi3D,sub_area3D(end:-1:1,end:-1:1,end:-1:1));
            %                             f= f + 1
            
            % catch
            %     continue
            % end
            
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
            % Shift per grid point in nm
            
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
end
