function [ccrot,ampPointrot,endPointrot] = process_coords(cc,ampPoint,endPoint,pxsize)
%ROTATE_CENTER Summary of this function goes here
%   inputs ampPoint and endPoint are [y,x] formats from image; 
%   outputs are 3d vectors [x,y,z] with z=0 for ampPointrot and endPointrot
    ampPoint_z = [ampPoint(2) ...
                ampPoint(1) 0];

    endPoint_z = [endPoint(2) ...
                endPoint(1) 0];

    
    %Prolif = FUCCI>=2;
    
    % this is the part that calculates the angle
    % we fit a line across the ray
    f1 = fit(cc(:,1),cc(:,2),'poly1');
    % calculate the slope of this line as a rotation angle
    th = atan(f1.p1);
    
    % we rotate the points, see documentation of rotz and Wikipedia
    % "Rotation matrix"
    R = rotz(-th*180/pi);
    ampPointrot = (R*ampPoint_z')'.*pxsize;
    endPointrot = (R*endPoint_z')'.*pxsize;%-ampPointrot;
    % Rotate coordinates and set zero at endPoint and set leftwards
    % positive and multiply by pixel size to get length in microns
    ccrot=(-(R*(cc'.*pxsize))'+repmat(endPointrot,[size(cc,1) 1]));
end

