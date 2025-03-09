function [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths)
%  Copyright (C) 2020  Alessandro De Simone
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%% It rotates each channel in the scale  

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;
   fieldName=paths.fieldName;

   basename=myScale.basename;

   ch=refCh;
   spath=[inFolder basename '_ch' num2str(refCh) '.tif'];
  
   mkdir(outFolder);
   
   myScale.(fieldName).folder=outFolder;
   
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   %using loadtiff
   stack=loadtiff(spath);
          
   % find angles
   [outStack,ths,rects,zPos] = rotateStack(stack,ths,rects,zPos); 
   myScale.(fieldName).rotAngles=ths;
   myScale.(fieldName).cropRects=rects;
   myScale.(fieldName).zPos=zPos;
   
   if(isfield(paths,'saveName'))
        saveName=paths.saveName;
   else       
       saveName = basename;
   end
   
   if(isfield(paths,'objSaveName'))
       objSaveName=paths.objSaveName;
   else       
       objSaveName = basename;
   end
   
   %% use ROI first to find crop_rect edited 072222
    roipath = dir([inFolder basename '_ROI' '.tif']);
    if ~isempty(roipath)
        roi_prerotate=loadtiff([roipath.folder filesep roipath.name]);
        [roi] = rotateStack(roi_prerotate,ths,rects,zPos);
        
        roi_ind = find(logical(roi));
        [roiIdx1,roiIdx2,roiIdx3] = ind2sub(size(roi),roi_ind);
        crop_rect = [max(min(roiIdx1)-30,1),min(max(roiIdx1)+30,size(roi,1));...
            max(min(roiIdx2)-30,1),min(max(roiIdx2)+30,size(roi,2));...
            max(min(roiIdx3)-2,1),min(max(roiIdx3)+2,size(roi,3))];
        
        [roi_outStack] = roi(crop_rect(1,1):crop_rect(1,2),...
            crop_rect(2,1):crop_rect(2,2),...
            crop_rect(3,1):crop_rect(3,2));
        
        %write channel
        wpath=[outFolder filesep roipath.name]; 
        saveastiff(roi_outStack, wpath, options);

        myScale.(fieldName).crop_rect = crop_rect;

    end
   %%
   if(any(chToRotate==refCh))
       % crop and save
        [outStack_crop] = outStack(crop_rect(1,1):crop_rect(1,2),...
            crop_rect(2,1):crop_rect(2,2),...
            crop_rect(3,1):crop_rect(3,2));
        clear("outStack");
        %write channel
        wpath=[outFolder saveName '_ch' num2str(ch) '.tif']; 
        myScale.(fieldName).paths{ch}=wpath;

        saveastiff(outStack_crop, wpath, options);
          
        %write projections
        wpath=[outFolder saveName '_ch' num2str(ch) 'maxproj.tif'];
        saveastiff(max(outStack_crop,[],3),wpath,options); 
        clear("outStack_crop");
   
   end
   
   for ch=chToRotate
      
       if(ch==refCh)
          continue; 
       end
       spath=[inFolder basename '_ch' num2str(ch) '.tif'];      
       infoStack=imfinfo(spath);
       
       %using loadtiff
       stack=loadtiff(spath);
   
       [outStack] = rotateStack(stack,ths,rects,zPos);
        clear("stack");
       % crop and save
        [outStack_crop] = outStack(crop_rect(1,1):crop_rect(1,2),...
            crop_rect(2,1):crop_rect(2,2),...
            crop_rect(3,1):crop_rect(3,2));
        clear("outStack");
       %write channel
       wpath=[outFolder saveName '_ch' num2str(ch) '.tif']; 
       myScale.(fieldName).paths{ch}=wpath;
        
       saveastiff(outStack_crop, wpath, options);
          
       %write projections
       wpath=[outFolder saveName '_ch' num2str(ch) 'maxproj.tif'];   
       saveastiff(max(outStack_crop,[],3),wpath,options); 
       clear("outStack_crop")

   end

    % crop ROI(s)
    st_dir = dir([inFolder basename '_ROI*' '.tif']);
    if ~isempty(st_dir)
        for i=st_dir' % edited 052723 to allow ROItight
            
            spath=[inFolder i.name];      
            
            %using loadtiff
            stack=loadtiff(spath);
            
            [outStack] = rotateStack(stack,ths,rects,zPos);
            clear("stack");
           % crop and save
            [outStack_crop] = outStack(crop_rect(1,1):crop_rect(1,2),...
                crop_rect(2,1):crop_rect(2,2),...
                crop_rect(3,1):crop_rect(3,2));
            clear("outStack");

            %write channel
            wpath=[outFolder filesep i.name]; 
            
            saveastiff(outStack_crop, wpath, options);
        
        end
    end

   %save mat file
   wpathmat=[objFolder objSaveName '.mat'];  
   save(wpathmat,'myScale');

end

