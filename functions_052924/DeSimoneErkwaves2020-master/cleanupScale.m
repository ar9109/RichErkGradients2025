function [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts)
%     Copyright (C) 2020  Alessandro De Simone
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts)
% 
% myScale is (obviously) the myScale object
% refCh is the reference channel to use for cleanup
% chToClean chooses what channels to clean
% roi is the roi to use for cleanup (no re-calculation of the roi)

%%

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;

   basename=myScale.basename;
   opts.basename = basename; % edited 032624

   mkdir(outFolder);
   
   ch=refCh;
   
   % parsing opts field
   if(isfield(opts,'type'))
       type=opts.type;
   else
       type='manual';
   end
       
   if(isfield(opts,'segOpts'))
        segOpts=opts.segOpts;
   else
         segOpts=[];
   end
   
   if(isfield(opts,'verbose'))
       verbose=opts.verbose;
   else
       verbose=[]; 
   end
   
% whether to use BF, edited 110922
    if(isfield(opts,'BF')) 
        if ~isempty(opts.BF)
            opts.BFpath = [inFolder basename '_ch' num2str(opts.BF) 'maxproj.tif'];
%             opts.BFpath = [inFolder basename '_ch' num2str(opts.BF)];
        end
    end
% whether to create roiTight, edited 052523
    if isfield(opts,'roiTight') && ~isempty(opts.roiTight)
        for c = opts.roiTight
            opts.roiTightpath{c} = [inFolder basename '_ch' num2str(c) 'maxproj.tif'];
        end
    end
    
   % writing options
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   % read old roi
   if(isfield(opts,'readROI'))
       if(~isempty(opts.readROI))
            roipath=[paths.(opts.readROI) basename '_ROI'  '.tif'];
            if(exist(roipath,'file'))
                oldroi= loadtiff(roipath);
                oldroi=logical(oldroi);
            else
                oldroi=[];
            end
       end
   else
        oldroi=[];
   end
      
   
   if(isempty(roi))
       if numel(ch)>1 || any(ch==0) % sum all channels % edited 032624, allow multiple channels
           if any(ch==0) % edited 032624, allow multiple channels
               ch = 1:myScale.nch;
           end

           for c = ch % edited 032624, allow multiple channels

               spath=[inFolder basename '_ch' num2str(c) '.tif'];

                %using loadtiff
%                 stack_ch = loadtiff(spath); 
%                 if(c==1)
%                     stack=stack_ch;
%                 else
%                     stack=max(stack,stack_ch);
%                 end

                % new way to make it brighter edited 032724
                stack_ch3d=loadtiff(spath); 
                if(c==1)
                    stack.YX = imadjust(squeeze(max(stack_ch3d,[],3)));
                    stack.YZ = imadjust(squeeze(max(stack_ch3d,[],2)));
                    stack.ZX = imadjust(squeeze(max(stack_ch3d,[],1))');
                else
                    stack.YX = max(stack.YX, imadjust(squeeze(max(stack_ch3d,[],3))));
                    stack.YZ = max(stack.YZ, imadjust(squeeze(max(stack_ch3d,[],2))));
                    stack.ZX = max(stack.ZX, imadjust(squeeze(max(stack_ch3d,[],1))'));
                end
           end
           
           if(~isempty(oldroi))
                stack(~oldroi)=0;
           end
           
           if(strcmp(type,'manual'))
               [stack,roi,roiTight]=cleanupStack(stack,roi,opts); 
           elseif(strcmp(type,'KTR'))
               [stack,roi]=cleanupScaleSegKTR(stack,segOpts,verbose); 
           end
       else

           spath=[inFolder basename '_ch' num2str(ch) '.tif'];

           myScale.cleanedScale=struct;
           myScale.cleanedScale.folder=outFolder;

           
           %using loadtiff
           tload = tic;
           stack=loadtiff(spath);    
           disp(['Image loading lapse = ' num2str(toc(tload))]);

           
           if(~isempty(oldroi))
                stack(~oldroi)=0;
           end
           
           
           if(strcmp(type,'manual'))
               [stack,roi,roiTight]=cleanupStack(stack,oldroi,opts); 
           elseif(strcmp(type,'KTR'))
               [stack,roi]=cleanupScaleSegKTR(stack,segOpts,verbose); 
           end

           if(any(chToClean==refCh))

               %write channel
               wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
               myScale.cleanedScale.paths{ch}=wpath;

               saveastiff(stack, wpath, options);

               %write projections
               wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
               saveastiff(max(stack,[],3),wpath,options); 
           end
       end
   end

   newroi=roi;
   clear('roi'); clear('oldroi');
   %write roi
   if isstruct(newroi) % edited 120822 to allow using faster code
       tsave = tic;
       mkdir(paths.fastFolder);
       % BF
       if isfield(newroi,'BF')
           wpath = [paths.fastFolder filesep basename '_ROIprojBF'  '.tif']; 
           saveastiff(im2uint8(newroi.BF), wpath, options);
       end
       % YX
       wpath = [paths.fastFolder filesep basename '_ROIprojYX'  '.tif']; 
       saveastiff(im2uint8(newroi.YX), wpath, options);
       % YZ
       wpath = [paths.fastFolder filesep basename '_ROIprojYZ'  '.tif']; 
       saveastiff(im2uint8(newroi.YZ), wpath, options);
       % ZX
       wpath = [paths.fastFolder filesep basename '_ROIprojZX'  '.tif']; 
       saveastiff(im2uint8(newroi.ZX), wpath, options);
       disp(['Image saving lapse = ' num2str(toc(tsave))]);
   else
       wpath=[outFolder basename '_ROI'  '.tif']; 
       myScale.cleanedScale.roipath=wpath;
       saveastiff(im2uint8(newroi), wpath, options);
   end
  %write roiTight, edited 071922
    if exist('roiTight','var')
        if ~isempty(roiTight)
            for i = 1:numel(roiTight)
            wpath=[paths.fastFolder basename '_ROIprojYX_tight_c' num2str(roiTight(i).channel) '.tif']; 
            saveastiff(im2uint8(roiTight(i).mask), wpath, options);
            end
        end
    end
  
   for ch=chToClean
      
%        if(ch==refCh)
%           continue; 
%        end
       
       spath=[inFolder basename '_ch' num2str(ch) '.tif'];
       infoStack=imfinfo(spath);
       
       %using loadtiff
       stack=loadtiff(spath);
       stack(~newroi)=0;   

       %write channel
       wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
       myScale.cleanedScale.paths{ch}=wpath;
        
       saveastiff(stack, wpath, options);
          
       %write projections
       wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
       saveastiff(max(stack,[],3),wpath,options); 

   end
   
   %save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');

end

