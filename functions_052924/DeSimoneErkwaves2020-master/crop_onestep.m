function [myScale] = crop_onestep(myScale,chToCrop,paths)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
inFolder = paths.inFolder;
outFolder= paths.outFolder;
objFolder= paths.objFolder;
fieldName=paths.fieldName;

basename=myScale.basename;

mkdir(outFolder);
myScale.(fieldName).cropfolder=outFolder;
options=[];
options.overwrite='true';
options.compress='lzw';

crop_rect = myScale.(fieldName).crop_rect;

for ch=chToCrop
    
    spath=[inFolder basename '_ch' num2str(ch) '.tif'];      
    
    %using loadtiff
    stack=loadtiff(spath);
    
    [outStack] = stack(crop_rect(1,1):crop_rect(1,2),...
        crop_rect(2,1):crop_rect(2,2),...
        crop_rect(3,1):crop_rect(3,2));
    
    %write channel
    wpath=[outFolder filesep basename '_ch' num2str(ch) '.tif']; 
    
    saveastiff(outStack, wpath, options);
      
    %write projections
    wpath=[outFolder filesep basename '_ch' num2str(ch) 'maxproj.tif'];   
    saveastiff(max(outStack,[],3),wpath,options); 
end

% crop ROI(s)
st_dir = dir([inFolder basename '_ROI*' '.tif']);
if ~isempty(st_dir)

    for i=st_dir'
        
        spath=[inFolder i.name];      
        
        %using loadtiff
        stack=loadtiff(spath);
        
        [outStack] = stack(crop_rect(1,1):crop_rect(1,2),...
            crop_rect(2,1):crop_rect(2,2),...
            crop_rect(3,1):crop_rect(3,2));
        
        %write channel
        wpath=[outFolder filesep i.name]; 
        
        saveastiff(outStack, wpath, options);
    
    end
end
%save mat file
wpathmat=[objFolder basename '.mat'];  
save(wpathmat,'myScale');
end

