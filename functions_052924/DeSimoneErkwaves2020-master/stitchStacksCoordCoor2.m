function [stitched,offsets] = stitchStacksCoordCoor2(mys,opts,paths)%refCh,targetVox,cropImage,verbose)
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
%
%

vox1=mys.metadata{1}.voxel;
%targetVox=vox1;

bitDepth=mys.metadata{1}.BitDepth;

% read options
% ref channel
refCh=opts.refCh;
if numel(refCh) < size(mys.dir_paths,1) % make it a list corresponding to each position, edited 091323
    refCh = [refCh, repmat(refCh(end), [1,size(mys.dir_paths,1)-numel(refCh)])];
end
% Use metadata position
if isfield(opts,'meta_pos')
    meta_pos = opts.meta_pos;
else
    meta_pos = 0;
end
if numel(meta_pos) < size(mys.dir_paths,1) % make it a list corresponding to each position, edited 020524
    meta_pos = [meta_pos, repmat(meta_pos(end), [1,size(mys.dir_paths,1)-numel(meta_pos)])];
end
% for to go leftwards % edited 020924
if isfield(opts,'force_leftwards')
    force_leftwards = opts.force_leftwards;
else
    force_leftwards = 1;
end
if numel(force_leftwards) < size(mys.dir_paths,1) % make it a list corresponding to each position
    force_leftwards = [force_leftwards, repmat(force_leftwards(end), [1,size(mys.dir_paths,1)-numel(force_leftwards)])];
end
% others
targetVox=opts.targetVox;
cropImage=opts.cropImage;
verbose=opts.verbose;

if(isfield(opts,'resizeFactors'))
    resizeFactors=opts.resizeFactors;
else
    resizeFactors=[1 1];
end


if(isfield(opts,'overlaps'))
    overlaps=opts.overlaps;
else
    overlaps=[];
end

% build stack1
k=1;

if(size(mys.dir_paths{k,refCh(k)},2)==1)
        stack1=im2uint8(loadtiff(mys.dir_paths{k,refCh(k)}{1}));
else
    stack1=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
    for n=1:size(mys.dir_paths{k,refCh(k)},2)
        img_temp=loadtiff(mys.dir_paths{k,refCh(k)}{n});
        stack1(:,:,n)=im2uint8(img_temp);
    end
end

% we rescale the stack
stack1=rescaleStack(stack1,vox1,targetVox,verbose);

%the first image fills the entire stack
subImages=[1 1 1 size(stack1)];

% we convert positions in pixels - posZ now is the first stack
if true
% if mys.metadata{1}.FlipX
    X1_1 = mys.metadata{1}.PosY/targetVox(1);
else
    X1_1 = -mys.metadata{1}.PosY/targetVox(1);
end
if true
% if mys.metadata{1}.FlipY
    X1_2 = -mys.metadata{1}.PosX/targetVox(2);
else
    X1_2 = mys.metadata{1}.PosX/targetVox(2);
end
X1=[X1_1 X1_2 mys.metadata{1}.PosZ/targetVox(3)];

stitched{refCh(k)}=stack1;

for ch=1:size(mys.dir_paths,2)
    if(ch==refCh(k))
        continue
    else
        
        %build stack1
        if(size(mys.dir_paths{1,ch},2)==1)
             stack_ch=im2uint8(loadtiff(mys.dir_paths{1,ch}{1}));
        else
             stack_ch=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
             for n=1:size(mys.dir_paths{1,refCh(k)},2)
                img_temp=loadtiff(mys.dir_paths{1,ch}{n});
                stack_ch(:,:,n)=im2uint8(img_temp);
             end
        end
        % we rescale the stack
        stack_ch=rescaleStack(stack_ch,vox1,targetVox,verbose);
        stitched{ch}=stack_ch;
    end
end

offsets=[0 0 0];

for k=2:size(mys.dir_paths,1)

    if verbose
    display(['-----------']);
    display(['Now stitching Position ' num2str(mys.Position_list(k))]);
    end
    stack1=stitched{refCh(k)};
    
    vox2=mys.metadata{k}.voxel;
     
    % build stack2
    if(size(mys.dir_paths{k,refCh(k)},2)==1)
        stack2=im2uint8(loadtiff(mys.dir_paths{k,refCh(k)}{1}));
    else
        stack2=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
        for n=1:size(mys.dir_paths{k,refCh(k)},2)
            img_temp=loadtiff(mys.dir_paths{k,refCh(k)}{n});
            stack2(:,:,n)=im2uint8(img_temp);
        end
    end

    % we rescale the stack
    stack2=rescaleStack(stack2,vox2,targetVox,verbose);
   
    if(~mys.metadata{k}.voxel==vox1)
        error('Images do not have the same voxel size!');
    end
    if true
%     if mys.metadata{k}.FlipX
        X2_1 = mys.metadata{k}.PosY/targetVox(1);
    else
        X2_1 = -mys.metadata{k}.PosY/targetVox(1);
    end
    if true
%     if mys.metadata{k}.FlipY
        X2_2 = -mys.metadata{k}.PosX/targetVox(2);
    else
        X2_2 = mys.metadata{k}.PosX/targetVox(2);
    end
    X2=[X2_1 X2_2 mys.metadata{k}.PosZ/targetVox(3)];
    size1=size(stack1);
    size2=size(stack2);
   
    %offsets using coordinates
    %yOffset=round(X2(1)-X1(1));
    %xOffset=round(X2(2)-X1(2));
    %Offset=round(X2(3)-X1(3));
    
    % we transform in discrete coordinates
    X1=floor(X1);X2=floor(X2);
    
    if(~isempty(overlaps))
        if(~isempty(overlaps{k}))
            XStart1=overlaps{k}.XStart1;
            XEnd1=overlaps{k}.XEnd1;
            XStart2=overlaps{k}.XStart2;
            XEnd2=overlaps{k}.XEnd2;
            % TO DO: modify so that you provide XStart2 and XEnd2 and the
            % subImage to use
        else
            % we calulate the overlap between stack2 and a subImage in stack1
            [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage,verbose);
            if verbose
            disp(['Subimage used: Position' num2str(mys.Position_list(usedSub))]);
            end
        end
    else
            if opts.allSubs==0 % & isempty(mys.metadata{k}.zBest) % allow zBest, edited 081622
                % we calulate the overlap between stack2 and a subImage in stack1
                [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage,verbose);
                if verbose
                disp(['Subimage used: Position' num2str(mys.Position_list(usedSub))]);
                end

                if(isempty(XStart1))
                    if verbose
                    disp('There is no overlap, so I check the overlap with all subs');
                    end
                end
            else
                XStart1=[];XEnd1=[];
                XStart2=[];XEnd2=[];
            end
            
            if(~isempty(mys.metadata{k}.zBest))
                if verbose
                   disp('You have provided the best z'); %, so I check the overlap with all subs');
                end
            end
        
    end
    
%     if false
    if (~isempty(XStart1))&&meta_pos(k) && all(abs(XEnd2(1:2)-XStart2(1:2))>=30) % edited 122022, correct for too small overlap

        yStart1=XStart1(1);    xStart1=XStart1(2);    
%         zStart1=XStart1(3);
        yStart2=XStart2(1);    xStart2=XStart2(2);    
%         zStart2=XStart2(3);
        zStart1 = 1;
        zStart2 = 1;

        yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
        yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);

        % dont crop in z, edited 080322
        zEnd1 = size1(3);
        zEnd2 = size2(3);

        %we carve out the templates
%         template1=stack1(yStart1:yEnd1,xStart1:xEnd1,zStart1:zEnd1);
%         template2=stack2(yStart2:yEnd2,xStart2:xEnd2,zStart2:zEnd2);
        template1=stack1(yStart1:yEnd1,xStart1:xEnd1,:); % prevent wrong z-pos remove too much from overlap, edited 080122
        template2=stack2(yStart2:yEnd2,xStart2:xEnd2,:); % prevent wrong z-pos remove too much from overlap, edited 080122

        %first round of alignment using a coarse-grained image

        template1small=imresize(template1,resizeFactors(1));
        template2small=imresize(template2,resizeFactors(1));
        [relOffset,corrValue,zbest2]=correlateStacks(template1small,template2small,mys.metadata{k}.zBest,verbose);
        relOffset=[relOffset(1)/resizeFactors(1) relOffset(2)/resizeFactors(1) relOffset(3)];
        
        % second round using only zbest1 and zbest2 and using the finer
        % movements
        zbest1=zbest2+relOffset(3);
        corrhere=normxcorr2(template2(:,:,zbest2),template1(:,:,zbest1));
        % use maxproj for xy stitching
%         corrhere=normxcorr2(max(template2,[],3),max(template1,[],3));

        [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
        relOffset(1) = ypeak-size(template2,1);
        relOffset(2) = xpeak-size(template2,2);
        corrValue=[corrhere(ypeak,xpeak)];

        % debug
        if k==2
            1;
        end

    else
        
       corrValue=-Inf;
       
       % subimages to loop through, forcing always use the last position
        if force_leftwards(k)
            sub_list = size(subImages,1);
        else
            sub_list = 1:size(subImages,1);
        end

       for sub=sub_list 
          
            subHereRel=subImages(sub,:);
            
            XStart1= subHereRel(1:3);
            XEnd1  = subHereRel(4:6);
            
            XStart2= [1 1 1];
%             XStart2= [1 size2(2)-100 1]; % Only use the right most 100 px region, edited 072922
            if force_leftwards(k)
                XStart2= [1 size2(2)-200 1]; % fix stitching, edited 031023
            end
            XEnd2  = size2;

%             XEnd2  = size2-[1 100 1];
             
            yStart1=XStart1(1);    xStart1=XStart1(2);    zStart1=XStart1(3);
            yStart2=XStart2(1);    xStart2=XStart2(2);    zStart2=XStart2(3);

            yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
            yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);

            %we carve out the templates
            template1=stack1(yStart1:yEnd1,xStart1:xEnd1,zStart1:zEnd1);
            template2=stack2(yStart2:yEnd2,xStart2:xEnd2,zStart2:zEnd2);
            
            template1small=imresize(template1,resizeFactors(2));
            template2small=imresize(template2,resizeFactors(2));
            [relOffsetHere,corrValueHere,zbest2]=correlateStacks(template1small,template2small,mys.metadata{k}.zBest,verbose);

            relOffsetHere=[relOffsetHere(1)/resizeFactors(2) relOffsetHere(2)/resizeFactors(2) relOffsetHere(3)]; 
            
            zbest1=zbest2+relOffsetHere(3);
            corrhere=normxcorr2(template2(:,:,zbest2),template1(:,:,zbest1));
            % use maxproj for xy stitching
%             corrhere=normxcorr2(max(template2,[],3),max(template1,[],3));

            if isfield(opts,'laplacian_corr')&&opts.laplacian_corr
%                 corrhere = -del2(corrhere); % using laplacian of the corr matrix instead; edited 031423
                H = fspecial('log',10,3);
                corrhere = abs(imfilter(corrhere,H,'replicate')); % using log instead; edited 020624
            end
            
            if force_leftwards(k)
                [ypeak, xpeak] = find(corrhere==max(corrhere(200:end-200,1:800),[],'all'),1,'first');
%                 [ypeak, xpeak] = find(corrhere==max(corrhere(:),[],'all'),1,'first');
            else
                [ypeak, xpeak] = find(corrhere==max(corrhere(:),[],'all'),1,'first');
            end

            relOffsetHere(1) = ypeak-size(template2,1);
            relOffsetHere(2) = xpeak-size(template2,2);
            corrValueHere=[corrhere(ypeak,xpeak)];
            
            % debug
            if k==2
                1;
            end


            if verbose
            sub
            relOffsetHere
            corrValueHere
            end
            
            if(corrValueHere > corrValue)
                chosenSub = sub;
                relOffset = relOffsetHere;
                corrValue =corrValueHere;
            end           
       end 
       
       if verbose
       disp(['Subimage used: Position' num2str(mys.Position_list(chosenSub))]);
       end
       subHereRel=subImages(chosenSub,:);        
       XStart1= subHereRel(1:3); XEnd1  = subHereRel(4:6);     
       XStart2= [1 1 1];
%        XStart2= [1 size2(2)-100 1]; % Only use the right most 100 px region, edited 072922
       if force_leftwards(k)
           XStart2= [1 size2(2)-200 1]; % fix stitching, edited 031023
       end
        XEnd2  = size2;

%             XEnd2  = size2-[1 100 1];

       yStart1=XStart1(1);    xStart1=XStart1(2);    zStart1=XStart1(3);
       yStart2=XStart2(1);    xStart2=XStart2(2);    zStart2=XStart2(3);

       yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
       yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);
       
    end
    clear('template1');clear('template2');
    
    if verbose
    disp('Stitching ...');
    end

    % we calculate the offset with respect to the total stack1
    yOffset=relOffset(1)+yStart1-1-(yStart2-1);
    xOffset=relOffset(2)+xStart1-1-(xStart2-1);
    zOffset=relOffset(3)+zStart1-1-(zStart2-1);
    
    overwrite1 = nan; % edited 032624 overwrite1
    [newImage,X1new,subImagesNew,overwrite1]=stitch(stitched{refCh(k)},stack2,[yOffset,xOffset,zOffset],X1,subImages,overwrite1); % edited 032624 overwrite1
    stitched{refCh(k)}=newImage;
    
    % time to stitch images (finally)  
    for ch=1:size(mys.dir_paths,2)
        if(ch==refCh(k))
            continue
        end
        
        if(size(mys.dir_paths{k,refCh(k)})==1)
            stack2stitch=im2uint8(loadtiff(mys.dir_paths{k,ch}{1}));
        else
            stack2stitch=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
            for n=1:size(mys.dir_paths{k,refCh(k)},2)
                img_temp=loadtiff(mys.dir_paths{k,ch}{n});
                stack2stitch(:,:,n)=im2uint8(img_temp);
            end
        end
        
        stack2stitch=rescaleStack(stack2stitch,vox2,targetVox,verbose);
        [newImage,X1new,subImagesNew]=stitch(stitched{ch},stack2stitch,[yOffset,xOffset,zOffset],X1,subImages,overwrite1); % edited 032624 overwrite1
        stitched{ch}=newImage;
    end

    % update coordinates and subImages
    offsets(k,:)=[yOffset,xOffset,zOffset];
    X1=X1new;
    subImages= subImagesNew;

    %with saveastiff save projections, edited 083022
    proj_xy = max(stitched{refCh(k)},[],3);
    proj_yz = squeeze(max(stitched{refCh(k)},[],2)); 
    options=[];
    options.overwrite='true';
    options.compress='lzw';
    wpath_xy = [paths.rayFolder filesep 'position_' num2str(mys.Position_list(k)) '_projxy.tif'];
    saveastiff(proj_xy,wpath_xy,options);
    wpath_yz = [paths.rayFolder filesep 'position_' num2str(mys.Position_list(k)) '_projyz.tif'];
    saveastiff(proj_yz,wpath_yz,options);
    %
   %
    if(verbose==1)
      h1=figure;
      imshow(proj_xy,[]);
      h2=figure;
      imshow(proj_yz,[]);
      %waitforbuttonpress;
      %close(h1);close(h2);
    end
end
    
if(verbose==2)
      h1=figure;
      imshow(max(stitched{refCh(k)},[],3),[]);
      h2=figure;
      imshow(squeeze(max(stitched{refCh(k)},[],2)),[]);
      %waitforbuttonpress;
      %close(h1);close(h2);
end

end

function [newStack12,newCoord,subImages,overwrite1]=stitch(stack1,stack2,offsets,coord,subImages,overwrite1) % edited 032624 overwrite1

    yOffset=offsets(1);xOffset=offsets(2);zOffset=offsets(3);
     
    size1=size(stack1);
    size2=size(stack2);
    
    yStart=min(1,1+yOffset);
    xStart=min(1,1+xOffset);
    zStart=min(1,1+zOffset);
    
    yEnd=max(size1(1),size2(1)+yOffset);
    xEnd=max(size1(2),size2(2)+xOffset); 
    zEnd=max(size1(3),size2(3)+zOffset);
    
    newCoord=coord+[yStart-1,xStart-1,zStart-1];
    sizeNew=[yEnd-yStart+1,xEnd-xStart+1,zEnd-zStart+1];   
    
    % subImages - [yStartSubIm1 xStartSubIm1 zStartSubIm1 yEndSubIm1 xEndSubIm1 zEndSubIm1; ...] 
    % we move subImages of stack1
    subImages(:,1)=subImages(:,1)+1-yStart;subImages(:,4)=subImages(:,4)+1-yStart;
    subImages(:,2)=subImages(:,2)+1-xStart;subImages(:,5)=subImages(:,5)+1-xStart;
    subImages(:,3)=subImages(:,3)+1-zStart;subImages(:,6)=subImages(:,6)+1-zStart;
    
    % overlap , edited 080122
    y_overlap1 = [max(1,1+yOffset),min(size1(1),size2(1)+yOffset)];
    x_overlap1 = [max(1,1+xOffset),min(size1(2),size2(2)+xOffset)];
    z_overlap1 = [max(1,1+zOffset),min(size1(3),size2(3)+zOffset)];

    overlap1 = stack1(y_overlap1(1):y_overlap1(2),...
        x_overlap1(1):x_overlap1(2),...
        z_overlap1(1):z_overlap1(2));

    overlap2 = stack2(max(1,1-yOffset):min(size1(1)-yOffset,size2(1)),...
        max(1,1-xOffset):min(size1(2)-xOffset,size2(2)),...
        max(1,1-zOffset):min(size1(3)-zOffset,size2(3)));

    if isnan(overwrite1) % edited 032624 overwrite1
        overlap1_mean = double(mean(overlap1,[1,3]));
        overlap2_mean = double(mean(overlap2,[1,3]));
        overlap_diff = movmean(overlap2_mean-overlap1_mean,0.5.*numel(overlap1_mean));
    %     min_idx = find(overlap_diff==min(overlap_diff));
    %     split_pos = min_idx(end);
        overlap_neg = find(overlap_diff<0);
        if isempty(overlap_neg)
            overwrite1 = [];
        elseif mean(overlap_neg)>=round(numel(overlap_diff)/2)
            overwrite1 = overlap_neg(1):numel(overlap_diff);
        else
            overwrite1 = 1:overlap_neg(end);
        end
    end


%     %we create a logical image to know where the subimages are
%     isNew=ones(sizeNew,'logical');
%     
%     for sub=1:size(subImages,1)
%         isNew(subImages(sub,1):subImages(sub,4),subImages(sub,2):subImages(sub,5),subImages(sub,3):subImages(sub,6))=0;
%     end

%     newStack=zeros(sizeNew,'uint8');
    newStack12=zeros(sizeNew,'uint8');
    
%     %we position stack2
%     newStack(1+yOffset-yStart+1:size2(1)+yOffset-yStart+1,1+xOffset-xStart+1:size2(2)+xOffset-xStart+1,1+zOffset-zStart+1:size2(3)+zOffset-zStart+1)=stack2;
%     % we position stack1 - it will make some black regions
%     newStack(1-yStart+1:size1(1)-yStart+1,1-xStart+1:size1(2)-xStart+1,1-zStart+1:size1(3)-zStart+1)=stack1;
   
    % we position stack1 - it has the problem of bleaching  
    newStack12(1-yStart+1:size1(1)-yStart+1,1-xStart+1:size1(2)-xStart+1,1-zStart+1:size1(3)-zStart+1)=stack1;
    %we position stack2
    newStack12(1+yOffset-yStart+1:size2(1)+yOffset-yStart+1,1+xOffset-xStart+1:size2(2)+xOffset-xStart+1,1+zOffset-zStart+1:size2(3)+zOffset-zStart+1)=stack2;
    
    % add back the overlap rigion (by split_pos) from stack1, edited 080122
    new_overwrite1 = 1-xStart+x_overlap1(1)+overwrite1-1;

    newStack12(y_overlap1(1)+1-yStart:y_overlap1(2)+1-yStart,...
            new_overwrite1,...
            z_overlap1(1)+1-zStart:z_overlap1(2)+1-zStart)...
            = overlap1(:,overwrite1,:);

%     new_split_pos = 1-xStart+x_overlap1(1)+split_pos-1;
%     if xOffset<0
%         newStack12(y_overlap1(1)+1-yStart:y_overlap1(2)+1-yStart,...
%             new_split_pos:x_overlap1(2)+1-xStart,...
%             z_overlap1(1)+1-zStart:z_overlap1(2)+1-zStart)...
%             = overlap1(:,split_pos:end,:);
%     else
%         newStack12(y_overlap1(1)+1-yStart:y_overlap1(2)+1-yStart,...
%             x_overlap1(1)+1-xStart:new_split_pos,...
%             z_overlap1(1)+1-zStart:z_overlap1(2)+1-zStart)...
%             = overlap1(:,1:split_pos,:);
%     end
   

    % we use newStack, but substitute the new regions with newStack12
%     newStack(isNew)=newStack12(isNew);
%     clear('newStack12'); clear('isNew');

    %we attach the subImage of stack2
    newSubImage=[1+yOffset-yStart+1, 1+xOffset-xStart+1,1+zOffset-zStart+1];
    newSubImage=[newSubImage, size2(1)+yOffset-yStart+1,size2(2)+xOffset-xStart+1, size2(3)+zOffset-zStart+1];
    subImages=vertcat(subImages,newSubImage);
    
end


function [offset,corrValue,zbest2]=correlateStacks(stack1,stack2,zbest2,verbose)
   
    if(isempty(zbest2))
        
%         % we need to choose which plane to use to compare with stack1
%         % we believe that stack2 is overlapping with stack 1
% 
        % aligns the projections
%         s1proj=max(stack1,[],3);
%         s2proj=max(stack2,[],3);
% 
%         corrhere=normxcorr2(s2proj,s1proj);
%         [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
%         xOffset = xpeak-size(s2proj,2);
%         xOverStart= max(1-xOffset,1); % this is the ref system of stack2
%         xOverEnd  = min(size(s1proj,2)-xOffset,size(s2proj,2));
% 
%         xbest = round((xOverStart+xOverEnd)./2);
% 
%         % this method aligns in z
%         s1_xz=squeeze(stack1(:,xbest+xOffset,:));
%         s2_xz=squeeze(stack2(:,xbest,:));    
% 
%         if(size(s1_xz,2)>size(s2_xz,2))
%             corrhere=normxcorr2(s2_xz,s1_xz);
%             [ypeak, zpeak] = find(corrhere==max(corrhere(:)),1,'first');
%             zOffset = zpeak-size(s2_xz,2);
%         else
%             corrhere=normxcorr2(s1_xz,s2_xz);
%             [ypeak, zpeak] = find(corrhere==max(corrhere(:)),1,'first');
%             zOffset = -zpeak+size(s1_xz,2);
%         end
% 
%         zOverStart= max(1-zOffset,1); % this is the ref system of stack2
%         zOverEnd  = min(size(s1_xz,2)-zOffset,size(s2_xz,2));
%      
%         % choose the brightest plane in the overlap
%         [m,zbest2]=max(squeeze(mean(mean(stack2(:,:,zOverStart:zOverEnd),1),2)));
%         zbest2 =zbest2+zOverStart-1;

          [m,zbest2]=max(squeeze(mean(stack2(:,:,:),[1,2])));
%           zbest2 = floor(size(stack2,3)./2);

    end
     
    best2=stack2(:,:,zbest2);
    
    if(std(im2double(best2(:)))==0)
       error('The selected portion of stack2 is completely flat!');
    end
    
   for z=1:size(stack1,3)
        if(std(im2double(stack1(:,:,z)))==0)
            if verbose
            display(['We skip z= ' num2str(z) 'because it is flat']);
            end
            continue;
        end
        
        corrhere=normxcorr2(best2,stack1(:,:,z));
        [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
        yOffset = ypeak-size(best2,1);
        xOffset = xpeak-size(best2,2);
        corrValueZ(z,:)=[yOffset xOffset corrhere(ypeak,xpeak)];
        
    end
    
    [corrValue,zpeak]=max(corrValueZ(:,3));
    zOffset=zpeak-1-zbest2+1;
    yOffset=corrValueZ(zpeak,1); xOffset=corrValueZ(zpeak,2);
    
    offset=[yOffset,xOffset,zOffset];

end

function [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage,verbose)
% X1, X2 are the absolute coordinates of the two stacks
% subImages contains the relative coordinates of the subImages

    XStart1=[];
    XEnd1=[];

    XStart2=[];
    XEnd2=[];

    usedSub=[];

    for sub=size(subImages,1):-1:1

        subHereRel=subImages(sub,:);
       
        %we calculate absolute coordinates subImage   
        subHereAbs(1:3)=subHereRel(1:3)-1+X1(1:3);
        subHereAbs(4:6)=subHereRel(4:6)-1+X1(1:3);

        % lets calculate overlap with subImages in absolute coordinates
        yStart=max(subHereAbs(1),X2(1));
        xStart=max(subHereAbs(2),X2(2));
        zStart=max(subHereAbs(3),X2(3));
    
        yEnd=min(subHereAbs(4),X2(1)+size2(1)-1);
        xEnd=min(subHereAbs(5),X2(2)+size2(2)-1);
        zEnd=min(subHereAbs(6),X2(3)+size2(3)-1);    
        
        if ((yEnd-yStart+1)>0)&&((xEnd-xStart+1)>0) %&((zEnd-zStart+1)>0)) % dont crop in z, edited 080322
            
            if(isempty(XStart1)||(yEnd-yStart+1)*(xEnd-xStart+1)>(XEnd2(1)-XStart2(1)+1)*(XEnd2(2)-XStart2(2)+1)) % corr with subimages with the most overlap
  
                    %we convert in relative coordinates in stack1 and stack2
                    yStart1=yStart-X1(1)+1;xStart1=xStart-X1(2)+1;zStart1=zStart-X1(3)+1;
                    yStart2=yStart-X2(1)+1;xStart2=xStart-X2(2)+1;zStart2=zStart-X2(3)+1;

                    yEnd1=yEnd-X1(1)+1;xEnd1=xEnd-X1(2)+1;zEnd1=zEnd-X1(3)+1;
                    yEnd2=yEnd-X2(1)+1;xEnd2=xEnd-X2(2)+1;zEnd2=zEnd-X2(3)+1; 
                    
                    XStart1=[yStart1,xStart1,zStart1];
                    XEnd1=[yEnd1,xEnd1,zEnd1];

                    XStart2=[yStart2,xStart2,zStart2];
                    XEnd2=[yEnd2,xEnd2,zEnd2];
                    
                    if(cropImage==0)
                        XStart1= [subHereRel(1:3)];
                        XEnd1  = [subHereRel(4:6)];
                    end
                    
                    usedSub=sub;
            end   
        end
    end
    if(isempty(XStart1))
        if verbose
       display('This position has no overlap with the other position! '); 
        end
    end

end


function [newStack]=rescaleStack(stack,voxSize,targetVoxSize,verbose)
   
 rescale_factor=voxSize./targetVoxSize;
  
 if(any(abs((rescale_factor-[1 1 1]))>0.05))
     if verbose
    display('Now resampling image');
    display(['Rescale factor:'  num2str(rescale_factor)])
    display(['Size before:'  num2str(size(stack))])
     end
    tform=affine3d(diag([rescale_factor 1]));
    newStack=imwarp(stack,tform,'nearest','SmoothEdges',true);
    if verbose
    display(['Size after:'  num2str(size(newStack))]);
    end
%     newStack=newStack(2:end-1,2:end-1,2:end-1);
    newStack=newStack; % edited 012524, size doesnt fit corr2norm when becomes 1022 instead of 1024
 else
    newStack=stack; 
 end

end





