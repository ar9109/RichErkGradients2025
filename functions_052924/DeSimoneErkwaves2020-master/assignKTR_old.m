function [myScale,obj] = assignKTR(myScale,obj,svList,paths,opts)
%ANALYZESCALE Summary of this function goes here
%   Detailed explanation goes here
   
   inFolder=paths.inFolder;
   basename=myScale.basename;
   
   infoName=opts.infoName;
   
   KTRsuffix=  paths.KTRsuffix;
   H2Asuffix=  paths.H2Asuffix;
   
   sktr= loadtiff([inFolder basename KTRsuffix '.tif']);
   sh2a= loadtiff([inFolder basename H2Asuffix '.tif']); 
   sktr=permute(sktr,[2 1 3]); sh2a=permute(sh2a,[2 1 3]);
   
    sizes=size(sktr);
   
  % run through objects and assigns colors
  verbose = opts.verbose;
  
  colorLevels=100;
  
  if(verbose)
      scolor_g=zeros(size(sktr));
      scolor_r=zeros(size(sktr));
      scolor_b=zeros(size(sktr));
      colors=parula(2*colorLevels);
      limLogAct=2;
  end

  cent=round(vertcat(obj.m));
  
       
  if(isfield(opts,'recalcNucleus') & opts.recalcNucleus)   
    
    unflattenFolder   =  paths.unflattenFolder;
    zendsFolder       =  paths.zendsFolder;
    zendssuffix       =  paths.zendssuffix;  
      
    zends=loadtiff([zendsFolder basename zendssuffix '.tif']);
    
    cUnflat= [cent(:,1) cent(:,2) 1-cent(:,3)+double(zends(sub2ind(size(zends),cent(:,2) ,cent(:,1))))];
    %cUnflat= [cent(:,1) cent(:,2) cent(:,3)];
    
    sktrUnflat = loadtiff([unflattenFolder basename paths.KTRsuffixUnflatten '.tif']);
    sh2aUnflat = loadtiff([unflattenFolder basename paths.H2AsuffixUnflatten '.tif']); 
    
    sktrUnflat=permute(sktrUnflat,[2 1 3]); 
    sh2aUnflat=permute(sh2aUnflat,[2 1 3]);
    %zends=permute(zends,[2 1]);
    
    sizesUnflat=size(sktrUnflat);    
    
  end
  
  if(isfield(obj,'area'))
    obj =rmfield(obj,'area');
  end
  
  if(isfield(obj,'volume'))
    obj =rmfield(obj,'volume');
  end
  
  if(isfield(obj,'volume'))
    obj =rmfield(obj,'ktr');
  end

  
  
  for c=1:numel(obj)    
     
     svIdx=obj(c).svIdx;
     cell_indices=[];
            
     for sv=1:numel(svIdx)
        cell_indices=vertcat(cell_indices,svList{svIdx(sv)+1}+1);
     end
         
     
     %% if needed, recalculate nucleus
     
     if(isfield(opts,'recalcNucleus') & opts.recalcNucleus) %ACTUALLY BROKEN 
        
         error('Recalc nucleus is not working at the time. Fix it or put recalcNucleus=0');
         
         smaskUnflat=zeros(size(sh2aUnflat),'logical');
         
         [ro,cl,pp]=ind2sub(sizes,cell_indices);
         cell_indices_sub = [ro cl];
         
         %cell_indices_sub(:,3)=pp;
         cell_indices_sub(:,3) = 1-pp+double(zends(sub2ind(size(zends),cl,ro)));
         
         cell_indicesUnflat= sub2ind(sizesUnflat,cell_indices_sub(:,1),cell_indices_sub(:,2),cell_indices_sub(:,3));
         smaskUnflat(cell_indicesUnflat)=1;

         dx = round(15); dy = round(15); dz = round(15);
         thisCentUnflat=cUnflat(c,:);
     
         ymin=max(thisCentUnflat(1)-dy,1);
         ymax=min(thisCentUnflat(1)+dy,sizesUnflat(1));
         xmin=max(thisCentUnflat(2)-dx,1);
         xmax=min(thisCentUnflat(2)+dx,sizesUnflat(2));
         zmin=max(thisCentUnflat(3)-dz,1);
         zmax=min(thisCentUnflat(3)+dz,sizesUnflat(3));
     
         mini_sh2aUnflat=sh2aUnflat(ymin:ymax,xmin:xmax,zmin:zmax);
         mini_sktrUnflat=sktrUnflat(ymin:ymax,xmin:xmax,zmin:zmax);
         mini_maskUnflat=smaskUnflat(ymin:ymax,xmin:xmax,zmin:zmax);
         
         % calculate cyt and nuc masks
         mm_nuc=imclose(mini_maskUnflat,strel('disk',3));
         %mm_cyt=logical(imdilate(minimask,strel('disk',3))-mm_nuc);
     
         % calculate signals
         %volNuc=sum(mm_nuc(:));
         %volCyt=sum(mm_cyt(:));
         %signalCyt=double(sum(sum(sum(mini_sktr(mm_cyt),1),2),3))./volCyt;
         %signalNuc=double(sum(sum(sum(mini_sktr(mm_nuc),1),2),3))./volNuc;
         %ktr=signalCyt./signalNuc;         
         
         % visualize masks
         if(verbose==2)   
             
             mini_sh2aUnflatShow=  permute(mini_sh2aUnflat,   [3 2 1]);
             mm_nucShow          = permute(mm_nuc          , [3 2 1]);

             imshow(max(mini_sh2aUnflatShow,[],3),[],'InitialMagnification',1000);

             sizexy=[size(mm_nucShow,1) size(mm_nucShow,2)];
             red =   cat(3, ones(sizexy),zeros(sizexy),zeros(sizexy));
             green = cat(3, zeros(sizexy),ones(sizexy),zeros(sizexy));     
             hold on 

             h = imshow(green); 
             set(h,'AlphaData',double(max(mm_nucShow,[],3)*0.1))

             %h2 = imshow(red); 
             %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.2)
             hold off
%             ktr
             waitforbuttonpress
         end

     else
         
     
         %% calculate cytoplasmic and nuclear signal ERK KTR

         % cut sub-windows
         smask=zeros(size(sh2a),'logical');
         smask(cell_indices)=1;
         dx = round(15); dy = round(15); dz = round(15);
         cent=round(obj(c).m);     

         ymin=max(cent(1)-dy,1);
         ymax=min(cent(1)+dy,size(smask,1));
         xmin=max(cent(2)-dx,1);
         xmax=min(cent(2)+dx,size(smask,2));
         zmin=max(cent(3)-dz,1);
         zmax=max(min(cent(3)+dz,size(smask,3)),1);

         minimask=  smask(ymin:ymax,xmin:xmax,zmin:zmax);
         mini_sh2a= double(sh2a(ymin:ymax,xmin:xmax,zmin:zmax));
         mini_sktr= double(sktr(ymin:ymax,xmin:xmax,zmin:zmax));
         
         area = squeeze(sum(sum(max(minimask,[],3),1),2));
         if(isempty(area))
             area=0;
         end

         [~,zbest]=max(squeeze(mean(mean(double(minimask),1),2)));
         minimask =  minimask(:,:,zbest);
         mini_sh2a=  mini_sh2a(:,:,zbest);
         mini_sktr=  mini_sktr(:,:,zbest);

         % calculate cyt and nuc masks
         mm_nuc=imclose(minimask,strel('disk',3));
         %mm_nuc = imerode(mm_nuc,strel('disk',3)); %ALE 30jan20
         mm_cyt=logical(imdilate(minimask,strel('disk',3))-mm_nuc);
     
         % calculate signals
         volNuc=sum(double(mm_nuc(:)));
         volCyt=sum(double(mm_cyt(:)));
         signalCyt=sum(mini_sktr(mm_cyt(:)))./volCyt;
         signalNuc=sum(mini_sktr(mm_nuc(:)))./volNuc; 
         averageH2ANuc=sum(mini_sh2a(mm_nuc(:)))./volNuc;
         averageH2ACyt=sum(mini_sh2a(mm_cyt(:)))./volCyt;
         ratioH2ANucCyt=averageH2ANuc./averageH2ACyt;
         ktr=signalCyt./signalNuc;
       
        % visualize masks
        if(verbose==2)   

            sizexy=[size(minimask,1) size(minimask,2)];
            red =   cat(3, ones(sizexy),zeros(sizexy),zeros(sizexy));
            green = cat(3, zeros(sizexy),ones(sizexy),zeros(sizexy));     
            
            figure;
            % cyto mask on nuclei
            imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
            hold on 
            %h = imshow(green); 
            %set(h,'AlphaData',double(max(mm_nuc,[],3))*0.1)
            h2 = imshow(red); 
            set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.5)
            hold off
            waitforbuttonpress
            close all

            
            figure;
            %  mask on nuclei 
            imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
            hold on 
            h = imshow(red); 
            set(h,'AlphaData',bwperim(imdilate(double(max(mm_nuc,[],3)),strel('disk',1)))*0.8);
            %h2 = imshow(red); 
            %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.1)
            hold off
            waitforbuttonpress
            close all
            
            figure;
            %  nuclear mask on sensor 
            imshow(max(mini_sktr,[],3),[],'InitialMagnification',1000);
            hold on 
            h = imshow(red); 
            set(h,'AlphaData',double(max(mm_nuc,[],3))*0.2)
            %h2 = imshow(red); 
            %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.1)
            hold off
            waitforbuttonpress
            close all
            
            figure;
            % cyto mask on sensor
            imshow(max(mini_sktr,[],3),[],'InitialMagnification',1000);
            hold on 
            %h = imshow(green); 
            %set(h,'AlphaData',double(max(mm_nuc,[],3))*0.1)
            h2 = imshow(green); 
            set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.2)
            hold off
            waitforbuttonpress
            close all

            ktr
        end
     
     end
         
%%     SAVE STUFF INTO OBJ!!!
     
%      obj(c).red=mean(double(sh2a(cell_indices)));
%      obj(c).green=mean(double(sktr(cell_indices)));  
      
       obj(c).volume=numel(cell_indices);
       obj(c).ktr=ktr;
       obj(c).signalCyt=signalCyt;
       obj(c).signalNuc=signalNuc;
       obj(c).area=area;
       obj(c).volNuc= volNuc;
       obj(c).volCyt= volCyt;
       obj(c).averageH2ANuc=averageH2ANuc;
       obj(c).averageH2ACyt=averageH2ACyt;
       obj(c).ratioH2ANucCyt=ratioH2ANucCyt;

       %obj(c).volume=1./sqrt(det(obj(c).W));
% 
%      [I1,I2,I3] = ind2sub(size(sktr),cell_indices);
%      [C,ia,ic]=unique([I1 I2],'rows');
%      obj(c).area=numel(ia);
%              
%      if(verbose)
%             colorIdx=round(log(ktr)*colorLevels./limLogAct+colorLevels);
%             scolor_r(cell_indices)=colors(colorIdx,1);
%             scolor_g(cell_indices)=colors(colorIdx,2);
%             scolor_b(cell_indices)=colors(colorIdx,3);
%      end
     
        
  end   
   
  if(verbose)
    scolor_r_proj=max(scolor_r,[],3);
    scolor_g_proj=max(scolor_g,[],3);
    scolor_b_proj=max(scolor_b,[],3);
    imshow(cat(3,scolor_r_proj,scolor_g_proj,scolor_b_proj));
  end

  opts.saveNameCyt = [opts.saveName 'Cyt'];
  opts.saveNameNuc = [opts.saveName 'Nuc'];
  
   %myScale.(fName)=[];
   myScale.(infoName).centers = vertcat(obj.m);
   %myScale.(infoName).red              =    [obj.red];
   %myScale.(infoName).green            =    [obj.green];
   myScale.(infoName).volume           =    vertcat(obj.volume);
   myScale.(infoName).areaNucleus      =    vertcat(obj.volNuc);
   myScale.(infoName).areaCyt          =    vertcat(obj.volCyt);
   
   myScale.(infoName).(opts.saveName)  =    vertcat(obj.ktr);
   myScale.(infoName).(opts.saveNameCyt)  = vertcat(obj.signalCyt);
   myScale.(infoName).(opts.saveNameNuc)  = vertcat(obj.signalNuc);

   myScale.(infoName).area             =    vertcat(obj.area);
   
   myScale.(infoName).averageH2ANuc           =     vertcat(obj.averageH2ANuc);
   myScale.(infoName).averageH2ACyt           =     vertcat(obj.averageH2ACyt);
   myScale.(infoName).ratioH2ANucCyt           =    vertcat(obj.ratioH2ANucCyt);
   
   
end
  