function [myScale,obj] = assignKTR_wntaf(myScale,obj,svList,paths,opts)
%ANALYZESCALE Summary of this function goes here
%   Detailed explanation goes here
   
   inFolder_ch1=paths.inFolder_ch1;
   inFolder_ch2=paths.inFolder_ch2;
   basename=myScale.basename;
   
   infoName=opts.infoName;
   
   KTRsuffix=  paths.KTRsuffix;
   H2Asuffix=  paths.H2Asuffix;
   
   sktr= loadtiff([inFolder_ch1 basename KTRsuffix '.tif']);
   sh2a= loadtiff([inFolder_ch2 basename H2Asuffix '.tif']); 
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



  % generate full mask for entire image, edited 071322
    full_mask = zeros(size(sh2a),'logical');
    svList_idx_all = cat(2,obj.svIdx)'+1;
    cell_indices_all = cat(1,svList{svList_idx_all})+1;
    full_mask(cell_indices_all) = 1;

    % save nuclei mask for QC edited 072422
    if isfield(opts,'checksegmentation') && opts.checksegmentation
        mkdir([paths.plotFolder filesep basename])
        c_rand = randsample(numel(obj),10);
    end
    
  for c=1:numel(obj)    
     
     svIdx=obj(c).svIdx;
     cell_indices=[];
            
     for sv=1:numel(svIdx)
        cell_indices=vertcat(cell_indices,svList{svIdx(sv)+1}+1);
     end
         
     
     
     
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
         mini_full= full_mask(ymin:ymax,xmin:xmax,zmin:zmax); % edited 071322

         area = squeeze(sum(sum(max(minimask,[],3),1),2));
         if(isempty(area))
             area=0;
         end

         [~,zbest]=max(squeeze(mean(mean(double(minimask),1),2)));
         minimask =  minimask(:,:,zbest);
         mini_sh2a=  mini_sh2a(:,:,zbest);
         mini_sktr=  mini_sktr(:,:,zbest);
         mini_full = mini_full(:,:,zbest); % edited 071322

         % calculate cyt and nuc masks
         mm_nuc=imclose(minimask,strel('disk',3));
         %mm_nuc = imerode(mm_nuc,strel('disk',3)); %ALE 30jan20
         mm_cyt=logical(imdilate(minimask,strel('disk',3))-mm_nuc);
         mm_cyt_nonoverlap=mm_cyt&(~mini_full); % edited 071322

         % calculate signals
         volNuc=sum(double(mm_nuc(:)));
         volCyt=sum(double(mm_cyt(:)));
         signalCyt=sum(mini_sktr(mm_cyt(:)))./volCyt;
         signalNuc=sum(mini_sktr(mm_nuc(:)))./volNuc; 
         averageH2ANuc=sum(mini_sh2a(mm_nuc(:)))./volNuc;
         averageH2ACyt=sum(mini_sh2a(mm_cyt(:)))./volCyt;
         ratioH2ANucCyt=averageH2ANuc./averageH2ACyt;
         ktr=signalCyt./signalNuc;

        % edied 071322
        volCyt_nonoverlap = sum(double(mm_cyt_nonoverlap(:)));
        signalCyt_nonoverlap = sum(mini_sktr(mm_cyt_nonoverlap(:)))./volCyt_nonoverlap;
        averageH2ACyt_nonoverlap=sum(mini_sh2a(mm_cyt_nonoverlap(:)))./volCyt_nonoverlap;
        ratioH2ANucCyt_nonoverlap=averageH2ANuc./averageH2ACyt_nonoverlap;
        ktr_nonoverlap=signalCyt_nonoverlap./signalNuc;

        % erode mask for ktrQC, edited 020623
        if isfield(opts,'ktrQC')&&opts.ktrQC
            mm_nuc_erode = logical(imerode(minimask,strel('disk',2)));
            mm_cyt_erode = logical(mm_nuc-imerode(minimask,strel('disk',2)));
            volNuc_erode=sum(double(mm_nuc_erode(:)));
            volCyt_erode=sum(double(mm_cyt_erode(:)));
            signalCyt_erode=sum(mini_sktr(mm_cyt_erode(:)))./volCyt_erode;
            signalNuc_erode=sum(mini_sktr(mm_nuc_erode(:)))./volNuc_erode; 

            ktr_erode=signalCyt_erode./signalNuc_erode;

            px_nuc = mini_sktr(mm_nuc(:));
            px_cyt = mini_sktr(mm_cyt_nonoverlap(:));
            px_nuc_erode = mini_sktr(mm_nuc_erode);
        end

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

            % edited 071322
            figure;
            % cyto mask nonoverlap on nuclei
            imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
            hold on 
            %h = imshow(green); 
            %set(h,'AlphaData',double(max(mm_nuc,[],3))*0.1)
            h2 = imshow(red); 
            set(h2,'AlphaData',double(max(mm_cyt_nonoverlap,[],3))*0.5)
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
     
        % save nuclei mask for QC edited 072422
        if isfield(opts,'checksegmentation') && opts.checksegmentation
        
            if ismember(c,c_rand)
                f = figure(Visible="off");
                %  mask on nuclei 
                sizexy=[size(minimask,1) size(minimask,2)];
                red =   cat(3, ones(sizexy),zeros(sizexy),zeros(sizexy));
                imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
                hold on 
                h = imshow(red); 
                set(h,'AlphaData',bwperim(imdilate(double(max(mm_nuc,[],3)),strel('disk',1)))*0.8);
                %h2 = imshow(red); 
                %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.1)
                hold off
                
                exportgraphics(f,[paths.plotFolder filesep basename filesep 'obj_' num2str(c) '.png'],'BackgroundColor','black')
                close all;
            end
        end
         
%%     SAVE STUFF INTO OBJ!!!
     
%      obj(c).red=mean(double(sh2a(cell_indices)));
%      obj(c).green=mean(double(sktr(cell_indices)));  
    opts.saveNameCyt = [opts.saveName 'Cyt'];
    opts.saveNameNuc = [opts.saveName 'Nuc'];

       obj(c).(opts.saveName)=ktr;
       obj(c).(opts.saveNameCyt)=signalCyt;
       obj(c).(opts.saveNameNuc)=signalNuc;

       obj(c).volume=numel(cell_indices);
       obj(c).area=area;
       obj(c).volNuc= volNuc;
       obj(c).volCyt= volCyt;
       obj(c).averageH2ANuc=averageH2ANuc;
       obj(c).averageH2ACyt=averageH2ACyt;
       obj(c).ratioH2ANucCyt=ratioH2ANucCyt;

       % edited 071322
        obj(c).volCyt_nonoverlap= volCyt_nonoverlap;
        obj(c).averageH2ACyt_nonoverlap=averageH2ACyt_nonoverlap;
        obj(c).ratioH2ANucCyt_nonoverlap=ratioH2ANucCyt_nonoverlap;

        obj(c).([opts.saveNameCyt '_nonoverlap'])=signalCyt_nonoverlap;
        obj(c).([opts.saveName '_nonoverlap'])=ktr_nonoverlap;

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


        % erode mask for ktrQC, edited 020623
        if isfield(opts,'ktrQC')&&opts.ktrQC
            obj(c).mean_nuc = mean(px_nuc);
            obj(c).median_nuc = median(px_nuc);
            obj(c).CoV_nuc = std(px_nuc)./mean(px_nuc);

            obj(c).mean_nuccyt = mean([px_nuc;px_cyt]);
            obj(c).median_nuccyt = median([px_nuc;px_cyt]);
            obj(c).CoV_nuccyt = std([px_nuc;px_cyt])./mean([px_nuc;px_cyt]);

            obj(c).mean_nuc_erode = mean(px_nuc_erode);
            obj(c).median_nuc_erode = median(px_nuc_erode);
            obj(c).CoV_nuc_erode = std(px_nuc_erode)./mean(px_nuc_erode);
            obj(c).px_num_nuc_erode = numel(px_nuc_erode);

            obj(c).ktr_erode=ktr_erode;

        end
        
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
   myScale.(infoName).areaCyt_nonoverlap          =    vertcat(obj.volCyt_nonoverlap);

   
   myScale.(infoName).(opts.saveName)  =    vertcat(obj.(opts.saveName));
   myScale.(infoName).(opts.saveNameCyt)  = vertcat(obj.(opts.saveNameCyt));
   myScale.(infoName).(opts.saveNameNuc)  = vertcat(obj.(opts.saveNameNuc));

   myScale.(infoName).area             =    vertcat(obj.area);
   
   myScale.(infoName).averageH2ANuc           =     vertcat(obj.averageH2ANuc);
   myScale.(infoName).averageH2ACyt           =     vertcat(obj.averageH2ACyt);
   myScale.(infoName).ratioH2ANucCyt           =    vertcat(obj.ratioH2ANucCyt);
   
    % edited 071322
    myScale.(infoName).([opts.saveNameCyt '_nonoverlap'])= vertcat(obj.([opts.saveNameCyt '_nonoverlap']));
    myScale.(infoName).([opts.saveName '_nonoverlap'])= vertcat(obj.([opts.saveName '_nonoverlap']));
    
    myScale.(infoName).volCyt_nonoverlap= vertcat(obj.volCyt_nonoverlap);
    myScale.(infoName).averageH2ACyt_nonoverlap= vertcat(obj.averageH2ACyt_nonoverlap);
    myScale.(infoName).ratioH2ANucCyt_nonoverlap= vertcat(obj.ratioH2ANucCyt_nonoverlap);

    if isfield(opts,'ktrQC')&&opts.ktrQC
        myScale.(infoName).mean_nuc = vertcat(obj.mean_nuc);
        myScale.(infoName).median_nuc = vertcat(obj.median_nuc);
        myScale.(infoName).CoV_nuc = vertcat(obj.CoV_nuc);
        
        myScale.(infoName).mean_nuccyt = vertcat(obj.mean_nuccyt);
        myScale.(infoName).median_nuccyt = vertcat(obj.median_nuccyt);
        myScale.(infoName).CoV_nuccyt = vertcat(obj.CoV_nuccyt);
        
        myScale.(infoName).mean_nuc_erode = vertcat(obj.mean_nuc_erode);
        myScale.(infoName).median_nuc_erode = vertcat(obj.median_nuc_erode);
        myScale.(infoName).CoV_nuc_erode = vertcat(obj.CoV_nuc_erode);
        myScale.(infoName).px_num_nuc_erode = vertcat(obj.px_num_nuc_erode);
        
        myScale.(infoName).ktr_erode=vertcat(obj.ktr_erode);
    
    end

   
end
  