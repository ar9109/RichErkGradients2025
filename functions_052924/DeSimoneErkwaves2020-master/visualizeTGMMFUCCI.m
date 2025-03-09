function [myScale] = visualizeTGMMScale(myScale,paths,opts)
%ANALYZESCALE Summary of this function goes here
%   Detailed explanation goes here
   
      stackFolder=paths.stackFolder;
      basename=myScale.basename;
      suffix=paths.suffix;
        
      infoName=opts.infoName;
      
      
      %idx2choose=myScale.(opts.statName).idx2choose;
%     f_ell=myScale.(statName).fit_ellipse;
      
      centers=myScale.(infoName).centers;
      if(isfield(myScale.(infoName),'isGreen'))
        isGreen=myScale.(infoName).isGreen;
      else
          isGreen=[];
      end
      
      spath=[stackFolder basename suffix '.tif']
      s=loadtiff(spath); 
      s=max(s,[],3);
      
      figure;
      x=1:size(s,2); 
      %plot(x,(x-f_ell.X0_in)*tan(angle)+f_ell.Y0_in,'r','LineWidth',3);
      imshow(s,opts.adj); hold on;
      
      if(~isempty(isGreen))
        plot(centers(idx2choose(:)&isGreen(:),1),centers(idx2choose(:)&isGreen(:),2),'g.','MarkerSize',6); 
        plot(centers(idx2choose(:)&(~isGreen(:)),1),centers(idx2choose(:)&(~isGreen(:)),2),'r.','MarkerSize',6);
      else
        %plot(centers(idx2choose(:),1),centers(idx2choose(:),2),'r.','MarkerSize',2);
        plot(centers(:,1),centers(:,2),'r.','MarkerSize',10);

      end
       
%       th=0:0.05:2*pi;
%       bb=ell2cart(th,1,f_ell);
%       plot(bb(:,1),bb(:,2),'r');    

      %waitforbuttonpress;
end
   


