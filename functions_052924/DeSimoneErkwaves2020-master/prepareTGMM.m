function [] = prepareTGMM(fname,config_tpl,configfolder,datasource,datapathmac,respathmac,datapathwin,respathwin,bat_tpl,tgmmFolder,TGMMbuildFolder,bckgrd,tau)
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
arguments
    fname
    config_tpl
    configfolder
    datasource
    datapathmac
    respathmac
    datapathwin
    respathwin
    bat_tpl = [];
    tgmmFolder = [];
    TGMMbuildFolder = [];
    bckgrd = 10;
    tau = 8;
    

end

    basename=fname(1:end-4);

    copyfile([datasource filesep fname], [datapathmac filesep fname(1:end-4) '_t00.tif']);
    copyfile([datasource filesep fname], [datapathmac filesep fname(1:end-4) '_t01.tif']);

        
    folderpath=[respathmac basename];
%     folderpath=[respathmac basename '_bg' num2str(bckgrd) '_tau' num2str(tau)]; % modified for scanning parameters


    mkdir(folderpath);
    
    % edit TGMMconfig file
    fin  = fopen(config_tpl,'r');
    fout = fopen([configfolder filesep  basename '_config.txt'],'w');
%     fout = fopen([configfolder filesep  basename '_bg' num2str(bckgrd) '_tau' num2str(tau) '_config.txt'],'w'); % modified for scanning parameters

    idk=0;
    while ~feof(fin)
    
        idk=idk+1;
        s = fgetl(fin);

        if(idk==10)
            s=['imgFilePattern=' datapathwin basename '_t??'];

        end

        if(idk==13)
            s=['debugPathPrefix=' respathwin basename '\'];
%             s=['debugPathPrefix=' respathwin basename '_bg' num2str(bckgrd) '_tau' num2str(tau) '\']; % modified for scanning parameters

        end

        if(idk==22)
            s = ['backgroundThreshold=' num2str(bckgrd)];
        end

        if (idk==25)
            s = ['persistanceSegmentationTau=' num2str(tau)];
        end

        fprintf(fout,'%s\n',s);
            
     end
  
     fclose(fin);
     fclose(fout);    

     % edit TGMM .bat executable file
    fin  = fopen(bat_tpl,'r');
    fout = fopen([tgmmFolder filesep  'batch_TGMM.bat'],'w');
    
    s = TGMMbuildFolder(1:2);
    fprintf(fout,'%s\n',s);

    s = ['cd ' TGMMbuildFolder];
    fprintf(fout,'%s\n',s);

    idk=0;
    while ~feof(fin)
    
        idk=idk+1;
        s = fgetl(fin);

        if(idk==1)
            s=['for %%f in (' datapathwin(1:end-5) 'TGMMconfig\*.txt) do ('];
        end

        fprintf(fout,'%s\n',s);
            
     end
  
     fclose(fin);
     fclose(fout);    
        
end
    

