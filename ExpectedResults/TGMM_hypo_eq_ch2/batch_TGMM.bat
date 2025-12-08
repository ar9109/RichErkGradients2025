C:
cd C:\Users\User\Downloads\processDirectory\TGMM_Supplementary_Software_1_0\build\
for %%f in (C:\Users\User\Downloads\processDirectory\TGMM_hypo_eq_ch2\TGMMconfig\*.txt) do (

       cd nucleiChSvWshedPBC\Release
       start /wait ProcessStackBatchMulticore.exe %%f 0 1
       cd ../..
       cd Release
       start /wait TGMM.exe %%f 0 1
       cd ..
       del 
        
)
