for %%f in (C:\Users\Ashley\Desktop\4Feb21_tph1b_H2A3b_ERK37a\TGMM_hypo_eq_ch2_repeat\TGMMconfig\*.txt) do (

       cd nucleiChSvWshedPBC\Release
       start /wait ProcessStackBatchMulticore.exe %%f 0 1
       cd ../..
       cd Release
       start /wait TGMM.exe %%f 0 1
       cd ..
       del 
        
)
