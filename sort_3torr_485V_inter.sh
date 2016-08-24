time ./ana.exe data/run00587sub00000.mid.gz data/run00589sub00000.mid.gz data/run00592sub00000.mid.gz data/run00593sub00000.mid.gz 
hadd -f output/3torr_485V_inter.root output/emma_ana_00587.root output/emma_ana_00589.root output/emma_ana_00592.root output/emma_ana_00593.root
#time ./ana.exe data/run00587sub00000.mid.gz data/run00588sub000* data/run00589sub00000.mid.gz data/run00592sub00000.mid.gz data/run00593sub00000.mid.gz 
#hadd -f output/3torr_485V_inter.root output/emma_ana_00587.root output/emma_ana_00588.root output/emma_ana_00589.root output/emma_ana_00592.root output/emma_ana_00593.root
#run 588 uses ADC and has too low statistics for calibration