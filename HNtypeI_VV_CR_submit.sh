#!/bin/bash
############################################
### VV,VG CR
############################################

### MC ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/Dilepton_VV_2016.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/Dilepton_VV_2016.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/NoSkim_VV_CR.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &




### MC ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/NoSkim_VV_CR.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &




### MC ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/NoSkim_VV_CR.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &




############################################
### Test
############################################

#python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -i WZTo3LNu_powheg -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
