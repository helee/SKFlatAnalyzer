#!/bin/bash
############################################
### VV,VG CR
############################################

### MC ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/Dilepton_VV_CR_2016.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_SSHN --nmax 80 &



### MC ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/Dilepton_VV_CR_2017.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_SSHN --nmax 80 &



### MC ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/Dilepton_VV_CR_2017.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/NoSkim_VV_CR.txt -n 150 --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_VV_CR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
