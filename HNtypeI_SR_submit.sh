#!/bin/bash
############################################
### SR, CR
############################################

### MC ###
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_SR_2016.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/NoSkim_SR.txt -n 150 --nmax 80 &

### CF ###
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --nmax 80 &



### MC ###
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_SR_2017.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/NoSkim_SR.txt -n 150 --nmax 80 &

### CF ###
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &



### MC ###
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_SR_2017.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/NoSkim_SR.txt -n 150 --nmax 80 &

### CF ###
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake ###
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &

### DATA ###
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_SSHN --nmax 80 &
#python python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &



### Signal MC ###

###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Sch_MuMu.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Tch_MuMu.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Sch_EE.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Tch_EE.txt -n 30 --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/DYTypeI_2016.txt -n 2 --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/VBFTypeI.txt -n 2 --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/DYTypeI.txt -n 2 --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/VBFTypeI.txt -n 2 --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/DYTypeI.txt -n 2 --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/VBFTypeI.txt -n 2 --nmax 80 &
