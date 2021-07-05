#!/bin/bash
############################################
### Data, background (SS)
############################################

### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_VV_2016.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/NoSkim_SR_SS.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_VV_2016.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/NoSkim_SR_SS.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --nmax 80 &

### CF (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &



### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/NoSkim_SR_SS.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/NoSkim_SR_SS.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &

### CF (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &



### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/NoSkim_SR_SS.txt -n 150 --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_VV_2017.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_ttX.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/NoSkim_SR_SS.txt -n 150 --userflags RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --nmax 80 &

### CF (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunCF --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunFake --nmax 80 &



############################################
### Data, background (OS)
############################################

### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_SR_OS_2016.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Dilepton_SR_OS_2016.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS,RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i DoubleMuon:H -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR_2016H -y 2016 -i MuonEG:H -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &



### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_SR_OS_2017.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/Dilepton_SR_OS_2017.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS,RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_DoubleEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/2017_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &



### MC ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_SR_OS_2017.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS --nmax 80 &

### Fake from prompt (MC) ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/Dilepton_SR_OS_2017.txt -n 200 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/NoSkim_SR_OS.txt -n 150 --userflags RunOS,RunFake --nmax 80 &

### DATA ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS --nmax 80 &

### Fake (DATA) ###
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_DoubleMuon.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_EGamma.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/2018_MuonEG.txt -n 150 --skim SkimTree_Dilepton --userflags RunOS,RunFake --nmax 80 &



############################################
### Signal MC
############################################

###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Sch_MuMu.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Tch_MuMu.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Sch_EE.txt -n 30 --nmax 80 &
###python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/Signal_Tch_EE.txt -n 30 --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/DYTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/VBFTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/DYTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -l submitList/VBFTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/DYTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/VBFTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/DYTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -l submitList/VBFTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &

#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/DYTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/VBFTypeI_SS.txt -n 2 --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/DYTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2018 -l submitList/VBFTypeI_OS.txt -n 2 --userflags RunOS,RunSyst --nmax 80 &





############################################
### Test
############################################

#python/SKFlat.py -a HNtypeI_SR -y 2016 -i WZTo3LNu_powheg -n 150 --skim SkimTree_Dilepton --userflags RunSyst --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2016 -i WZTo3LNu_powheg -n 150 --skim SkimTree_Dilepton --nmax 80 &
#python/SKFlat.py -a HNtypeI_SR -y 2017 -i DYTypeI_SS_EE_M500 -n 2 --nmax 80 &
