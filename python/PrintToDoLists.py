import os

def PutColor(a,color):
  return '${'+color+'}'+a+'${Color_Off}'

SKFlat_WD=os.environ['SKFlat_WD']
TMPDIR = SKFlat_WD+'/tmp/'
os.system('mkdir -p '+TMPDIR)
out = open(TMPDIR+'/ToDoLists.sh','w')

tags = [
'FIXME',
'TODO'
]

dirnames = [
'include',
'src'
]

for tag in tags:

  for dirname in dirnames:

    out.write('echo ""\n')
    out.write('echo "'+'###################"\n')
    out.write('echo "'+tag+' in '+dirname+'"\n')
    out.write('echo "'+'###################"\n')

    os.system('grep -r "'+tag+'" '+dirname+'/ > '+TMPDIR+'tmp_'+tag+'_'+dirname+'.txt')
    lines = open(TMPDIR+'tmp_'+tag+'_'+dirname+'.txt').readlines()
    os.system('rm '+TMPDIR+'tmp_'+tag+'_'+dirname+'.txt')
    for line in lines:
      line = line.strip('\n')
      words = line.split(':')

      out.write('printf "'+PutColor(words[0]+': ','BGreen')+'\t'+line.replace(words[0]+': ','')+'\\n"\n')
    out.write('\n')


out.write('printf "\\n\\n"')

out.close()
