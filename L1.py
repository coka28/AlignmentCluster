# Layer1 server script
# app: alignmentCluster

import os, shutil, pickle, subprocess
from time import sleep

os.chdir(os.path.expanduser("~"))
appDir = os.getcwd() + f'/apps/aligner'

# setup app directory structure
if not os.path.exists(f'{appDir}/projects'):
    os.mkdir(f'{appDir}/projects')

projLib = f'{appDir}/projects'

def setup(project):
    projDir = f'{projLib}/{project}'
    os.mkdir(projDir)
    # setup directory structure
    shutil.move(f'{projLib}/{project}.project',f'{projDir}/{project}.project')
    os.mkdir(f'{projDir}/clients')
    os.mkdir(f'{projDir}/registrations')
    os.mkdir(f'{projDir}/backup')
    os.mkdir(f'{projDir}/results')
    os.mkdir(f'{projDir}/src')
    # extract project file
    with open(f'{projDir}/{project}.project','rb') as tmp:
        p = pickle.load(tmp)
    ## syntax of p list
    titles = p[0]
    sequences = p[1]
    partNrs = p[2]
    partLens = p[3]
    params = p[4]
    buildstr = f'{params[0]}\n{params[1]}\n{params[2]}\n{params[3]}\n'
    for n,s,q,l in zip(titles,sequences,partNrs,partLens):
        with open(f'{projDir}/{n}.seq','w') as seqFile:
            seqFile.write(s)
        buildstr += f'{n}.seq;{q}*{l}\n'
    ## create initial WLs file
    WLs = 1
    for q in p[2]:
        WLs *= q
    WLs = [i for i in range(WLs)]
    with open(f'{projDir}/backup/openWLs','wb') as tmp:
        pickle.dump(WLs,tmp)
    with open(f'{projDir}/backup/pendingWLs','wb') as tmp:
        pickle.dump([],tmp)
    with open(f'{projDir}/backup/assignmentTimes','wb') as tmp:
        pickle.dump([],tmp)
    del(WLs)
    ## existence of project.build file is criterium for client to recognize project > do last
    with open(f'{projDir}/project.build','w') as buildFile:
        buildFile.write(buildstr)
    print(f'created new project "{project}"')

active = []
while True:
    for n in os.listdir(projLib):
        if '.project' in n:
            print(f'project registration request found... ',end='')
            projName = n[:n.find('.project')]
            setup(projName)
        
    for n in os.listdir(projLib):
        if '__done__' in n:
            try: active.remove(n.replace('__done__',''))
            except: pass
        elif n not in active:
            active.append(n)
            subprocess.Popen(f'python3 {appDir}/L2.py {n}', shell=True)
        
    sleep(2.47)
