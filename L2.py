# Layer 2 server script
# project worker


'''-.
+#_pü'-.....
ö*+...:(loop):..............................................
m}°:                                                        \
€>!:   1. register clients                                   \
&w^:   2. distribute WLs and add them to pending              \
j/6:   3. move results to results dir                          \
@²%:   4. remove timed-out from pending and re-open them       :§
#ß$:   5. check if done                                        /
6@y:   6. backup and call htmlUpdate                          /
µ<§:                                                         /
%$":......................................................../
%&"$%!§.-´´´´
€$"!.-´
'''


import sys, os, pickle, shutil, htmlTool
from time import time, sleep

os.chdir(os.path.expanduser("~"))
project = sys.argv[-1]
projDir = f'apps/aligner/projects/{project}'
clientsDir = f'{projDir}/clients'
regDir = f'{projDir}/registrations'
backupDir = f'{projDir}/backup'
resDir = f'{projDir}/results'


def registerClient(ID):
    print(f'{project}: \tregistering new client with ID {ID}')
    os.mkdir(f'{clientsDir}/{ID}')
    os.mkdir(f'{clientsDir}/{ID}/res')
    os.mkdir(f'{clientsDir}/{ID}/res/done')
    with open(f'{clientsDir}/{ID}/res/done/done','wb') as doneFile:
        pickle.dump(0,doneFile)

def passWLs():
    global openWLs, pendingWLs
    clients = os.listdir(clientsDir)
    for n in clients:
        if os.path.exists(f'{clientsDir}/{n}/inactive'):
            clients.remove(n)
    
    for n in clients:
        if os.path.exists(f'{clientsDir}/{n}/WL'):
            if time()-os.path.getmtime(f'{clientsDir}/{n}/WL') > 3600:
                print(f'{project}: \tclient {n} did not retrieve their workload... reassigning and setting to inactive')
                wl = pickle.load(open(f'{clientsDir}/{n}/WL','rb'))
                os.remove(f'{clientsDir}/{n}/WL')
                for w in wl:
                    if w in pendingWLs:
                        i = pendingWLs.index(w)
                        del(pendingWLs[i])
                        del(assignmentTimes[i])
                        openWLs.insert(0,w)
                with open(f'{clientsDir}/{n}/inactive','w') as tmp: pass
        else:
            tmp = min(min(128, int(len(openWLs)/len(clients))*4+1),len(openWLs))
            if tmp > 0: print(f'{project}: \tassigned {tmp} workloads to client {n}')
            wl = [openWLs.pop(0) for i in range(tmp)]
            with open(f'{clientsDir}/{n}/WL_tmp','wb') as tmp:
                pickle.dump(wl,tmp)
            for i in wl:
                pendingWLs.append(i)
                assignmentTimes.append(time())
            os.rename(f'{clientsDir}/{n}/WL_tmp',f'{clientsDir}/{n}/WL')

def moveResults():
    clientDirs = os.listdir(clientsDir)
    stored = 0
    for n in clientDirs:
        resFiles = os.listdir(f'{clientsDir}/{n}/res')
        resFiles.remove('done')
        with open(f'{clientsDir}/{n}/res/done/done','rb') as doneFile:
            doneWLs = pickle.load(doneFile)
        for m in resFiles:
            if os.path.getsize(f'{clientsDir}/{n}/res/'+m) == 0 and time()-os.path.getmtime(f'{clientsDir}/{n}/res/'+m)<60:
                pass
            else:
                resIndex = int(m[m.find('.')+1:])
                if resIndex in pendingWLs:
                    i = pendingWLs.index(resIndex)
                    alList = open(f'{clientsDir}/{n}/res/{m}','r').read().split('\n\n')
                    alList = [i for i in alList if i!='']
                    alignments = []
                    for al in alList:
                        tmp = al.split('\n')
                        tmp = [tuple(int(k) for k in tmp[j].split(';') if tmp[j]!='')
                               for j in range(len(tmp))]
                        alignments.append(tmp)
                    doneWLs += 1
                    with open(resDir+'/'+str(resIndex),'wb') as tmp:
                        pickle.dump(alignments,tmp)
                    os.remove(f'{clientsDir}/{n}/res/{m}')
                    # shutil.move(f'{clientsDir}/{n}/res/{m}',f'{resDir}/{m}')
                    stored += 1
                    del(pendingWLs[i])
                    del(assignmentTimes[i])
                else:
                    os.remove(f'{clientsDir}/{n}/res/{m}')
        with open(f'{clientsDir}/{n}/res/done/done','wb') as doneFile:
            pickle.dump(doneWLs,doneFile)
    if stored > 0: print(f'{project}: \tstored {stored} alignment parts in /results')

def reopen():
    reNr = 0
    for i in range(len(pendingWLs)-1,-1,-1):
        if time()-assignmentTimes[i] > 1800:
            openWLs.insert(0,pendingWLs[i])
            del(pendingWLs[i])
            del(assignmentTimes[i])
            reNr += 1
    if reNr > 0: print(f'{project}: \treopened {reNr} timed-out workloads')

def checkDone():
    if len(pendingWLs) + len(openWLs) == 0:
        print(f'{project}: \tproject finished')
        return True
    else: return False

def backup():
    with open(f'{backupDir}/openWLs','w+b') as tmp:
        pickle.dump(openWLs,tmp)
    with open(f'{backupDir}/pendingWLs','w+b') as tmp:
        pickle.dump(pendingWLs,tmp)
    with open(f'{backupDir}/assignmentTimes','w+b') as tmp:
        pickle.dump(assignmentTimes,tmp)
    print(f'{project}: \tcreated backup')

# load from backup
with open(f'{backupDir}/openWLs','rb') as tmp:
    openWLs = pickle.load(tmp)
with open(f'{backupDir}/pendingWLs','rb') as tmp:
    pendingWLs = pickle.load(tmp)
with open(f'{backupDir}/assignmentTimes','rb') as tmp:
    assignmentTimes = pickle.load(tmp)
print(f'{project}: \tretrieved data from project backup (open: {len(openWLs)}; pending: {len(pendingWLs)})')

backup_counter = 0
done = False
while not done:
    # 1.
    for ID in os.listdir(regDir):
        registerClient(ID)
        os.remove(f'{regDir}/{ID}')

    # 2.
    passWLs()

    # 3.
    moveResults()

    # 4.
    reopen()

    # 5.
    if checkDone(): done = True

    # 6.
    if backup_counter == 100 or done:
        backup()
        try: htmlTool.update()
        except: pass
        backup_counter = 0
    if done:
        os.rename(projDir,f'{projDir}__done__')
    backup_counter += 1
    sleep(1.74)
