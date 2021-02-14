# client main menu

import ftplib, os, sys, pickle, subprocess, traceback, shutil
from time import sleep

appDir = os.getcwd()
projLib = f'{appDir}/projects'

dn = 'bioinformatik-aufgabe.ddns.net'
acc = ('coka','coka')
ID = ''

try:
    with open('ID','r') as IDfile:
        ID = IDfile.read()
except:
    with open('ID','w') as IDfile:
        from random import randint
        IDchars = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',
                   'r','s','t','u','v','w','x','y','z','1','2','3','4','5','6','7','8',
                   '9','0']

        for i in range(12):
            nextChar = IDchars[randint(0,33)]
            ID += nextChar

        IDfile.write(ID)

def registerNew():
    try:
        print('          ----- Project builder -----')
        projNameAlph = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',
                        'r','s','t','u','v','w','x','y','z','1','2','3','4','5','6','7','8',
                        '9','0','-']
        
        print('Specify the name of the project:')
        print('(Naming rules: max length 100, not empty;\n letters (capital and small), numbers and dashes "-" only)')

        while True:
            valid = True
            pname = input('Project name: ')
            if pname == '':
                print('Empty input... enter project name!')
                valid = False
            else:
                for char in pname:
                    if char.lower() not in projNameAlph:
                        print('Project name consists of invalid symbols... enter valid project name')
                        valid = False
                if len(pname) > 100:
                    print('Project name is too long (max 100 symbols)... enter valid project name')
                    valid = False
                
            if valid:
                print('Trying connection to server...')
                while True:
                    try:
                        ftp = ftplib.FTP(dn)
                        ftp.login(*acc)
                        break
                    except:
                        print('Server not available! Retrying...')
                        sleep(1)

                print('Connection established... checking existing project names...')
                availableProjects = []
                ftp.retrlines('NLST apps/aligner/projects/',availableProjects.append)
                for name in availableProjects:
                    if f'{pname}' in name:
                        print(f'Project with name "{pname}" already exists... choose another name!')
                        ftp.quit()
                        valid = False
                        break

                
                if valid:
                    ftp.quit()
                    break

        print('Project name available... set project parameters:\n')

        while True:
            try:
                mt = abs(float(input('\t                 Match score: ')))
                mm = abs(float(input('\t          Substitution score: ')))
                pt = abs(float(input('\t           Insertion penalty: ')))
                minLen = int(input('\t                     Overlap: '))+1
                if mt == 0 or mm == 0 or pt == 0:
                    print('Zero is not a valid parameter! Enter valid parameters:')
                    raise Exception
                if minLen < 1:
                    print('Negative minimum length? Fuck off... provide valid parameters:')
                    raise Exception
                break
            except:
                print('You entered invalid parameters!')
                print('Valid formats for parameters: floats')
                print('Valid formats for overlap: natural number\n')
                print('Set valid project parameters:\n')
                pass
                

        print('\nSet titles for sequences! Press ENTER without input if you are done!')
        print('(Naming rules: max length 100, not empty;\n letters (capital and small), numbers and dashes "-" only)\n')
        
        titles = []
        i = 1

        while i<100000:
            valid = True
            titles.append(input(f'\tTitle - sequence nr {i}:\t'))
            if titles[-1] in titles[:-1]:
                print('Duplicate title... enter valid sequence title!\n')
                valid = False
            else:
                if titles[-1] == "":
                    del(titles[-1])
                    break
                else:
                    for char in titles[-1]:
                        if char.lower() not in projNameAlph:
                            print('Sequence title consists of invalid symbols... enter valid sequence title!\n')
                            valid = False
                            break
                    if len(titles[-1]) > 100:
                        print('Sequence title is too long... enter valid sequence title!\n')
                        valid = False
            
            if valid: i += 1
            else: del(titles[-1])
            
        print('\nProvide paths to sequence text files:\n')
        paths = []
        i = 0

        while i < len(titles): 
            paths.append(input(f'\tPath to {titles[i]}: \t'))
            try:
                with open(paths[-1],'r') as tmp: i += 1
            except:
                print(f'File not found: {paths[-1]} ... repeat input!\n')
                del(paths[-1])

        sequences = []

        for p in paths:
            with open(p,'r') as tmp:
                sequences.append(tmp.read().replace('\n',''))    
                
        print('\nSegmenting sequences... ',end='')
        targetsize = 3000000
        s = []                              # segment sizes
        s_ = targetsize**(1/len(sequences)) # target segment size
        q = []                              # number of segments per sequence with appendages
        ol = minLen-1                       # overlap between segments
        q0 = (len(sequences[0])-ol)/(s_-ol)
        eps = (s_*q0-len(sequences[0]))/(q0-1)
        seq_alphabet = []
        for seq in sequences:
            for char in seq:
                if char not in seq_alphabet:
                    seq_alphabet.append(char)
        appendage_index = 48
        lower_length_seqs = 0

        from math import ceil
        for i in range(len(sequences)):
            q_ = (eps-len(sequences[i])/(eps-s_))
            possible_qs = [k for k in range(int(q_*0.75),int(q_*1.1)+1)]
            tmp = [(len(sequences[i])+ol*(k-1))%k for k in possible_qs]
            for k in range(len(tmp)):
                for j in range(len(tmp)-1):
                    if tmp[j] < tmp[j+1]:
                        tmp[j],tmp[j+1] = tmp[j+1],tmp[j]
                        possible_qs[j],possible_qs[j+1] = possible_qs[j+1],possible_qs[j]
            q.append(possible_qs[0])
            targetLength = ceil((len(sequences[i])+ol*(q[-1]-1))/q[-1])*q[-1]-ol*(q[-1]-1)
            del(tmp)
            del(possible_qs)

            if len(sequences[i])<s_*1.5:
                lower_length_seqs += 1
                targetLength = len(sequences[i])
                q[-1] = 1
                if len(sequences)-lower_length_seqs > 0:
                    s_ = (targetsize/len(sequences[i]))**(1/(len(sequences)-lower_length_seqs))

            while len(sequences[i]) < targetLength:
                if (tmp:=chr(appendage_index)) not in seq_alphabet:
                    sequences[i] += tmp
                appendage_index += 1

            s.append(int((len(sequences[i])+ol*(q[-1]-1))/q[-1]))
        print('done')

        build = [titles,sequences,q,s,[mt,mm,pt,minLen]]

        print('Sending project to server... ',end='')
        with open(f'{pname}_tmp','wb') as tmp:
            pickle.dump(build,tmp)
        with open(f'{pname}_tmp','rb') as tmp:
            while True:
                try:
                    ftp = ftplib.FTP(dn)
                    ftp.login(*acc)
                    ftp.storbinary(f'STOR apps/aligner/projects/{pname}.project',tmp)
                    ftp.quit()
                    break
                except:
                    try: ftp.quit()
                    except: pass
                    sleep(5)
        os.remove(f'{pname}_tmp')
        print('success!\n')
        print('Returning to main menu:\n')
    except KeyboardInterrupt:
        return

def browseExisting():
    print('          ----- Project browser -----\n')
    from time import time
    while True:
        while True:
            try:
                ftp = ftplib.FTP(dn)
                ftp.login(*acc)
                tmp = []
                ftp.retrlines('NLST apps/aligner/projects',tmp.append)
                projects = []
                for name in tmp:
                    if '/projects/' in name and '.project' not in name and '__done__' not in name:
                        projects.append(name[name.find('/projects/')+10:])
                        
                for p in projects[::-1]:
                    projContents = []
                    ftp.retrlines(f'NLST apps/aligner/projects/{p}',projContents.append)
                    if f'apps/aligner/projects/{p}/project.build' not in projContents:
                        projects.remove(p)
                        
                ftp.quit()
                break
            except:
                traceback.print_exc()
                print('No connection to server! Retrying...')
                sleep(5)

        active = []
        
        for name in os.listdir(projLib):
            active.append(name)
        
        print('\nProjects on the server:\n')
        for i in range(len(projects)):
            tmp = f'(active)' if projects[i] in active else ''
            print(f'\t({i+1}) {projects[i]} {tmp}')
        print(f'\n\t({len(projects)+1}) Return to main menu')
        print()

        query = input('Choose project/option: ')
        while query not in [str(i) for i in range(1,len(projects)+2)]:
            query = input('Invalid input, repeat: ')

        i = int(query)-1
        if i == len(projects):
            print('\nReturning to main menu\n')
            return
        print(f'\n\t{projects[i]}:\n')
        if projects[i] in active:
            print('\t(1) Delete project from local workspace')
            print('\t(2) Display project information')
            print('\t(3) Return to project browser')

            query2 = input('\nChoose option: ')
            print()
            while query2 not in ['1','2','3']:
                query2 = input('Invalid input, repeat: ')

            if query2 == '1':
                try:
                    shutil.rmtree(f'{projLib}/{projects[i]}')
                    active.remove(projects[i])
                    print('Removed project. Returning to project browser:\n')
                except:
                    print('Project is currently opened, try later or kill workers!')
                    print('Returning to project browser:\n')

            if query2 == '2':
                while True:
                    print(f'Retrieving information about project "{projects[i]}":\n')
                    try:
                        ftp = ftplib.FTP(dn)
                        ftp.login(*acc)

                        with open(f'{projects[i]}.tmp','wb') as build:
                            ftp.retrbinary(f'RETR apps/aligner/projects/{projects[i]}/project.build',build.write)

                        clientDirs = []
                        ftp.retrlines(f'NLST apps/aligner/projects/{projects[i]}/clients',
                                      clientDirs.append)

                        doneWLs = 0
                        for c in clientDirs:
                            with open('WLdoneTMP','wb') as wltmp:
                                ftp.retrbinary(f'RETR {c}/res/done/done',wltmp.write)
                            with open('WLdoneTMP','rb') as wltmp:
                                wltemp = pickle.load(wltmp)
                            os.remove('WLdoneTMP')
                            doneWLs += wltemp

                        buildstr = open(f'{projects[i]}.tmp','r').read().split('\n')[:-1]
                        os.remove(f'{projects[i]}.tmp')
                        print(f'\tSequences:\n')
                        partNr = 1
                        partSize = 1
                        for seq in buildstr[4:]:
                            ol = int(buildstr[3])-1
                            q = int(seq[seq.find(".seq;")+5:seq.find("*")])
                            partNr *= q
                            s = int(seq[seq.find("*")+1:])
                            partSize *= s
                            print(f'\t{seq[:seq.find(".seq;")]}\t\t'+
                                  f'{q} parts of length {s}\t'+
                                  f'(total length: {q*s-(q-1)*ol})')
                        print()
                        print(f'\t Overall: {partNr} matrix parts total with size of {partSize} each')
                        print(f'\tFinished: {" "*(len(str(partNr))-len(str(doneWLs)))}{doneWLs} matrix parts ({round(doneWLs/partNr,4)*100}%)')
                        print(f'\n\t    Alignment parameters:')
                        print(f'\t             Match score:\t{buildstr[0]}')
                        print(f'\t      Substitution score:\t{buildstr[1]}')
                        print(f'\t                 Penalty:\t{buildstr[2]}')
                        print(f'\t                 Overlap:\t{int(buildstr[3])-1}')
                        
                        input('\nPress ENTER to return to project browser')
                        print()
                        try: ftp.quit()
                        except: pass
                        break
                    except Exception as e:
                        try: ftp.quit()
                        except: pass
                        traceback.print_exc()
                        print(f'Server not available ({e})... retrying')
                        sleep(5)

            if query2 == '3':
                print('Returning to project browser:\n')

        else:
            print('\t(1) Register for project')
            print('\t(2) Display project information')
            print('\t(3) Return to project browser')

            query2 = input('\nChoose option: ')
            print()
            while query2 not in ['1','2','3']:
                query2 = input('Invalid input, repeat: ')

            if query2 == '1':
                print(f'Retrieving project "{projects[i]}" and registering client...')
                while True:
                    try:
                        ftp = ftplib.FTP(dn)
                        ftp.login(*acc)

                        projDir = f'{projLib}/{projects[i]}'
                        os.mkdir(projDir)
                        shutil.copytree(f'{appDir}/bin/core',
                                        f'{projDir}/core')
                        shutil.copyfile(f'{appDir}/src/sysresources.txt',
                                        f'{projDir}/sysresources.txt')
                        print('Set up local workspace...')
                        with open(f'{projDir}/project.build','wb') as build:
                            ftp.retrbinary(f'RETR apps/aligner/projects/{projects[i]}/project.build',build.write)
                        print('Retrieved project build')
                        paths = []
                        with open(f'{projDir}/project.build','r') as build:
                            while '' not in paths:
                                paths.append(build.readline())
                        paths = [p[:p.find(';')] for p in paths[4:-1]]
                        for p in paths:
                            print(f'Retrieving {p}... ',end='')
                            with open(f'{projDir}/{p}','wb') as seqfile:
                                ftp.retrbinary(f'RETR apps/aligner/projects/{projects[i]}/{p}',seqfile.write)
                            print('done')
                        
                        registered = []
                        regDone = False
                        ftp.retrlines(f'NLST apps/aligner/projects/{projects[i]}/clients',
                                      registered.append)
                        for r in registered:
                            if ID in r: regDone = True

                        if not regDone:
                            with open('reg.tmp','wb') as tmp: pass
                            ftp.storbinary(f'STOR apps/aligner/projects/{projects[i]}/registrations/{ID}',
                                           open('reg.tmp','rb'))
                            os.remove('reg.tmp')
                            print('Client registration request sent to server...')
                        else:
                            print('Client already registered...')
                        
                        print('Project set up in local workspace... returning to project browser:\n')
                        ftp.quit()
                        break
                    except:
                        try: shutil.rmtree(projDir)
                        except:
                            traceback.print_exc()
                            pass
                        try: os.remove('reg.tmp')
                        except: pass
                        try: ftp.quit()
                        except: pass
                        print('Server not available... retrying')
                        sleep(5)

            if query2 == '2':
                while True:
                    print(f'Retrieving information about project "{projects[i]}":\n')
                    try:
                        ftp = ftplib.FTP(dn)
                        ftp.login(*acc)

                        with open(f'{projects[i]}.tmp','wb') as build:
                            ftp.retrbinary(f'RETR apps/aligner/projects/{projects[i]}/project.build',build.write)

                        clientDirs = []
                        ftp.retrlines(f'NLST apps/aligner/projects/{projects[i]}/clients',
                                      clientDirs.append)

                        doneWLs = 0
                        for c in clientDirs:
                            with open('WLdoneTMP','wb') as wltmp:
                                ftp.retrbinary(f'RETR {c}/res/done/done',wltmp.write)
                            with open('WLdoneTMP','rb') as wltmp:
                                wltemp = pickle.load(wltmp)
                            os.remove('WLdoneTMP')
                            doneWLs += wltemp

                        buildstr = open(f'{projects[i]}.tmp','r').read().split('\n')[:-1]
                        os.remove(f'{projects[i]}.tmp')
                        print(f'\tSequences:\n')
                        partNr = 1
                        partSize = 1
                        for seq in buildstr[4:]:
                            ol = int(buildstr[3])-1
                            q = int(seq[seq.find(".seq;")+5:seq.find("*")])
                            partNr *= q
                            s = int(seq[seq.find("*")+1:])
                            partSize *= s
                            print(f'\t{seq[:seq.find(".seq;")]}\t\t'+
                                  f'{q} parts of length {s}\t'+
                                  f'(total length: {q*s-(q-1)*ol})')
                        print()
                        print(f'\t Overall: {partNr} matrix parts total with size of {partSize} each')
                        print(f'\tFinished: {" "*(len(str(partNr))-len(str(doneWLs)))}{doneWLs} matrix parts ({round(doneWLs/partNr,4)*100}%)')
                        print(f'\n\t    Alignment parameters:')
                        print(f'\t             Match score:\t{buildstr[0]}')
                        print(f'\t      Substitution score:\t{buildstr[1]}')
                        print(f'\t                 Penalty:\t{buildstr[2]}')
                        print(f'\t                 Overlap:\t{int(buildstr[3])-1}')
                        
                        input('\nPress ENTER to return to project browser')
                        print()
                        try: ftp.quit()
                        except: pass
                        break
                    except Exception as e:
                        try: ftp.quit()
                        except: pass
                        print(f'Server not available ({e})... retrying')
                        sleep(5)
            if query2 == '3':
                print('Returning to project browser:\n')

def threadsMenu(fromMainMenu=True):
    print('Thread number menu:\n')
    print('Projects in local workspace:')
    projects = os.listdir(projLib)
    for i in range(len(projects)):
        print(f'\t({i+1}) {projects[i]}')
    print(f'\t({len(projects)+1}) Return to main menu')
    query = input('\nChoose project/option:')
    print()
    while query not in [str(i) for i in range(1,len(projects)+2)]:
        query = input('Invalid input, repeat: ')

    if query == str(len(projects)+1):
        if fromMainMenu: print('\nReturning to main menu\n')
        else: print('\n')
        return

    i = int(query)-1
    sysFile = open(f'projects/{projects[i]}/sysresources.txt','r').read()
    maxThreads = sysFile[sysFile.find('=')+1:sysFile.find('\n')]
    threadLine = sysFile.split('\n')[1]
    threadNr = threadLine[threadLine.find('threads=')+8:]
    print(f'\tNumber of logical processors: {maxThreads}')
    print(f'\t   Current number of threads: {threadNr}\n')
    query = input(f'\nSet new number of threads or leave empty to cancel: ')
    if query == '':
        print('\nReturning to thread menu:\n')
        threadsMenu()
        return
    
    while not query.isnumeric():
        query = input('Invalid input, repeat: ')

    sysFile = sysFile.replace(f'threads={threadNr}',f'threads={query}')
    with open(f'projects/{projects[i]}/sysresources.txt','w+') as tmp:
        tmp.write(sysFile)

    print(f'\nChanged number of workers for project {projects[i]}')
    print('Returning to thread menu:\n')
    threadsMenu()
    return

def startWorker():
    global workers, workerNr
    allProjects = os.listdir(projLib)
    print('Getting thread numbers...:')
    threads = {}
    wls = {}
    maxStrLen = 0
    with open('src/sysresources.txt','r') as sysFile:
        sysFile = sysFile.read()
        maxThreads = int(sysFile[sysFile.find('=')+1:sysFile.find('\n')])

    for proj in allProjects:
        with open(f'{projLib}/{proj}/sysresources.txt','r') as sysFile:
            text = sysFile.read().split('\n')[1]
            i = text.find('=')
            threadNr = int(text[i+1:])
            threads[proj] = threadNr
        if len(proj) > maxStrLen: maxStrLen = len(proj)
    allThreads = 0
    print()
    for n,i in zip(threads,range(1,len(threads)+1)):
        allThreads += threads[n]
        print(f'\t({i}) {n}: {" "*(maxStrLen-len(n))}{threads[n]}')
    print(f'\n\t({len(threads)+1}) Return to main menu')

    projectNr = input('\nChoose option/worker to start: ')
    while not projectNr.isnumeric() or int(projectNr) not in range(1,len(threads)+2):
        projectNr = input('Invalid input, repeat: ')

    if projectNr == str(len(threads)+1):
        print('\nReturning to main menu\n')
        return

    project = allProjects[int(projectNr)-1]
    query = input('\nDo you want to change thread settings now? (y/N) ').lower()
    while query not in ['y','n']:
        query = input('Invalid input, repeat: ').lower()
    if query == 'y':
        print()
        threadsMenu(False)

    print()
    input('Press ENTER to start worker or CTRL+C to abort!')

    print(f'\nStarting worker for project "{project}"...\n')
    
    cmd = ['python' if sys.platform.find('win') > -1 else 'python3',
            'main.py', project, str(workerNr)]
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = subprocess.SW_HIDE
    workers.append((project,subprocess.Popen(cmd, startupinfo=startupinfo),workerNr))
    workerNr += 1

    print('\nReturning to main menu... worker can continue in the background\n')

def stopWorker(stopAll=False):
    global workers, workerNr
    if not stopAll:
        if len(workers) == 0:
            print('Returning to main menu\n')
            return
        print('Active workers:\n')
        for i in range(len(workers)):
            print(f'\t({i+1}) {workers[i][0]}')
        print(f'\n\t({len(workers)+1}) Stop all workers')
        print(f'\n\t({len(workers)+2}) Return to main menu')
        query = input('\nChoose option/worker to stop: ')
        while query not in [str(i) for i in range(1,len(workers)+3)]:
            query = input('Invalid input, repeat: ')

        if query == f'{len(workers)+2}':
            print('\nReturning to main menu\n')
            return
    else: query = ''
    if query == f'{len(workers)+1}' or stopAll:
        print('\nStopping all workers.',end='')
        for i in range(len(workers)-1,-1,-1):
            workers[i][1].kill()
            while True:
                sleep(0.7)
                try:
                    shutil.rmtree(f'projects/{workers[i][0]}/wrktmp{workers[i][2]}')
                    break
                except: print('.',end='')
            del(workers[i])
        if not stopAll: print(' success! Returning to main menu\n')
        else: print(' workers stopped! See you soon!\n')
        return
    i = int(query)-1
    print('\nStopping worker.',end='')
    workers[i][1].kill()
    while True:
        sleep(0.7)
        try:
            shutil.rmtree(f'projects/{workers[i][0]}/wrktmp{workers[i][2]}')
            break
        except: print('.',end='')
    del(workers[i])
    print(' success! ',end='')
    stopWorker()
    return

global workers, workerNr
workers = []
workerNr = 0

def main(defOption=None):
    print('          ----- Alignment Cluster: Main menu -----\n')
    query = ""

    if defOption == 'start':
        startWorking(autoStart=True)
        sys.exit(0)
    
    print()
    while True:
        try:
            print('\t(1) Register new project')
            print('\t(2) Browse existing projects')
            print('\t(3) Set thread number for active projects')
            print('\t(4) Start worker')
            if workers:
                print('\t(5) Stop worker\n')
            else: print()
            print('Press CTRL+C to quit\n')
            query = input('Type number of preferred option: ')
            validInputs = ['1','2','3','4'] if not workers else ['1','2','3','4','5']
            while query not in validInputs:
                print('\nInvalid input...')
                query = input('Type number of preferred option: ')

            if query == '1':
                print('\n\n')
                registerNew()
                print()
            if query == '2':
                print('\n\n')
                browseExisting()
                print()
            if query == '3':
                print('\n\n')
                threadsMenu()
                print()
            if query == '4':
                print('\n\n')
                startWorker()
                print()
            if query == '5':
                print('\n\n')
                stopWorker()
                print()
                
        except KeyboardInterrupt:
            stopWorker(True)
            print()
            return

if not sys.argv[-1].isnumeric():
    main()
    input('\nPress ENTER to close window')
    sys.exit(0)

project = sys.argv[-2]
tmp = sys.argv[-1]
projDir = f'projects/{project}'
ID = ''

try:
    with open('ID','r') as IDfile:
        ID = IDfile.read()
except:
    with open('ID','w') as IDfile:
        from random import randint
        IDchars = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',
                   'r','s','t','u','v','w','x','y','z','1','2','3','4','5','6','7','8',
                   '9','0']

        for i in range(12):
            nextChar = IDchars[randint(0,33)]
            ID += nextChar

        IDfile.write(ID)

try: os.mkdir(f'{projDir}/wrktmp{tmp}')
except:
    print("Program didn't clean up properly on last exit\nTrying to clean up now... ",end='')
    try:
        shutil.rmtree(f'{projDir}/wrktmp{tmp}')
        print('success!')
    except Exception as e:
        print('error:',e)
        print('Close all other instances and retry!')
        sys.exit(0)
    
shutil.copyfile(f'{projDir}/project.build',f'{projDir}/wrktmp{tmp}/project.build')
shutil.copytree(f'{projDir}/core',f'{projDir}/wrktmp{tmp}/core')
for n in os.listdir(projDir):
    if '.seq' in n:
        shutil.copyfile(f'{projDir}/{n}',f'{projDir}/wrktmp{tmp}/{n}')

tmpDir = f'{projDir}/wrktmp{tmp}'
sysFilePath = os.getcwd()+f'/{projDir}'

os.chdir(tmpDir)
threadNr = None

while True:
    try:
        ftp = ftplib.FTP(dn)
        ftp.login(*acc)
        try:
            ftp.rename(f'apps/aligner/projects/{project}__done__',f'apps/aligner/projects/{project}__done__')
            break
        except: pass
        try:
            ftp.delete(f'apps/aligner/projects/{project}/clients/{ID}/inactive')
            sleep(10)
        except: pass
        revert = False
        with open(f'wl','wb') as wlFile:
            try:
                ftp.retrbinary(f'RETR apps/aligner/projects/{project}/clients/{ID}/WL',
                               wlFile.write)
                ftp.delete(f'apps/aligner/projects/{project}/clients/{ID}/WL')
            except:
                revert = True
        if revert:
            os.remove('wl')
            ftp.quit()
            raise Exception
        else:
            worklist = pickle.load(open('wl','rb'))
            os.remove('wl')
            ftp.quit()

        while len(worklist) > 0:
            with open(f'{sysFilePath}/sysresources.txt','r') as sysFile:
                tmp = sysFile.read().split('\n')[1]
                threadNr = int(tmp[tmp.find('=')+1:])
            worklist,nextUp = worklist[:-threadNr],worklist[-threadNr:]
            args = [str(s) for s in nextUp]
            cmd = ['java',f'core.Exec']+args
            worker = subprocess.Popen(cmd)
            worker.communicate()
            worker.wait()
            try:
                ftp = ftplib.FTP(dn)
                ftp.login(*acc)
                for file in os.listdir():
                    if 'alignments.' in file:
                        ftp.storbinary(f'STOR apps/aligner/projects/{project}/clients/{ID}/res/{file}',
                                       open(file,'rb'))
                        os.remove(file)
                ftp.quit()
            except:
                try: ftp.quit()
                except: pass

    except KeyboardInterrupt:
        for file in os.listdir():
            if 'alignments.' in file:
                os.remove(file)
        ftp.quit()
        break
    except:
        for file in os.listdir():
            if 'alignments.' in file:
                os.remove(file)
        sleep(5)
