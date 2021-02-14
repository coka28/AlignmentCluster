# Layer 0 server script:
    # starts python project scripts

import os, subprocess
from time import sleep

# setup directory structure
if not os.path.exists('apps'):
    os.mkdir('apps')

if not os.path.exists('src'):
    os.mkdir('src')

active = []

while True:
    # start apps
    var = os.listdir('apps')
    for n in var:
         if n not in active and os.path.exists(f'apps/{n}/L1.py'):
             active.append(n)
             subprocess.Popen(f'python3 {os.getcwd()}/apps/{n}/L1.py', shell=True)
             print(f'started app "{n}"')

    # extract installers
    var = os.listdir('src')
    for n in var:
        name = n[:n.find('.tar.gz')]
        os.system(f'tar -xzvf src/{n} -C {os.getcwd()}/apps/{name}')
        print(f'L0: created app "{n}" from installer')

    sleep(8.9)
