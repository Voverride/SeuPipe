import os
import shutil

def makeDocs(folder, docRoot):
    docPath = os.path.join(docRoot, folder)
    if os.path.exists(docPath):
        shutil.rmtree(docPath)
    os.mkdir(docPath)
    files = os.listdir(folder)
    for file in files:
        if file.startswith('__'):
            continue
        fullPath = os.path.join(folder, file)
        if os.path.isdir(fullPath):
            makeDocs(fullPath, docRoot)
        elif file.endswith('.py'):
            mdPath = os.path.join(docRoot, fullPath.replace('.py', '.md'))
            with open(mdPath, 'w') as f:
                f.write('::: '+fullPath.replace('/', '.').replace('.py', ''))

def makeNav(folder, docRoot='docs', prefix='  '):
    docPath = os.path.join(docRoot, folder)
    header = prefix+'- '+folder.capitalize()+':\n'
    mark = False
    files = os.listdir(docPath)
    for file in files:
        if file.startswith('__'):
            continue
        constantPath = os.path.join(folder, file)
        fullPath = os.path.join(docRoot, constantPath)
        if os.path.isdir(fullPath):
            items = makeNav(constantPath, docRoot, prefix+'  ')
            if items!='':
                mark = True
                header+=items
        else:
            mark = True
            item = prefix+prefix+'- '+file.replace('.md', '')+': '+constantPath+'\n'
            header+=item
    return header if mark else ''

detected = ['controller', 'dataManager', 'api']
docRoot = 'docs'
for folder in detected:
    makeDocs(folder, docRoot)


config = 'mkdocs.yml'
templete = ''
with open(config, 'r') as f:
    for line in f:
        templete+=line
        if ' # auto generate' in line:
            break

nav = ''

for folder in detected:
    navItem = makeNav(folder)
    nav+=navItem

with open(config, 'w') as f:
    f.write(templete+'\n')
    f.write(nav)

print('auto generated docs successfully !')