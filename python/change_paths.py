#!/usr/bin/env python

import fnmatch,os,string,sys

def add_directory(dirname):
    parent = os.path.dirname(dirname)
    if not os.path.exists(parent) and parent != '':
        add_directory(parent)
    if not os.path.exists(dirname):
        print 'creating directory',dirname
        os.makedirs(dirname)
        os.system('svn add %s'%dirname)

def remove_directory(dirname):
    if os.path.exists(dirname):
        print 'removing directory',dirname
        os.system('svn rm %s'%dirname)
#        os.system('rm -rf %s'%dirname)

def move_file(source_filename,target_filename):
    print 'moving %s to %s' % (source_filename,target_filename)
    os.system('svn mv %s %s'%(source_filename,target_filename))    
#    os.system('mv %s %s'%(source_filename,target_filename))    

target_public = 'public'
target_private = 'private'

base_excludes = ['.svn','archive','tribble','integrationtests','settings',target_public,target_private]
private_paths = ['playground','oneoffprojects','oneoffs','archive','analysis','alignment','bwa','c','lua','matlab','perl','python','ruby','R','shell']
paths_to_trim = ['playground','oneoffprojects','oneoffs']
source_extensions = ['*.java','*.scala']

def intersect(a,b):
    return list(set(a) & set(b))

def is_source_file(file):
    for source_extension in source_extensions:
        if fnmatch.fnmatch(file,source_extension):
            return True
    return False

def modify_path(path):
    tokens = path.split('/')
    # compute proper target path: public or private?
    if(intersect(tokens,private_paths)):
        # kill private directory indicator only if private directory indicator is not the first element in the path.
        modified_tokens = [token for token in tokens if token not in paths_to_trim]
        return string.join(modified_tokens,'/').replace('./',target_private+'/')
    else:
        return path.replace('./',target_public+'/')

add_directory(target_public)
add_directory(target_private)

# just move archive wholesale; don't worry about processing at this point.
move_file('archive',target_private)

for root,dirs,files in os.walk('.'):
    # filter out non-processed files from root directory and base_excludes
    tokens = string.split(root,'/')
    if len(tokens) == 1:
        continue
    if len(intersect(tokens,base_excludes)) > 0:
        continue

    # handle file move
    for file in ["%s/%s"%(root,file) for file in files]:
        modified_path = modify_path(file)
        dirname = os.path.dirname(modified_path)
        add_directory(dirname)
        move_file(file,dirname)

# handle source code modification
for root,dirs,files in os.walk('.'):
    # process only public and private directories
    if not (root.startswith('./'+target_public) or root.startswith('./'+target_private)):
        continue
    for file in ["%s/%s"%(root,file) for file in files if is_source_file(file)]:
        f = open(file,'r')
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('package') or line.startswith('import'):
                tokens = line.split('.')
                if intersect(tokens,private_paths):
                    tokens = [token for token in tokens if token not in paths_to_trim]
                    modified_line = string.join(tokens,'.')
                    print "%s: '%s' => '%s'" % (file,line.rstrip(),modified_line.rstrip())
                    lines[i] = modified_line
        f.close()
        f = open(file,'w')
        f.writelines(lines)
        f.close()

for file in os.listdir('.'):
    if os.path.isdir(file) and not file in base_excludes and not file.startswith(target_public) and not file.startswith(target_private):
        remove_directory(file)
