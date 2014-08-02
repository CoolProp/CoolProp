
if __name__=='__main__':

    import shutil, os, sys, subprocess
    
    subprocess.check_call(['python','generate_headers.py'], shell = True, cwd=os.path.join('..','..','..','dev'))
    name = 'CoolProp5'
    
    import tempfile
    tmp = tempfile.mkdtemp(dir = '.')
    
    def collect():
        print('copying sources')
        shutil.copytree(os.path.join('..','..','..','src'),os.path.join(tmp,'src'))
        print('copying include')
        shutil.copytree(os.path.join('..','..','..','include'),os.path.join(tmp,'include'))
        print('copying externals')
        shutil.copytree(os.path.join('..','..','..','externals'),os.path.join(tmp,'externals'))
        print('copying python files')
        shutil.copytree(os.path.join('..',name),os.path.join(tmp,name))
        print('copying MANIFEST.in')
        shutil.copy2(os.path.join('..','MANIFEST.in'),os.path.join(tmp,'MANIFEST.in'))
        print('copying .version')
        shutil.copy2(os.path.join('..','..','..','.version'),os.path.join(tmp,'.version'))
        print('copying setup.py')
        shutil.copy2(os.path.join('..','setup.py'),os.path.join(tmp,'setup.py'))
        
    try:
       collect()
    except BaseException as B:
       #shutil.rmtree('collect')
       raise
    else:
        pass
       #shutil.rmtree('collect')