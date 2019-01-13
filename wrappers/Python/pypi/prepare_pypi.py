
def collect(tmp):
    print('copying sources')
    shutil.copytree(os.path.join('..', '..', '..', 'src'), os.path.join(tmp, 'src'))
    print('copying include')
    shutil.copytree(os.path.join('..', '..', '..', 'include'), os.path.join(tmp, 'include'))
    print('copying externals')
    shutil.copytree(os.path.join('..', '..', '..', 'externals'), os.path.join(tmp, 'externals'))
    print('copying python files')
    shutil.copytree(os.path.join('..', name), os.path.join(tmp, name))
    print('copying MANIFEST.in')
    shutil.copy2(os.path.join('..', 'MANIFEST.in'), os.path.join(tmp, 'MANIFEST.in'))
    print('copying .version')
    shutil.copy2(os.path.join('..', '..', '..', '.version'), os.path.join(tmp, '.version'))
    print('copying setup.py')
    shutil.copy2(os.path.join('..', 'setup.py'), os.path.join(tmp, 'setup.py'))
    print('touching .build_without_cython')
    fp = open(os.path.join(tmp, '.build_without_cython'), 'w'); fp.close()
    print('touching .use_this_directory_as_root')
    fp = open(os.path.join(tmp, '.use_this_directory_as_root'), 'w'); fp.close()


if __name__ == '__main__':

    import shutil, os, sys, subprocess, glob

    subprocess.check_call('python generate_headers.py', shell=True, cwd=os.path.join('..', '..', '..', 'dev'), stdout=sys.stdout, stderr=sys.stderr)
    subprocess.check_call('python generate_constants_module.py', shell=True, cwd='..', stdout=sys.stdout, stderr=sys.stderr)
    for pyx in ['CoolProp.pyx', '_constants.pyx']:
        subprocess.check_call('cython --cplus ' + os.path.split(pyx)[1], shell=True, cwd=os.path.join('..', 'CoolProp'), stdout=sys.stdout, stderr=sys.stderr)
    name = 'CoolProp'

    # Make a temporary directory in this folder
    import tempfile
    tmp = tempfile.mkdtemp(dir='.')

    try:
        collect(tmp)

        # Make the source distro in this folder
        subprocess.check_call(' '.join(['python', 'setup.py', 'sdist'] + sys.argv[1::]), shell=True, cwd=tmp, stdout=sys.stdout, stderr=sys.stderr)

    except BaseException as B:

        shutil.rmtree(tmp)
        raise
    else:
        pass
        shutil.rmtree(tmp)
