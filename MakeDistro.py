from __future__ import print_function
import subprocess,os,shutil

#These should be paths to python executables that you want want use to build versions of CoolProp
PYTHONVERSIONS=['python.exe', #This is python 2.7 on my computer
                'c:\\python\\py27_x64\\python.exe',
                'c:\\python\\py32\\python.exe',
                'c:\\python\\py32_x64\\python.exe',
                'c:\\python\\py33\\python.exe',
                'c:\\python\\py33_x64\\python.exe',
                ]

if not os.path.exists('_deps'):
    os.mkdir('_deps')
        
def InstallPrereqs():
    """ Get the requirements for CoolProp """
    #Collect the source for Cython and put in _deps/cython-master
    import urllib,zipfile
    print('getting cython sources')
    urllib.urlretrieve('https://github.com/cython/cython/archive/master.zip', filename = 'master.zip')
    with zipfile.ZipFile('master.zip', 'r') as myzip:
        myzip.extractall(path='_deps')
    os.remove('master.zip')
    for python_install in PYTHONVERSIONS:
        for cwd in ['_deps/cython-master']:
            print(subprocess.check_output([python_install, 'setup.py', 'install'], cwd = cwd))
            
    
def PYPI():
    subprocess.call(['python','setup.py','sdist','upload'],cwd=os.path.join('wrappers','Python'))
    
def Source():
    print(subprocess.check_output(['python','setup.py','sdist','--dist-dir=../../dist_temp/Python'],shell=True,cwd=os.path.join('wrappers','Python')))

def DLL_and_Excel():
    """ Build a DLL using __stdcall calling convention """
    subprocess.check_output(['BuildDLL'],shell=True,cwd=os.path.join('wrappers','Excel'))
    subprocess.check_output(['BuildDLLx64'],shell=True,cwd=os.path.join('wrappers','Excel'))
    #Collect the zip file and p
    try:
        os.makedirs(os.path.join('dist_temp','Excel and DLL'))
    except os.error:
        pass
    
    shutil.copy2(os.path.join('CoolProp','CoolProp.h'),os.path.join('dist_temp','Excel and DLL','CoolProp.h'))    
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.dll'),os.path.join('dist_temp','Excel and DLL','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp_x64.dll'),os.path.join('dist_temp','Excel and DLL','CoolProp_x64.dll'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.xlam'),os.path.join('dist_temp','Excel and DLL','CoolProp.xlam'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.xla'),os.path.join('dist_temp','Excel and DLL','CoolProp.xla'))
    shutil.copy2(os.path.join('wrappers','Excel','TestExcel.xlsx'),os.path.join('dist_temp','Excel and DLL','TestExcel.xlsx'))
    shutil.copy2(os.path.join('wrappers','Excel','README.rst'),os.path.join('dist_temp','Excel and DLL','README.rst'))
    
def Octave():
    try:
        os.makedirs(os.path.join('dist_temp','Octave'))
        os.makedirs(os.path.join('dist_temp','Octave','3.6.1'))
        os.makedirs(os.path.join('dist_temp','Octave','3.6.2'))
        os.makedirs(os.path.join('dist_temp','Octave','3.6.4'))
    except os.error: pass
        
    subprocess.check_output(['OctaveBuilder.bat'],shell=True,cwd=os.path.join('wrappers','Octave'))
    shutil.copy2(os.path.join('wrappers','Octave','3.6.1','CoolProp.oct'),os.path.join('dist_temp','Octave','3.6.1','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','3.6.2','CoolProp.oct'),os.path.join('dist_temp','Octave','3.6.2','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','3.6.4','CoolProp.oct'),os.path.join('dist_temp','Octave','3.6.4','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','example.m'),os.path.join('dist_temp','Octave','example.m'))
    shutil.copy2(os.path.join('wrappers','Octave','README.rst'),os.path.join('dist_temp','Octave','README.rst'))
    
def Csharp():
    try:
        os.makedirs(os.path.join('dist_temp','C#'))
    except os.error: pass
        
    subprocess.check_output(['BuildCsharpDLL.bat'],shell=True,cwd=os.path.join('wrappers','C#'))
    shutil.copy2(os.path.join('wrappers','C#','readme.txt'),os.path.join('dist_temp','C#','readme.txt'))
    shutil.copy2(os.path.join('wrappers','C#','Csharp.7z'),os.path.join('dist_temp','C#','Csharp.7z'))
    
def MATLAB():
    print('MATLAB')
    try:
        os.makedirs(os.path.join('dist_temp','MATLAB'))
    except os.error: pass
        
    process = subprocess.Popen(['C:\\MATLAB_32bit\\bin\\matlab','-wait','-nodesktop','-nosplash','-nojvm','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    process = subprocess.Popen(['matlab','-nojvm','-nodesktop','-nosplash','-wait','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw64'),os.path.join('dist_temp','MATLAB','Props.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw64'),os.path.join('dist_temp','MATLAB','HAProps.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw32'),os.path.join('dist_temp','MATLAB','Props.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw32'),os.path.join('dist_temp','MATLAB','HAProps.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','PropsSI.mexw32'),os.path.join('dist_temp','MATLAB','PropsSI.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','PropsSI.mexw64'),os.path.join('dist_temp','MATLAB','PropsSI.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','README.rst'),os.path.join('dist_temp','MATLAB','README.rst'))
    shutil.copy2(os.path.join('wrappers','MATLAB','MATLAB_sample.m'),os.path.join('dist_temp','MATLAB','example.m'))
    
def Labview():
    print('Labview')
    import CoolProp
    version = CoolProp.__version__
    try:
        os.makedirs(os.path.join('dist_temp','Labview'))
    except os.error: pass
    
    process = subprocess.Popen(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','Labview'))
    process.wait()
    shutil.copy2(os.path.join('wrappers','Labview','CoolProp.dll'),os.path.join('dist_temp','Labview','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Labview','CoolProp.llb'),os.path.join('dist_temp','Labview','CoolProp.llb'))
    shutil.copy2(os.path.join('wrappers','Labview','CoolProp.vi'),os.path.join('dist_temp','Labview','CoolProp.vi'))
    shutil.copy2(os.path.join('wrappers','Labview','README.rst'),os.path.join('dist_temp','Labview','README.rst'))

def EES():
    print('EES')
    import CoolProp
    version = CoolProp.__version__
    try:
        os.makedirs(os.path.join('dist_temp','EES'))
    except os.error: pass
        
    subprocess.check_output(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    subprocess.check_output(['BuildDLF.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    subprocess.check_output(['BuildMSI.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    subprocess.check_output('C:\EES32\ees.exe C:\EES32\Userlib\CoolProp_EES\CoolProp_EES_Sample.ees /Solve',shell=True)
    
    shutil.copy2(os.path.join('wrappers','EES','Debug','CoolProp_EES_installer.msi'),os.path.join('dist_temp','EES','CoolProp_EES_installer.msi'))
    shutil.copy2(os.path.join('wrappers','EES','CoolProp.htm'),os.path.join('dist_temp','EES','CoolProp.htm'))
    shutil.copy2(os.path.join('wrappers','EES','README.rst'),os.path.join('dist_temp','EES','README.rst'))
    
def Javascript():
    print('Javascript')
    import CoolProp
    version = CoolProp.__version__
    try:
        os.makedirs(os.path.join('dist_temp','Javascript'))
    except os.error: pass
        
    subprocess.check_output('python build.py',shell=True,cwd=os.path.join('wrappers','Javascript'))

    shutil.copy2(os.path.join('wrappers','Javascript','index.html'),os.path.join('dist_temp','Javascript','index.html'))
    shutil.copy2(os.path.join('wrappers','Javascript','coolprop.js'),os.path.join('dist_temp','Javascript','coolprop.js'))
    shutil.copy2(os.path.join('wrappers','Javascript','README.rst'),os.path.join('dist_temp','Javascript','README.rst'))
    
def Java():
    print('Java')
    import CoolProp
    version = CoolProp.__version__
    try: 
        os.makedirs(os.path.join('dist_temp','Java','win32')) 
    except os.error as E:  print(E)
    try: 
        os.makedirs(os.path.join('dist_temp','Java','x64')) 
    except os.error as E:  print(E)
        
    subprocess.check_output('build_win32.bat',shell=True,cwd=os.path.join('wrappers','Java'))
    subprocess.check_output('build_x64.bat',shell=True,cwd=os.path.join('wrappers','Java'))

    subprocess.check_call(['7z','a','-r','dist_temp/Java/sources.zip','wrappers/Java/*.java'])
    shutil.copy2(os.path.join('wrappers','Java','win32','CoolProp.dll'),os.path.join('dist_temp','Java','win32','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Java','x64','CoolProp.dll'),os.path.join('dist_temp','Java','x64','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Java','README.rst'),os.path.join('dist_temp','Java','README.rst'))
    
def Python():
    print('Python')
    subprocess.check_output(['python','setup.py','install'],shell=True,cwd=os.path.join('wrappers','Python'))
    for python_install in PYTHONVERSIONS:
        subprocess.check_output([python_install,'setup.py','bdist','--format=wininst','--dist-dir=../../dist_temp/Python'],shell=True,cwd=os.path.join('wrappers','Python'))

def Maple():
    print('Maple')
    try:
        os.makedirs(os.path.join('dist_temp','Maple'))
    except os.error: pass
        
    process = subprocess.check_output(['BuildDLLx64.bat'],shell=True,cwd=os.path.join('wrappers','Maple'))
    
    listing = ['Analysis of a Refrigeration Cycle with CoolProp.mw','sample_file.mw','CoolProp_x64.dll']
    
    for file in listing:
        shutil.copy2(os.path.join('wrappers','Maple',file),os.path.join('dist_temp','Maple',file))
    
def Mathematica():
    print('Mathematica')
    try:
        os.makedirs(os.path.join('dist_temp','Mathematica'))
    except os.error: pass
        
    process = subprocess.check_output(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','Mathematica'))
    
    listing = ['README.rst','example.nb','CoolProp.dll']
    
    for file in listing:
        shutil.copy2(os.path.join('wrappers','Mathematica',file),os.path.join('dist_temp','Mathematica',file))
        
def Scilab():
    print('Scilab')
    try:
        os.makedirs(os.path.join('dist_temp','Scilab'))
    except os.error: pass
        
    process = subprocess.check_output(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','Scilab'))
    
    listing = ['README.rst','sample.sce','CoolProp_x64.dll','CoolProp.dll']
    
    for file in listing:
        shutil.copy2(os.path.join('wrappers','Scilab',file),os.path.join('dist_temp','Scilab',file))
    
def MathCAD():
    print('MathCAD')
    try:
        os.makedirs(os.path.join('dist_temp','MathCAD','Prime'))
    except os.error: pass
        
    process = subprocess.check_output(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','MathCAD'))
    process = subprocess.check_output(['BuildDLL.bat'],shell=True,cwd=os.path.join('wrappers','MathCAD','Prime'))
        
    shutil.copy2(os.path.join('wrappers','MathCAD','CoolPropMathcadWrapper.dll'),os.path.join('dist_temp','MathCAD','CoolPropMathcadWrapper.dll'))
    shutil.copy2(os.path.join('wrappers','MathCAD','CoolPropFluidProperties.xmcd'),os.path.join('dist_temp','MathCAD','CoolPropFluidProperties.xmcd'))
    shutil.copy2(os.path.join('wrappers','MathCAD','README.rst'),os.path.join('dist_temp','MathCAD','README.rst'))
    
    shutil.copy2(os.path.join('wrappers','MathCAD','Prime','CoolPropMathcadWrapper.dll'),os.path.join('dist_temp','MathCAD','Prime','CoolPropMathcadWrapper.dll'))
    shutil.copy2(os.path.join('wrappers','MathCAD','Prime','CoolPropFluidProperties.mcdx'),os.path.join('dist_temp','MathCAD','Prime','CoolPropFluidProperties.xmcd'))
    shutil.copy2(os.path.join('wrappers','MathCAD','Prime','README.rst'),os.path.join('dist_temp','MathCAD','Prime','README.rst'))
    
def Modelica():
    print('Modelica')
    try:
        os.makedirs(os.path.join('dist_temp','Modelica'))
    except os.error: pass
        
    process = subprocess.Popen(['BuildLIB-VS2008.bat'],shell=True,cwd=os.path.join('wrappers','Modelica')); process.wait()
    process = subprocess.Popen(['BuildLIB-VS2010.bat'],shell=True,cwd=os.path.join('wrappers','Modelica')); process.wait()
        
    shutil.copy2(os.path.join('wrappers','Modelica','README.rst'),os.path.join('dist_temp','Modelica','README.rst'))
    #shutil.copy2(os.path.join('wrappers','Modelica','src_modelica','CoolProp2Modelica.mo'),os.path.join('dist_temp','Modelica','CoolProp2Modelica.mo'))
    shutil.copy2(os.path.join('wrappers','Modelica','src','CoolPropLib.h'),os.path.join('dist_temp','Modelica','CoolPropLib.h'))
    shutil.copytree(os.path.join('wrappers','Modelica','bin','VS2008'),os.path.join('dist_temp','Modelica','VS2008'))
    shutil.copytree(os.path.join('wrappers','Modelica','bin','VS2010'),os.path.join('dist_temp','Modelica','VS2010'))
    
def UploadSourceForge():
    #Rename folder to version number
    import CoolProp
    try:
        shutil.copytree('dist_temp',CoolProp.__version__)
    except WindowsError: pass
    
    call_str = ['pscp','README.txt','ibell,coolprop@frs.sf.net:/home/pfs/project/c/co/coolprop/CoolProp/']
    print('Calling: '+' '.join(call_str))
    subprocess.check_output(call_str,shell=True)
    
    call_str = ['pscp','-r','-v',CoolProp.__version__,'ibell,coolprop@frs.sf.net:/home/pfs/project/c/co/coolprop/CoolProp/']
    print('Calling: '+' '.join(call_str))
    subprocess.check_output(call_str,shell=True)
    
    
    
#######  DOCUMENTATION STUFF ################
def Doxygen():
    # Open Doxyfile, and update the version number in the file
    lines = open('Doxyfile','r').readlines()
    import CoolProp
    for i in range(len(lines)):
        if lines[i].startswith('PROJECT_NUMBER'):
            line = lines[i].split('=')[0]+' = '+CoolProp.__version__+'\n'
            lines[i]=line
            break
    open('Doxyfile','w').write(''.join(lines))
    
    subprocess.check_output(['doxygen','Doxyfile'],shell=True)
    
def RunExamples():
    subprocess.check_output(['run_examples.bat'],shell=True,cwd='Web/examples')
    
def BuildDocs():
    
    # Inject the revision number into the docs main pages for the link
    lines = open('Web/_templates/index.html','r').readlines()
    import CoolProp
    languages = ['Python','Modelica','Labview','MATLAB','EES','Octave','Excel','C#','Java','Javascript','MathCAD','Maple','Mathematica','Scilab']
    for i in range(len(lines)):
        if (lines[i].find('http://sourceforge.net/projects/coolprop/files/CoolProp/') > -1
            and any([lines[i].find(a) > -1 for a in languages])
                ):
            oldVersion = lines[i].split('http://sourceforge.net/projects/coolprop/files/CoolProp/')[1].split('/',1)[0]
            lines[i] = lines[i][:].replace(oldVersion,CoolProp.__version__)
    open('Web/_templates/index.html','w').write(''.join(lines))
    
    shutil.rmtree(os.path.join('Web','_build'),ignore_errors = True)
    subprocess.check_output(['BuildCPDocs.bat'],shell=True,cwd='Web')
    
def UploadDocs():
    call_str = ['pscp','-r','-v','Web/_build/html/*.*','ibell@web.sourceforge.net:/home/groups/coolprop/htdocs']
    print('Calling: '+' '.join(call_str))
    subprocess.check_output(call_str, shell = True)
    
def Superpacks():
    
    import CoolProp
    subprocess.check_call(['git','archive','-o','dist_temp/CoolProp-'+CoolProp.__version__+'-source_code.zip','HEAD'])

    ## Windows superpack
    try:
        os.mkdir(os.path.join('dist_temp','windows_superpack'))
    except WindowsError:
        pass
    shutil.copy2(os.path.join('dist_temp/CoolProp-'+CoolProp.__version__+'-source_code.zip'),os.path.join('dist_temp','windows_superpack','CoolProp-'+CoolProp.__version__+'-source_code.zip'))
    
    for folder in ['Excel and DLL','Python','C#','Octave','MATLAB','EES','Labview','Modelica','Maple','Scilab','Mathematica','Java','Javascript','MathCAD']:
        shutil.copytree(os.path.join('dist_temp',folder), os.path.join('dist_temp','windows_superpack',folder))
    
    subprocess.check_call(['7z','a','-r','dist_temp/CoolProp-'+CoolProp.__version__+'-windows_superpack.zip','dist_temp/windows_superpack/*.*'])
    shutil.rmtree(os.path.join('dist_temp','windows_superpack'))
    
if __name__=='__main__':
    
#     InstallPrereqs()  #This is optional if you think any of the pre-reqs have been updated

    ## Always run the examples first.  This serves as an integration test of
    ## CoolProp
    RunExamples()

#     Python() # This one must be first to ensure that version gets rebuild properly
#     DLL_and_Excel()
#     Source()
#     Csharp()
#     Octave()
#     MATLAB()
#     Maple()
#     Mathematica()
#     Scilab()
#     EES()
#     Javascript()
#     Java()
#     MathCAD()
#     Labview()
#     Modelica()
#     Superpacks()

#     PYPI()
#     UploadSourceForge()
    
#     Doxygen()
#     BuildDocs()
#     UploadDocs()