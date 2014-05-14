
import subprocess, os
import glob2 as glob

exports = ['-s','EXPORTED_FUNCTIONS=\"[\'_main\',\'_F2K\',\'_PropsSI\',\'_get_global_param_string\']\"']
optimization = '-D__ISLINUX__ -O2 -s DISABLE_EXCEPTION_CATCHING=0'

def compile_sources():
    for f in glob.glob(os.path.join('..','..','src','**','*.cpp')):
        
        call = [r'em++.bat',optimization,f,'-I../../include','-c','-DEXTERNC']+ exports
        print 'Calling:',' '.join(call)
        subprocess.check_output(' '.join(call), shell = True)

def link():
    call = [r'C:\Users\Belli\Downloads\emsdk-1.16.0-portable-64bit\emscripten\1.16.0\em++',optimization,'-v','-o','coolprop.js']+glob.glob('*.o')+['-DEXTERNC']  +  exports
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def closure_compiler():
    call = ['java','-Xmx1024m','-jar','compiler.jar','--js','coolprop.js','--js_output_file','coolprop2.js','--compilation_level','ADVANCED_OPTIMIZATIONS','--language_in','ECMASCRIPT5']
    print 'Using the closure compiler, this will take a while...   (from https://developers.google.com/closure/compiler/)'
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def cleanup():
    for file in glob.glob('*.o'):
        print 'removing',file
        os.remove(file)

def run():
    os.startfile('index.html')

if __name__=='__main__':
    compile_sources()
    link()
    #closure_compiler()
    cleanup()
    #run()