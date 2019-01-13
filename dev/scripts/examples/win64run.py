import subprocess, os
from example_generator import *
import requests


def DownloadFile(url, outdir):
    fname = url.rsplit('/', 1)[1]
    r = requests.get(url)
    f = open(os.path.join(outdir, fname), 'wb')
    for chunk in r.iter_content(chunk_size=512 * 1024):
        if chunk:  # filter out keep-alive new chunks
            f.write(chunk)
    f.close()
    return


if __name__ == '__main__':

    P = Python()
    code = P.parse()
    P.write('Example.py', code)
    subprocess.call(r'python Example.py', stderr=sys.stderr, stdout=sys.stdout, cwd='.')

    if not os.path.exists('Octave'): os.mkdir('Octave')
    O = Octave()
    code = O.parse()
    O.write('Octave\Example.m', code)
    DownloadFile('http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/Octave/Octave3.8.2_Windows_32bit/CoolProp.oct', 'Octave')
    subprocess.call(r'c:\octave-3.8.2\bin\octave Example.m', stderr=sys.stderr, stdout=sys.stdout, cwd='Octave')

    if not os.path.exists('Java'): os.mkdir('Java')
    J = Java()
    code = J.parse()
    J.write('Java\Example.java', code)
    DownloadFile('http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/Java/platform-independent.7z', 'Java')
    subprocess.call(r'7z  platform-independent.7z', stderr=sys.stderr, stdout=sys.stdout, cwd='Java')

    DownloadFile('http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/Java/platform-independent.7z', 'Java')

    if not os.path.exists('Csharp'): os.mkdir('Csharp')
    C = Csharp()
    code = C.parse()
    C.write('Java\Example.cs', code)

    # ~ M = MATLAB()
    # ~ code = M.parse()
    #~ M.write('Example.m', code)
    #~ DownloadFile('http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/MATLAB/CoolPropMATLAB_wrap.mexw64', 'Octave')
