import subprocess

with open('stdout','w') as fp:
    subprocess.call('./CatchTestRunner --list-tags', shell=True, stdout=fp)

with open('stdout') as fp:
    tags = []
    for line in fp.readlines()[2::]:
        if line.strip():
            tag = [t for t in line.split(' ') if t][1].strip()
            tags.append(tag)

for tag in tags:
    subprocess.call('valgrind --tool=memcheck --error-limit=no --track-origins=yes --leak-check=full --show-leak-kinds=all ./CatchTestRunner '+tag, shell=True, cwd='.')
