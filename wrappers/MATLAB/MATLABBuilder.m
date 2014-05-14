
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%All the include folders we need
include_string = [' -I',path_to_src];

if isempty(strfind(mexext(),'32'))
    mexopts_string = ' '
else
    mexopts_string = ' -f mexopts_w32.bat -DCONVENTION=__cdecl '
end
 
%List of files to be compiled to object files
main_files = dir([path_to_src,'*.cpp']);
main_files = cellfun(@(x) x, {main_files.name}, 'uniformoutput', false)';

%Append path to source to the list of the CoolProp main files
for i=1:size(main_files,1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files=[main_files];
o_files = '';
cpp_files = '';

for i=1:size(files,1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.obj');
	o_files = [o_files, ' ', o_file];
	cpp_files = [cpp_files, ' ',file];
    eval(['mex -c', include_string,mexopts_string,' -DCOOLPROP_LIB -outdir . ',file])
    disp(file)
end

%Build the MEX files
eval(['mex ', include_string,mexopts_string,' Props.c *.obj'])
eval(['mex ', include_string,mexopts_string,' PropsSI.c *.obj'])
eval(['mex ', include_string,mexopts_string,' HAProps.c *.obj'])

%Clean up - delete the obj files
delete('*.obj')
