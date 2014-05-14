
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%All the include folders we need
include_string = [' -I',path_to_src]; 
% Add ",' -ldl'" to the above in order to compile with refprop support on 
% Linux: include_string = [' -I',path_to_src,' -ldl']; 
 
%List of files to be compiled to object files
bare_files = dir([path_to_src,'*.cpp']);
bare_files = cellfun(@(x) x, {bare_files.name}, 'uniformoutput', false)';

% A copy
main_files = bare_files;

%Append path to source to the list of the CoolProp main files
for i=1:size(main_files,1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files = main_files;
o_files = '';

for i=1:size(files,1)
	file = files{i,1};
	o_file = strrep(bare_files{i,1},'.cpp','.o');
    o_files = [o_files,' ',o_file];
    disp(file);
    eval(['mex -DEXTERNC -c', include_string,' -outdir . ',file])
end

%Build the MEX files
eval(['mex -v ', include_string,' Props.c ', o_files])
eval(['mex -v ', include_string,' HAProps.c ', o_files])

%Clean up - delete the obj files
delete('*.o')

%Quit MATLAB
%quit
