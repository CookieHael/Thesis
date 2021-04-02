function uncommentBlocks(nbBlocks)
h=load_system("MMADC");
for i = 1:nbBlocks
    basepath = 'MMADC/Sensing Channel ';
    path = append(basepath, num2str(i));
    set_param(path,'commented','off')
end

% for i = 1:nbBlocks
%     basepath = 'MMADC/LNA ';
%     path = append(basepath, num2str(i));
%     set_param(path,'commented','off')
% end
end