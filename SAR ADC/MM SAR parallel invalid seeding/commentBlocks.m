function commentBlocks(start, stopnb)
%COMMENTBLOCKS Comment out sensing channels starting from start to end
    for i = start:stopnb
        basepath = 'MMADC/Sensing Channel ';
        path = append(basepath, num2str(i));
        set_param(path,'commented','on')
    end
%     for i = start:stopnb
%     basepath = 'MMADC/LNA ';
%     path = append(basepath, num2str(i));
%     set_param(path,'commented','on')
%     end
end