function [matrix1, matrix2] = makeTertiary(matrix1, matrix2)

for i=1:size(matrix1,1)
    idx = find(matrix1(i,:));
    if (length(idx)>5)
        matrix1(i,idx(end-1)) = -1*matrix1(i,idx(end-1));
        matrix2(i,idx(end-1)) = -1*matrix2(i,idx(end-1));
    end
    if (length(idx)>10)
        matrix1(i,idx(end-2)) = -1*matrix1(i,idx(end-2));
        matrix2(i,idx(end-2)) = -1*matrix2(i,idx(end-2));
        matrix1(i,idx(end-3)) = -1*matrix1(i,idx(end-3));
        matrix2(i,idx(end-3)) = -1*matrix2(i,idx(end-3));
    end
end
        

end

