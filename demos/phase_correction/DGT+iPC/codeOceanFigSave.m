function codeOceanFigSave(figHandle,folderName,fileName)
% This function is written only for Code Ocean.
if ~usejava('desktop') % detecting Code Ocean
    if exist(folderName,'dir')
        saveas(figHandle,[folderName '/' fileName])
    else
        warning 'Something wrong...?'
    end
end
end