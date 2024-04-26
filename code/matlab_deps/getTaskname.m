function [taskName] = getTaskname(fileName)
% extracts taskNames correctly given different file namings in soaring eagle dataset
% CSB Aug 29 2022

underscoreLoc = strfind(fileName,'_');
if contains(fileName,'PSNS')
    restOfString  = fileName((underscoreLoc(end)+1):end);
    temp = strfind(restOfString,' ')-1;
    taskName = restOfString(1:temp);
else
    restOfString  = fileName((underscoreLoc(1)+1):end);
    temp = strfind(restOfString,' ')-1;
    taskName = restOfString(1:temp);
    if contains(taskName,'_')
        underscoreLoc2 = strfind(taskName,'_');
        taskName = taskName(1:underscoreLoc2-1);
    end
end


end