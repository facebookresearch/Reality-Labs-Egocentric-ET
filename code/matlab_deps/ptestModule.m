function [p] = ptestModule(group1,group2,iterations,tail)
% CSB - july 7 2021
% use: group1 > group2 by default, i.e., tail = 1
% tail = -1, group2>group1
% tail = 0, two-tailed test

if isempty(tail)
    tail = 1;
end

testStatDist = NaN(iterations,1);
for ii = 1:iterations
    concat1 = [group1(:); group2(:)];
    shuff1 = concat1(randperm(length(concat1)));
    newGroup1 = shuff1(1:length(group1));
    newGroup2 = shuff1(length(group1)+1:end);
    
    if tail == 1
        testStatDist(ii) = nanmean(newGroup1)-nanmean(newGroup2);
    elseif tail == -1
        testStatDist(ii) = nanmean(newGroup2)-nanmean(newGroup1);
    elseif tail == 0
        testStatDist(ii) = abs(nanmean(newGroup1)-nanmean(newGroup2));
    end
    
end

if tail == 0
    trueTestStat = abs(nanmean(group1) - nanmean(group2));
    p = (sum(abs(trueTestStat)<=abs(testStatDist)) +1)./(iterations+1); % two-tailed p-value; corrected for finite number of permutations
elseif tail == 1
    trueTestStat = nanmean(group1) - nanmean(group2);
    p = (sum(trueTestStat<=testStatDist) +1)./(iterations+1); % two-tailed p-value; corrected for finite number of permutations
elseif tail == -1
     trueTestStat = nanmean(group2) - nanmean(group1);
    p = (sum(trueTestStat<=testStatDist) +1)./(iterations+1); % two-tailed p-value; corrected for finite number of permutations
end


end