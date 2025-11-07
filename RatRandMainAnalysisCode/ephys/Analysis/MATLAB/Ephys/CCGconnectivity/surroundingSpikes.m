function [T3pre,T2pre,T1pre,T1post,T2post,T3post] = surroundingSpikes(R,T)
    
    T3preTmp  = [];
    T2preTmp  = [];
    T1preTmp  = [];
    T1postTmp = [];
    T2postTmp = [];
    T3postTmp = [];

    for i = 1:length(R)
        diffTimes = T-R(i);
        [~,idx]   = min(abs(diffTimes));
        
        if (idx > 4) && (idx < (length(T)-4))
            if (diffTimes(idx) < 0)
                T3preTmp  = [T3preTmp  T(idx-2)];
                T2preTmp  = [T2preTmp  T(idx-1)];
                T1preTmp  = [T1preTmp  T(idx)];
                T1postTmp = [T1postTmp T(idx+1)];
                T2postTmp = [T2postTmp T(idx+2)];
                T3postTmp = [T3postTmp T(idx+3)];
            elseif diffTimes(idx) > 0
                T3preTmp  = [T3preTmp  T(idx-3)];
                T2preTmp  = [T2preTmp  T(idx-2)];
                T1preTmp  = [T1preTmp  T(idx-1)];
                T1postTmp = [T1postTmp T(idx)];
                T2postTmp = [T2postTmp T(idx+1)];
                T3postTmp = [T3postTmp T(idx+2)];
            end
        end
    end

    % remove duplicates

    T3pre  = unique(T3preTmp);
    T2pre  = unique(T2preTmp);
    T1pre  = unique(T1preTmp);
    T1post = unique(T1postTmp);
    T2post = unique(T2postTmp);
    T3post = unique(T3postTmp);

end