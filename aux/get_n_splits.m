function [n_split, part_start, part_end, n_part] = get_n_splits(ctrl,n_mc,n_meas,poolsize)
%& GET_N_SPLIT

%% currently not possible to evlauate poolsize on a worker
% % Get number of workers
% poolsize = matlabpool('size');
% if poolsize == 0
%     poolsize =1;
% end
% Get system memory
if ~isfield(ctrl,'sys')
    [~, memory] = system('cat /proc/meminfo | grep MemTotal');
    memory = str2double(memory(14:end-3))*1000; % in bytes
else
    memory = ctrl.sys.memory;
end

% 8 byte per double
% 4 empiric number to deal with the offset memory use (MIGHT BE IMPROVED)

n_split = ceil((8* poolsize * n_mc * n_meas*6)/ memory);

if n_split>1
    split = fix(n_meas/n_split);
    part_start(1) = 1;
    part_end(1)   = split;
    n_part(1)     = split;
    for  t = 2:n_split
        part_start(t) = part_end(t-1)+1;
        if t == n_split ,
            part_end(t)   = n_meas;
        else
            part_end(t) = part_end(t-1)+split;
        end;
        n_part(t) =  part_end(t) - part_start(t) +1;
    end
    
else
    part_start = 1;
    part_end   = n_meas;
end

