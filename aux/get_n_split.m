function n_split = get_n_split(ctrl,n_mc,n_meas,poolsize)
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

end
