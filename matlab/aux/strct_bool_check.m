function bool =  strct_flag_check(struct,field)
%% Function to check the existens and that the value is true
% A Geiges 09-13

bool = isfield(struct,field) && struct.(field);
