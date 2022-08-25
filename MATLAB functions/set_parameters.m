function param_map=set_parameters(param_ids,param_values)
    % set_parameters Stores the given imaging parameters in a map container for easy access in the pipelines.
    %   param_map=set_parameters(param_ids,param_values)
    %   
    %   param_ids     should be a 1-by-C cell array of 1-by-C character arrays (cellstr), containing the parameter identifiers.
    %   param_values  should be a 1-by-C cell array, containing the corresponding parameter values.
    %   param_map     is a R-by-1 map container, which maps param_ids to param_values.
    %   
    %   See also read_ids, sort_participants, move_unselected, unzip_all, unzip_files.
    
    % check input arguments
    if ~iscellstr(param_ids) || ~isrow(param_ids); error('param_ids should be a 1-by-C cell array of 1-by-C character arrays (cellstr).'); end %#ok<ISCLSTR>
    if ~all(cellfun(@isrow,param_ids)); error('All elements of param_ids should be 1-by-C character arrays.'); end
    if ~iscell(param_values) || ~isrow(param_values); error('param_values should be a 1-by-C cell array.'); end
    if length(param_ids) ~= length(param_values); error('param_ids and param_values should be of the same length.'); end
    
    param_map = containers.Map(param_ids,param_values); % create map container
    
end
