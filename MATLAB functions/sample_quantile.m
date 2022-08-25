function quantiles=sample_quantile(data,q)
    % sample_quantile Computes the q sample quantiles of data.
    %   quantiles=number_summary(data,q)
    %  
    %  data       should be a R-by-C numerical array containing the data.
    %  q          should be a 1-by-C numerical array containing proportions.
    %  quantiles  is a R-by-C numerical array containing the q-th quantiles.
    
    % check input arguments
    if ~isnumeric(data) || ~ismatrix(data); error('dat should be a R-by-C numeric array.'); end
    if size(data,1) < size(data,2); error('The data should be presented in columns.'); end
    if ~isnumeric(q) || ~isrow(q); error('q should be a 1-by-C numeric array.'); end
    if ~all(0 <= q & q <= 1); error('All elements in q should lie in the interval [0 1].'); end
    
    N = length(data);
    data_sorted = sort(data);
    q = sort(q);
    idx = N*q;
    if idx(1) == 0; idx(1) = 1; end
    
    quantiles = zeros(length(q),size(data,2));
    for i = 1:length(q)
        if mod(idx(i),1) ~= 0 || q(i) == 0 || q(i) == 1
            quantiles(i,:) = data_sorted(ceil(idx(i)),:);
        else
            quantiles(i,:) = mean([data_sorted(idx(i),:); data_sorted(idx(i)+1,:)]);
        end
    end
    
end
