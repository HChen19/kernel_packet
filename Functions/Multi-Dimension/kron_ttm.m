function kron_rs = kron_ttm(As, m)
    
    %Computes matrix-matrix multiplication between
    %matrix kron(As{1}, As{2}, ..., As{N}) and matrix
    %m by Tucker operator.
    
    [nrows,ncols] = cellfun(@size,As);
    dims = [size(m,2),nrows];
    dims_flip = fliplr(dims);
    m_ts = tensor(reshape(m, dims_flip));
    for i = 1:length(As)
        m_ts = ttm(m_ts,As{i},length(As)-i+1);
    end
    kron_rs = double(tenmat(m_ts,length(As)+1)');
end