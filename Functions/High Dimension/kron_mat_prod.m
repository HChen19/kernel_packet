function vt_ravel = kron_mat_prod(As, v)
    
    %Computes matrix-matrix multiplication between
    %matrix kron(As{1}, As{2}, ..., As{N}) and matrix
    %v without forming the full kronecker product.
    
    [nrows,ncols] = cellfun(@size,As);
    dims = [size(v,2),nrows];
    dims_flip = fliplr(dims);
    vt = reshape(v, dims_flip);
    for i = 1:length(As)
        vt = refold(As{i}*unfold(vt, i+1, dims_flip), i+1, dims);
    end
    vt_ravel = reshape(vt,[],size(v,2));
end

function ts_unfold = unfold(tens_flip,mode,dims_flip)

    %Unfolds tensor into matrix.
    %Parameters
    %----------
    %tens : ndarray, tensor with shape == dims
    %mode : int, which axis to move to the front
    %dims : list, holds tensor shape
    %Returns
    %-------
    %matrix : ndarray, shape (dims[mode], prod(dims[/mode]))
    dims = fliplr(dims_flip);
    if mode == 1
        ts_unfold = reshape(reshape(tens_flip,1,[]),[],dims(1)).';
    elseif mode == length(dims)
        ts_moveaxis = permute(tens_flip,[2:length(dims),1]);
        ts_unfold = reshape(reshape(ts_moveaxis,1,[]),[],dims(mode)).';
    else
        ts_moveaxis = permute(tens_flip,[1:length(dims)-mode,length(dims)-mode+2:length(dims),length(dims)-mode+1]);
        ts_unfold = reshape(reshape(ts_moveaxis,1,[]),[],dims(mode)).';
    end
end


function ts_refold_flip = refold(vec, mode, dims)

    %Refolds vector into tensor.
    %Parameters
    %----------
    %vec : ndarray, tensor with len == prod(dims)
    %mode : int, which axis was unfolded along.
    %dims : list, holds tensor shape
    %Returns
    %-------
    %tens : ndarray, tensor with shape == dims
    dims_flip = fliplr(dims);
    if mode == 1
        ts_refold_flip = reshape(reshape(vec.',1,[]),dims_flip);
    elseif mode == length(dims)
        tens = reshape(reshape(vec.',1,[]),fliplr([dims(mode),dims(1:mode-1)]));
        ts_refold_flip = permute(tens,[length(dims),1:length(dims)-1]);
    else
        % Reshape and then move dims[mode] back to its
        % appropriate spot (undoing the `unfold` operation).
        tens = reshape(reshape(vec.',1,[]),fliplr([dims(mode),dims(1:mode-1),dims(mode+1:length(dims))]));
        ts_refold_flip = permute(tens,[1:length(dims)-mode,length(dims),length(dims)-mode+1:length(dims)-1]);
    end
end