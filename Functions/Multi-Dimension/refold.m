function ts_refold_flip = refold(vec, mode, dims)

    %Refolds vector into tensor.
    %Parameters
    %----------
    %vec : ndarray, matrix with shape (dims[mode], prod(dims[/mode]))
    %mode : int, which axis was unfolded along.
    %dims : list, holds tensor shape
    %
    %Returns
    %-------
    %tensor_flip : ndarray, tensor with shape == dims_flip
    
    %dims_flip = fliplr(dims_col);
    if mode == 1
        ts_refold_flip = reshape(reshape(vec.',1,[]),[dims(end:-1:2), dims(mode)]);
    elseif mode == length(dims)
        tens = reshape(reshape(vec.',1,[]),fliplr([dims(mode),dims(1:mode-1)]));
        ts_refold_flip = permute(tens,[length(dims),1:length(dims)-1]);
    else
        % Reshape and then move dims[mode] back to its
        % appropriate spot (undoing the `unfold` operation).
        tens = reshape(reshape(vec.',1,[]),fliplr([dims(mode),dims(1:mode-1),dims(mode+1:length(dims))]));
        ts_refold_flip = permute(tens,[1:length(dims)-mode,length(dims),length(dims)-mode+1:length(dims)-1]);
    end
end%end refold function