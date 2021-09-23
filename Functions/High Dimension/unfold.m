function ts_unfold = unfold(tens_flip,mode,dims_flip)

    %Unfolds tensor into matrix.
    %Parameters
    %----------
    %tens_flip : ndarray, tensor with shape == dims_flip
    %mode : int, which axis to move to the front
    %dims_flip : list, holds tensor shape with flip(dims)
    %
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
end%end unfold function