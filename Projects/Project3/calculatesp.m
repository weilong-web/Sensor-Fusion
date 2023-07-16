function [x, p] = calculatesp(h, R, sp, W)
    n = size(sp,1);
    x = zeros(n,1);
    for i=1:numel(W)
        x = x + h(sp(:,i)) * W(i);
    end   
    p = R; %zeros(n,n);
    for i=1:numel(W)
        p = p + (h(sp(:,i))-x)*(h(sp(:,i))-x).' * W(i);
    end
end