function y = Eval_Spline(Tvec, Mat_S, tq)
N = length(Tvec)-1;
if tq >= Tvec(end)
    k = N;
elseif tq <= Tvec(1)
    k = 1;
else
    k = find(Tvec <= tq, 1, 'last');
    if k > N
        k = N;
    end
end
dx = tq - Tvec(k);
y = Mat_S(k,1)*dx^3 + Mat_S(k,2)*dx^2 + Mat_S(k,3)*dx + Mat_S(k,4);
end