function H = minus(H1, H2)
%MINUS  H1 - H2 for two HSS matrices on the same cluster tree (see plus).
if ~(isa(H1, 'hss') && isa(H2, 'hss'))
    error('hss:minus:operand', ['H - X is only supported for two HSS matrices; ' ...
        'for a dense X use full(H) - X.']);
end
H = plus(H1, uminus(H2));
end
