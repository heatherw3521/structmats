classdef FactorCache < handle
%FACTORCACHE  Stored factorizations of one HSS matrix (used by @hss).
%   A handle object kept in the hidden property H.factorcache of a root hss
%   object, so that H\b, minnorm and tikhonov can reuse factors although
%   hss is a value class. Not meant to be used directly.
%
%   ulv       factors of H itself (square or wide), from hss_ulvminnormsolve(H)
%   adjoint   factors of H' (tall H), used by H\b and B/H
%   weighted  struct('key', {weight}, 'F', factors of H*inv(L)) - last weight used
%   tikhonov  struct('key', {lambda, L, S}, 'F', ..., 'ix', ..., 'is', ...) - last call
    properties
        ulv = []
        adjoint = []
        weighted = []
        tikhonov = []
    end
    methods
        function reset(obj)
            obj.ulv = []; obj.adjoint = []; obj.weighted = []; obj.tikhonov = [];
        end
    end
end
