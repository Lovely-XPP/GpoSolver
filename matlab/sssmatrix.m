classdef sssmatrix  < handle
    % Sparse Symbolic Symmetric Matrix
    
    properties
        dim
        idx
        vals
        isnum
        nvals
    end
    
    methods
        function obj = sssmatrix(s)
            if (numel(s) ~= 1 || round(s) ~= s)
                error('Argument must be an integer');
            end
            
            obj.dim = s;
            obj.idx = zeros(0, 2);
            obj.vals = [];
            obj.isnum = [];
            obj.nvals = 0;
        end
                
        function  obj = push(obj, i, j, val)            
            obj.idx = [obj.idx; i, j];
            obj.vals = [obj.vals; val];
            obj.isnum = [obj.isnum; isnumber(val)]; 
            obj.nvals = obj.nvals + 1;
        end
        
        function [idx, vals, isnum] = get_values(obj)
            % Return non-symmetric non-zero values of the symmetric matrix
            idx = obj.idx;
            vals = obj.vals;
            isnum = obj.isnum;
        end        
        
        function A = get_real_matrix(obj, tvals_map)
            s = [obj.dim, obj.dim];
            is = obj.idx(:,1);
            js = obj.idx(:,2);
            ks1 = sub2ind(s, is, js);
            ks2 = sub2ind(s, js, is);
            
            vs = obj.vals;
            vs_vals = zeros(size(vs));
            
            for i = 1:numel(ks1)
                if (obj.isnum(i))
                    vs_vals(i) = double(vs(i));
                else
                    pars = symvar(vs(i));
                    pars_vals = arrayfun(@(x) tvals_map(char(x)), pars);
                    vs_vals(i) = double(subs(vs(i), pars, pars_vals));
                end
            end
            
            A = sparse(zeros(obj.dim));
            A(ks1) = vs_vals;
            A(ks2) = vs_vals;
        end
        
    end
    
end