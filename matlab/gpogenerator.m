function gpogenerator(prob, params)
% GPOGENERATOR - problem C++ source code generator
%
% Given problem class parametrization PROB and parameter structure PARAMS,
% generates problem C++ source code
%
% GPOGENERATOR(PROB, PARAMS)
%
% The function takes two mandatory parameters. The first parameter is the 
% problem class definition PROBLEM, the second parameter is a data structure PARAMS.
%
% PARAMS structure can have five optional fields:
%
% filename    : A string containing the full path to the resulting C++ code.
%               For example, if set to 'problem', the gposolver will
%               generate two files problem.cpp and problem.h into the
%               current working directory. If set to '/home/user/p1.cpp',
%               files p1.cpp and p1.h will be generated into the /home/user/
%               directory. The default value is 'gpoproblem'.
% relax_order : Relaxation order of the generated LMI relaxation. 
%               If not set, gposolver will choose the minimal possible relaxation.
% classname   : A string containing the name of the generated C++ class. 
%               Default value is 'GpoProblem'.
% cext        : A string containing the file extension of the generated 
%               source file. Default value is 'cpp'.
% hext        : A string containing the file extension of the generated 
%               header file. Default value is 'h'.
%
%   See also GPORELAX, GPOSOLVE.

% (c) J. Heller, 2014-2015


% Test inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = default_param(params, 'filename', 'gpoproblem');
params = default_param(params, 'classname', 'GpoProblem');
params = default_param(params, 'hext', 'h');
params = default_param(params, 'cext', 'cpp');

% H file macro
[fpath, fname, ~] = fileparts(params.filename);
if (numel(fpath) > 0)
    fpath = [fpath filesep];
end

hlabel = sprintf('%s_H', upper(fname));
classname = params.classname;

% LMI relaxation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build relaxation if not already done so
if (~isfield(prob, 'Fs'))
    if (isfield(params, 'relax_order'))
        prob = gporelax(prob, params.relax_order);
    else
        prob = gporelax(prob);
    end        
end

% SDP Objective
obj_data = [prob.cobj_const, prob.cobj'];
[noflag, zoflag] = isnumber(obj_data);

% SDP constraints
[no_cons, no_mons] = size(prob.Fs);
no_con_data = 0;
no_nz_blocks = 0;
con_vars = [];

for i = 1:no_cons
    for j = 1:no_mons        
        con_vars = [con_vars; prob.Fs{i,j}.vals(~prob.Fs{i,j}.isnum)]; %#ok<AGROW>
        no_con_data = no_con_data + prob.Fs{i,j}.nvals;
        if (prob.Fs{i,j}.nvals ~= 0)
            no_nz_blocks = no_nz_blocks + 1;
        end
    end
end
no_con_vars = numel(con_vars);
    
% Export header file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hfile = [fname '.' params.hext];
hfname = [fpath hfile];
[fid, msg] = fopen(hfname, 'w');

fprintf(1, ['Generating ' hfile '... ']);

if (fid == -1)
    error(['Cannot open file (' hfname '): ' msg]);
end

export_h_file(fid, classname, hlabel, prob);

fclose(fid);
fprintf(1, 'done\n');

% Export source file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfile = [fname '.' params.cext];
cfname = [fpath cfile];
[fid, msg] = fopen(cfname, 'w');

fprintf(1, ['Generating ' cfile '... ']);

if (fid == -1)
    error(['Cannot open file (' cfname '): ' msg]);
end

no_aux = numel(keys(prob.tvars_map));

export_c_header(fid, hfile, prob);
export_tvars(fid, classname, prob);
export_polyobj(fid, classname, prob);
export_polycons(fid, classname, no_aux, prob);
export_objvars(fid, classname, noflag, prob);
export_convars(fid, classname, con_vars, no_aux);

export_probdims(fid, classname, prob);
export_objdata(fid, classname, noflag, zoflag, prob);

export_condata(fid, classname, no_con_data, no_con_vars, prob);

fclose(fid);
fprintf(1, 'done\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exp = default_param(exp, param_name, param_value)
    if (~isfield(exp, param_name))
        exp = setfield(exp, param_name, param_value);
    end
end

function export_c_header(fid, hfile, prob)
    header_str = ...
['/* \n' ...
 '   This file was automatically generated on %s\n' ...
 '   by GpoSolver, (c) 2014--2015, Jan Heller, hellej1@cmp.felk.cvut.cz,\n' ...
 '   Czech Technical University, Prague, Czech Republic \n' ...
 '\n' ...
 '   This solver takes:\n' ...
 '     %3d residual parameters: %s\n' ...
 '     %3d problem parameters : %s\n' ...
 '*/\n\n' ...
 '#include "%s"\n\n' ...
 'using namespace GpoSolver;\n' ...
 'using namespace Eigen;\n\n'];

    rpars_str = char(prob.rpars);
    if (numel(prob.rpars) > 1)
        rpars_str = rpars_str(9:end - 2);
    elseif (numel(prob.rpars) == 1)
        rpars_str = ['[' rpars_str ']'];
    end
    
    ppars_str = char(prob.ppars);        
    if (numel(prob.ppars) > 1)
        ppars_str = ppars_str(9:end - 2);
    elseif (numel(prob.ppars) == 1)
        ppars_str = ['[' ppars_str ']'];
    end

    fprintf(fid, header_str, datestr(now,'mmmm dd, yyyy HH:MM:SS'), ...
        numel(prob.rpars), rpars_str, numel(prob.ppars), ppars_str, hfile);
end

function export_h_file(fid, class_name, hlabel, prob)
    header_str = ...
['/* \n' ...
 '   This file was automatically generated on %s\n' ...
 '   by GpoSolver, (c) 2014--2015, Jan Heller, hellej1@cmp.felk.cvut.cz,\n' ...
 '   Czech Technical University, Prague, Czech Republic \n' ...
 '\n' ...
 '   This solver takes:\n' ...
 '     %3d residual parameters: %s\n' ...
 '     %3d problem parameters : %s\n' ... 
 '*/\n\n' ...
 '#ifndef %s\n' ...
 '#define %s\n\n' ...
 '#include <gposolver/gpoproblem.h>\n\n' ...
 'namespace GpoSolver {\n\n'];
 
    class_str = ...
['class %s : public GpoProblem {\n' ...
 '  private:\n' ...
 '    static const int num_blocks;\n' ...
 '    static const int num_cons;\n' ...
 '    static const int num_pcons;\n' ...
 '    static const int num_pvars;\n' ...
 '    static const int num_rpars;\n' ...
 '    static const int num_ppars;\n' ...
 '    static const int problem_size;\n' ...
 '    static const int relax_order;\n' ...
 '    static const int max_degree;\n' ...
 '    static const int rank_shift;\n' ... 
 '    static const int block_sizes[];\n' ...
 ' \n' ...
 '    static const int num_obj_data;\n' ...
 '    static const int num_obj_vars;\n' ... 
 '    static const objData obj_data[];\n' ... 
 '    static const objVar obj_vars[];\n' ...
 '\n' ...
 '    static const int num_con_data;\n' ...
 '    static const int num_con_vars;\n' ... 
 '    static const int num_con_blocks;\n' ...  
 '    static const conData con_data[];\n' ... 
 '    static const conVar con_vars[];\n' ...
 '    static const conBlock con_blocks[];\n' ...  
 '\n' ...
 '  public:\n' ...
 '  \n' ...
 '    int getNumBlocks(void) const {\n' ...
 '      return num_blocks;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumConstraints(void) const {\n' ...
 '      return num_cons;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumProblemVariables(void) const {\n' ...
 '      return num_pvars;\n' ...
 '    }\n' ...
 ' \n' ... 
 '    int getNumPolyConstraints(void) const {\n' ...
 '      return num_pcons;\n' ...
 '    }\n' ...
 ' \n' ...  
 '    int getNumObjectiveParams(void) const {\n' ...
 '      return num_rpars;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumConstraintsParams(void) const {\n' ...
 '      return num_ppars;\n' ...
 '    }\n' ...
 ' \n' ... 
 '    int getProblemSize(void) const {\n' ...
 '      return problem_size;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getRelaxationOrder(void) const {\n' ...
 '      return relax_order;\n' ...
 '    }\n' ...
 ' \n' ... 
 '    int getMaxDegree(void) const {\n' ...
 '      return max_degree;\n' ...
 '    }\n' ...
 ' \n' ... 
 '    int getRankShift(void) const {\n' ...
 '      return rank_shift;\n' ...
 '    }\n' ...
 ' \n' ...  
 '    const int * getBlockSizes(void) const {\n' ...
 '      return block_sizes;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getBlockSize(int i) const {\n' ...
 '      return block_sizes[i];\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumSdpObjectiveData(void) const {\n' ...
 '      return num_obj_data;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumSdpObjectiveVariables(void) const {\n' ...
 '      return num_obj_vars;\n' ...
 '    }\n' ...   
 ' \n' ...
 '    const objData * getSdpObjectiveData(void) const {\n' ...
 '      return obj_data;\n' ...
 '    }\n' ...
 ' \n' ...
 '    const objVar * getSdpObjectiveVariables(void) const {\n' ...
 '      return obj_vars;\n' ...
 '    }\n' ...
 ' \n' ...
 '    const conData * getSdpConstraintsData(void) const {\n' ...
 '      return con_data;\n' ...
 '    }\n' ...
 ' \n' ...
 '    const conVar * getSdpConstraintsVariables(void) const {\n' ...
 '      return con_vars;\n' ...
 '    }\n' ...
 '    int getNumSdpConstraintsData(void) const {\n' ...
 '      return num_con_data;\n' ...
 '    }\n' ...
 ' \n' ...
 '    int getNumSdpConstraintsVariables(void) const {\n' ...
 '      return num_con_vars;\n' ...
 '    }\n' ...  
 ' \n' ...
 '    void updateAuxiliaries(const double *, double *) const;\n' ...
 '    double evalPolyObjective(const double *, const double *, const double *) const;\n' ...
 '    void evalPolyConstraints(const double *, const double *, double *) const;\n' ...
 '    void evalSdpObjectiveVariables(const double *, const double *, double **) const;\n' ...
 '    void evalSdpConstraintsVariables(const double *, const double *, double **) const;\n' ...
 ];

    footer_str = ...
['\n' ...
 '}; // %s\n', ...
 '\n' ...
 '} // GPOSolver namespace \n' ...
 '\n' ...
 '#endif // %s\n\n'];

    rpars_str = char(prob.rpars);
    if (numel(prob.rpars) > 1)
        rpars_str = rpars_str(9:end - 2);
    elseif (numel(prob.rpars) == 1)
        rpars_str = ['[' rpars_str ']'];
    end
    
    ppars_str = char(prob.ppars);        
    if (numel(prob.ppars) > 1)
        ppars_str = ppars_str(9:end - 2);
    elseif (numel(prob.ppars) == 1)
        ppars_str = ['[' ppars_str ']'];
    end
   
    fprintf(fid, header_str, datestr(now,'mmmm dd, yyyy HH:MM:SS'), ...
        numel(prob.rpars), rpars_str, numel(prob.ppars), ppars_str, hlabel, hlabel);
    
    %fprintf(fid, header_str, datestr(now,'mmmm dd, yyyy HH:MM:SS'), 
    fprintf(fid, class_str, class_name);    
    fprintf(fid, footer_str, class_name, hlabel); 
end

function export_probdims(fid, class_name, prob)
    sz = size(prob.Fs);
    
    fprintf(fid, '\n');
    fprintf(fid, 'const int %s::num_blocks = %d;\n', class_name, sz(1));
    fprintf(fid, 'const int %s::num_cons = %d;\n', class_name, sz(2) - 1);
    fprintf(fid, 'const int %s::num_pcons = %d;\n', class_name, numel(prob.ccons));     
    fprintf(fid, 'const int %s::num_pvars = %d;\n', class_name, numel(prob.vars));    
    fprintf(fid, 'const int %s::num_rpars = %d;\n', class_name, numel(prob.rpars));
    fprintf(fid, 'const int %s::num_ppars = %d;\n', class_name, numel(prob.ppars));        
    fprintf(fid, 'const int %s::problem_size = %d;\n', class_name, sum(prob.cs_sizes));    
    fprintf(fid, 'const int %s::relax_order = %d;\n', class_name, prob.order);        
    fprintf(fid, 'const int %s::max_degree = %d;\n', class_name, prob.max_degree);        
    fprintf(fid, 'const int %s::rank_shift = %d;\n', class_name, prob.rankshift);          
    fprintf(fid, 'const int %s::block_sizes[] = {', class_name);    
    
    for i = 1:(numel(prob.cs_sizes) - 1)
        fprintf(fid, '%d, ', prob.cs_sizes(i));
    end
    fprintf(fid, '%d};\n\n', prob.cs_sizes(end));
end

function export_objdata(fid, class_name, sflag, zflag, prob)    
    obj_data = [prob.cobj_const, prob.cobj'];       
    no_data = numel(obj_data);
    nonzeros = ~zflag;
    sflag = ~sflag;

    fprintf(fid, 'const int %s::num_obj_data = %d;\n', class_name, sum(nonzeros));
    fprintf(fid, 'const int %s::num_obj_vars = %d;\n', class_name, sum(sflag));    

    if (any(nonzeros))    
        idx = 1:no_data;
        idx = idx(nonzeros);
        k = idx(end);
        c = 1;
        
        fprintf(fid, 'const %s::objData  %s::obj_data[] = {\n', class_name, class_name);
        for i = idx            
            if (c > 6)
                fprintf(fid, '\n');
                c = 1;
            end
            c = c + 1;
            
            if (sflag(i))
                data_str = '0.0';
            else
                data_str = char(obj_data(i));
            end
            
            if (i < k)
                fprintf(fid, '  {%d, %s},', i - 1, data_str);
            else
                fprintf(fid, '  {%d, %s}\n};\n', i - 1, data_str);
            end
        end    
    else
        fprintf(fid, 'const %s::objData * %s::obj_data[] = {};\n', class_name, class_name);
    end
    
    if (any(sflag))
        idx = 1:no_data;
        idx = idx(sflag);
        k = idx(end);
        c = 1; 

        fprintf(fid, 'const %s::objVar %s::obj_vars[] = {\n', class_name, class_name);
        for i = idx            
            if (c > 6)
                fprintf(fid, '\n');
                c = 1;
            end
            c = c + 1;
            
            if (i < k)
                fprintf(fid, '  {%d},', i - 1);
            else
                fprintf(fid, '  {%d}\n};\n', i - 1);
            end
        end    
    else
        fprintf(fid, 'const %s::objVar %s::obj_vars[] = {};\n', class_name, class_name);
    end
    
    fprintf(fid, '\n');
end

function export_objvars(fid, class_name, sflag, prob)  
    fprintf(fid, '\nvoid %s::evalSdpObjectiveVariables(const double *rpars, const double *ppars, double **vars) const {\n', class_name); 
    obj_data = [prob.cobj_const, prob.cobj'];
    obj_data = obj_data(~sflag);
    
    nrpars = getsyms('rpars%d', numel(prob.rpars)); 
    nppars = getsyms('ppars%d', numel(prob.ppars)); 
    obj_data = subs(obj_data, prob.rpars, nrpars);
    obj_data = subs(obj_data, prob.ppars, nppars);
        
    if (isempty(obj_data))
        fprintf(fid, '}\n'); 
        return;
    end        
    
    % generate C code
    [fid_c, tmp_file] = genccode(obj_data', 1);        
    
    line = fgets(fid_c);
    while (ischar(line))
        line = regexprep(line, 'rpars(\d+)', 'rpars[$1]');
        line = regexprep(line, 'ppars(\d+)', 'ppars[$1]');
        line = regexprep(line, '\s+(t\d+ = .*)', '  double $1');
        line = regexprep(line, '\s+A0\[(\d+)\]\[\d+\] = (.*)', '  *(vars[$1]) = $2');     
        
        if ((numel(line) > 150) && (strcmp(line(1:8), '  *(vars')))
            split_line(line, fid, '');
        else
            fprintf(fid, line);         
        end
        
        line = fgets(fid_c);
    end
    
    fdelete(fid_c, tmp_file);    
 
    fprintf(fid, '}\n');   
end

function export_condata(fid, class_name, no_con_data, no_con_vars, prob)
    % Size
    fprintf(fid, 'const int %s::num_con_data = %d;\n', class_name, no_con_data);
    fprintf(fid, 'const int %s::num_con_vars = %d;\n', class_name, no_con_vars);
    [no_cons, no_mons] = size(prob.Fs);
    
    % Data    
    if (no_con_data)
        fprintf(fid, 'const %s::conData %s::con_data[] = {\n', class_name, class_name);
        
        l = 1;
        for i = 1:no_cons
            for j = 1:no_mons
                for k = 1:prob.Fs{i,j}.nvals
                    if (~mod(l, 6))
                        fprintf(fid, '\n');
                    end
                    
                    if (~prob.Fs{i,j}.isnum(k))
                        data_str = '0.0';
                    else
                        data_str = [char(prob.Fs{i,j}.vals(k)) '.0'];
                    end
                    
                    if (l < no_con_data)
                        fprintf(fid, '  {%d, %d, %d, %d, %s},', i - 1, j - 1, ...
                            prob.Fs{i,j}.idx(k,1) - 1, prob.Fs{i,j}.idx(k,2) - 1, data_str);
                    else
                        fprintf(fid, '  {%d, %d, %d, %d, %s}\n};\n', i - 1, j - 1, ...
                            prob.Fs{i,j}.idx(k,1) - 1, prob.Fs{i,j}.idx(k,2) - 1, data_str);
                    end
                    l = l + 1;
                end
            end
        end
    else
        % Bogus value to fool msvc
        fprintf(fid, 'const %s::conData %s::con_data[] = {-1};\n', class_name, class_name);        
    end
    
    % Variables
    if (no_con_vars)
        fprintf(fid, 'const %s::conVar %s::con_vars[] = {\n', class_name, class_name);
        
        l = 1;
        for i = 1:no_cons
            for j = 1:no_mons
                for k = 1:prob.Fs{i,j}.nvals
                    if (prob.Fs{i,j}.isnum(k))
                        continue;
                    end
                    
                    if (~mod(l, 6))
                        fprintf(fid, '\n');
                    end
                    
                    if (l < no_con_vars)
                        fprintf(fid, '  {%d, %d, %d},', i - 1, j - 1, k - 1);
                    else
                        fprintf(fid, '  {%d, %d, %d}\n};\n', i - 1, j - 1, k - 1);
                    end
                    l = l + 1;
                end
            end
        end
    else
        % Bogus value to fool msvc
        fprintf(fid, 'const %s::conVar %s::con_vars[] = {-1};\n', class_name, class_name);
    end
end

function export_convars(fid, class_name, con_vars, no_aux)
    fprintf(fid, '\nvoid %s::evalSdpConstraintsVariables(const double *pars, const double *vmul, double **vars) const {\n', class_name); 
    fprintf(fid, '\n  double aux[%d];\n', no_aux + 1); 
    fprintf(fid, '\n  updateAuxiliaries(pars, aux);\n\n');    

    if (isempty(con_vars))
        fprintf(fid, '}\n'); 
        return;
    end    
    
    % generate C code
    [fid_c, tmp_file] = genccode(con_vars, 1);        
    
    c = 1;
    line = fgets(fid_c);
    while (ischar(line))
        [start_index, end_index] = regexp(line, 'aux(\d+)');
        index = cell(size(start_index));
        for j = 1:length(start_index)
            index{j} = line((start_index(j)+3):end_index(j));
        end
        for j = 1:length(start_index)
            line = regexprep(line, sprintf('aux(%s)', index{j}), sprintf('aux[%d]', str2double(index{j})));
        end
        line = regexprep(line, '\s+(t\d+ = .*)', '  double $1');
        line = regexprep(line, '\s+A0\[(\d+)\]\[\d+\] = ([^;]*)(;.*)', '  *(vars[$1]) = vmul[$1] * ($2)$3');
        
        if (mod(c, 4))
            line = regexprep(line, '(\n|\c)+', '');
        end
        c = c + 1;
        
        fprintf(fid, line);         
        line = fgets(fid_c);
    end
    
    fdelete(fid_c, tmp_file);    
 
    if (mod(c - 1,6))
        fprintf(fid, '\n');   
    end
    
    fprintf(fid, '}\n');   
end

function export_tvars(fid, class_name, prob)
    fprintf(fid, '\nvoid %s::updateAuxiliaries(const double *pars, double *aux) const {\n', class_name); 

    tvars_exp = values(prob.tvars_map);
    tmp = [tvars_exp{:}]';
    % tmp = [];
    % for i = 1:length(tvars_exp)
    %     tmp = [tmp, prob.tvars_map(sprintf('aux%d', i-1))];
    % end
    nppars = getsyms('pars%d', numel(prob.ppars)); 
    tmp = subs(tmp, prob.ppars, nppars);

    if (isempty(tmp))
        fprintf(fid, '}\n'); 
        return;
    end
    
    % generate C code
    [fid_c, tmp_file] = genccode(tmp, 1);    
    
    line = fgets(fid_c);
    while (ischar(line))
        line = regexprep(line, 'pars(\d+)', 'pars[$1]');        
        line = regexprep(line, '\s+(t\d+ = .*)', '  double $1');
        line = regexprep(line, '\s+A0\[(\d+)\]\[\d+\] = (.*)', '  aux[$1] = $2');
        
        fprintf(fid, line);        
        line = fgets(fid_c);
    end
    
    fdelete(fid_c, tmp_file);
    
    fprintf(fid, '}\n'); 
end

function export_polycons(fid, class_name, no_aux, prob)
    fprintf(fid, '\nvoid %s::evalPolyConstraints(const double *vars, const double *pars, double *cons) const {\n', class_name); 
    
    fprintf(fid, '\n  double aux[%d];\n', no_aux + 1); 
    fprintf(fid, '\n  updateAuxiliaries(pars, aux);\n\n'); 

    dims = cellfun(@(x) size(x, 1), prob.ccons, 'UniformOutput', true);
    dims1 = find(dims == 1);
    dims2 = find(dims ~= 1);
    
    if (isempty(dims))
        fprintf(fid, '}\n'); 
        return;
    end
    
    nvars = getsyms('vars%d', numel(prob.vars));
            
    % generate C code
    if (~isempty(dims1))
        
        tmp = prob.ccons(dims == 1);
        ccons = subs(tmp{1}, prob.vars, nvars);
        
        [fid_c, tmp_file] = genccode(ccons, 1);
        
        line = fgets(fid_c);
        c = 1;
        while (ischar(line))
            [start_index, end_index] = regexp(line, 'aux(\d+)');
            index = cell(size(start_index));
            for j = 1:length(start_index)
                index{j} = line((start_index(j)+3):end_index(j));
            end
            for j = 1:length(start_index)
                line = regexprep(line, sprintf('aux(%s)', index{j}), sprintf('aux[%d]', str2double(index{j})));
            end
            line = regexprep(line, 'vars(\d+)', 'vars[$1]');
            line = regexprep(line, '\s+(t\d+ = .*)', '  double $1');
            if (~isempty(regexp(line, '\s+A0.*', 'match')))
                line = regexprep(line, '\s+A0\[\d+\]\[(\d+)\] = (.*)', ['  cons[' sprintf('%d', dims1(c) - 1) '] = $2']);
                c = c + 1;
            end
            
            fprintf(fid, line);
            line = fgets(fid_c);
        end
        
        fdelete(fid_c, tmp_file);
    end
    
    if (~isempty(dims2)) 
        fprintf(fid, '  MatrixXd M;\n');
        
        for i = dims2
            fprintf(fid, '\n  M.resize(%d, %d);\n', dims(i), dims(i));    
            
            M = prob.ccons{i};
            M = subs(M(:), prob.vars, nvars);
            [fid_c, tmp_file] = genccode(M(:), 1);

            line = fgets(fid_c);
            c = 1;
            while (ischar(line))
                [start_index, end_index] = regexp(line, 'aux(\d+)');
                index = cell(size(start_index));
                for j = 1:length(start_index)
                    index{j} = line((start_index(j)+3):end_index(j));
                end
                for j = 1:length(start_index)
                    line = regexprep(line, sprintf('aux(%s)', index{j}), sprintf('aux[%d]', str2double(index{j})));
                end
                line = regexprep(line, 'vars(\d+)', 'vars[$1]');
                line = regexprep(line, '\s+(t\d+ = .*)', '  double $1');
                if (~isempty(regexp(line, '\s+A0.*', 'match')))
                    line = regexprep(line, '\s+A0\[(\d+)\]\[\d+\] = (.*)', '  *(M.data() + $1) = $2');     
                    c = c + 1;
                end

                fprintf(fid, line); 
                line = fgets(fid_c);
            end

            fdelete(fid_c, tmp_file);

            fprintf(fid, '  cons[%d] = M.eigenvalues().real().minCoeff();\n', i - 1);
        end
    end
    
    fprintf(fid, '}\n'); 
end

function export_polyobj(fid, class_name, prob)    
    fprintf(fid, '\ndouble %s::evalPolyObjective(const double *rpars, const double *ppars, const double *vars) const {\n', class_name); 

    nrpars = getsyms('rpars%d', numel(prob.rpars)); 
    nppars = getsyms('ppars%d', numel(prob.ppars)); 
    nvars = getsyms('vars%d', numel(prob.vars)); 

    obj = prob.obj;
    obj = subs(obj, prob.rpars, nrpars);
    obj = subs(obj, prob.ppars, nppars);
    obj = subs(obj, prob.vars, nvars);
    
    % generate C code
    [fid_c, tmp_file] = genccode(obj);
    
    line = fgets(fid_c);
    while (ischar(line))
        line = regexprep(line, 'rpars(\d+)', 'rpars[$1]');
        line = regexprep(line, 'ppars(\d+)', 'ppars[$1]');
        line = regexprep(line, 'vars(\d+)', 'vars[$1]');        
                
        if ((numel(line) > 150) && strcmp(line(1:7), '  t0 = '))
            split_line(line(3:end), fid, '  double ');
        else        
            fprintf(fid, '  double %s', line(3:end));
        end
        
        line = fgets(fid_c);
    end
    
    fdelete(fid_c, tmp_file);    
    
    fprintf(fid, '\n  return t0;\n'); 
    fprintf(fid, '}\n'); 
end

function [fid, tmp_file] = genccode(code, addline)
    tmp_file = [tempname '_gensolver'];
        
    ccode(code, 'file', tmp_file);
       
    if (numel(code) == 1 && nargin > 1)
        [fid, msg] = fopen(tmp_file, 'a');  

        if (fid == -1)
            error(['Cannot open file (' tmp_file '): ' msg]);
        end
        
        fprintf(fid, '  A0[0][0] = t0;\n');
        fclose(fid);
    end
    
    [fid, msg] = fopen(tmp_file, 'r');

    if (fid == -1)
        error(['Cannot open file (' tmp_file '): ' msg]);
    end
end

function split_line(line, fid, prefix)
    [t1, t2] = regexp(line, '.{80,100}(\+|\*|-)','match', 'split');

    fprintf(fid, prefix);

    if (numel(t1))
        fprintf(fid, '%s%s\n', t1{1}, t2{1});
    else
        fprintf(fid, '%s\n', t2{1});        
    end
    
    for i = 2:numel(t1)
        fprintf(fid, '         %s%s\n', t1{i}, t2{i});
    end

    if (numel(t2) > numel(t1) && (~strcmp(t2{end}, '')))
        fprintf(fid, '         %s', t2{end});
    end
end

function fdelete(fid, fname)
    fclose(fid);

    if (exist(fname, 'file'))
        delete(fname);
    end
end

function svars  = getsyms(mask, num)
    svars = sym(zeros(1, num));
    
    for i = 1:num
        svars(i) = sym(sprintf(mask, i - 1), 'real');
    end
end
