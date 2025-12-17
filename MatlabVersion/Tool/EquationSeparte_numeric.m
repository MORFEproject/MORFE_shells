%% EquationSeparte_numeric.m
% This script performs the following tasks:
% 1. Reads state space equations from the text file "state_space_eq.txt" (e.g., "z1' = ...").
% 2. Extracts each term from the right-hand side polynomial and classifies them as constant, linear, or nonlinear.
% 3. Converts each term into a numeric vector of the form:
%       [coefficient, exp_z1, exp_z2, ..., exp_zN]
%    where N is the maximum index among the state variables in the system.
% 4. For each nonlinear term, the script computes its partial derivatives with respect to each state variable,
%    i.e. the tangent contributions required for the nonlinear internal force tangent matrix.
% 5. Saves the numeric results (and the tangent contributions) to both a MAT file and a text file for further numerical computations.

function EquationSeparte_numeric(FileName)
    inFilename = FileName;
    fid = fopen(inFilename, 'r');
    if fid == -1
        error('Unable to open file %s', inFilename);
    end
    
    lines = {}; 
    while ~feof(fid)
        tline = fgetl(fid);
        if ischar(tline) && ~isempty(strtrim(tline))
            lines{end+1} = tline;  
        end
    end
    fclose(fid);
    
    if isempty(lines)
        warning('The input file is empty!');
    end
    
    %%
    eqStruct = struct('derivative', {}, 'terms', {});
    for iLine = 1:length(lines)
        line = lines{iLine};
        disp(line)
        % Process only lines containing "="
        if contains(line, '=')
            parts = strsplit(line, '=');
            lhs = strtrim(parts{1});
            rhs = strtrim(parts{2});
            
            derivVar = regexprep(lhs, '''', '');
            
            % Preprocess the RHS:
            % 1. Remove all spaces
            rhs_nospace = regexprep(rhs, '\s+', '');
            % 2. Replace occurrences of "+-" or "-+" with "-"
            rhs_nospace = strrep(rhs_nospace, '+-', '-');
            rhs_nospace = strrep(rhs_nospace, '-+', '-');
            
            termPattern = '([+\-]?(?:\d+(?:\.\d+)?|\.\d+)(?:e[+\-]?\d+)?(?:\*z\d+\^\d+)*)';
            termMatches = regexp(rhs_nospace, termPattern, 'match');
            termMatches = termMatches(~cellfun(@isempty, termMatches));
            
            if isempty(termMatches)
                warning('No terms matched in line %d! Please check the format: %s', iLine, line);
            end
            
            % Initialize structure array for terms
            eqTerms = struct('coefficient', {}, 'factors', {}, 'totalDegree', {}, 'termType', {});
            for iTerm = 1:length(termMatches)
                termStr = termMatches{iTerm};
                % Split the term by "*" 
                partsTerm = strsplit(termStr, '*');
                coeffStr = partsTerm{1};
                coeff = str2double(coeffStr);
                if isnan(coeff)
                    warning('Unable to convert %s to a number!', coeffStr);
                end
                
                varList = {};   
                totalDegree = 0;
                for iFactor = 2:length(partsTerm)
                    factor = partsTerm{iFactor};
                    % Use regex to extract variable name and exponent (e.g., "z1^1")
                    token = regexp(factor, '(z\d+)\^(\d+)', 'tokens');
                    if ~isempty(token)
                        token = token{1};
                        varName = token{1};
                        expVal = str2double(token{2});
                        totalDegree = totalDegree + expVal;
                        varList{end+1} = struct('var', varName, 'exp', expVal);
                    else
                        warning('Unable to parse factor %s (in line %d)', factor, iLine);
                    end
                end
                
                % Classify the term:
                % No factors: constant term;
                % Single factor with exponent 1: linear term;
                % Otherwise: nonlinear term.
                if isempty(varList)
                    termType = 'constant';
                elseif (length(varList)==1) && (varList{1}.exp == 1)
                    termType = 'linear';
                else
                    termType = 'nonlinear';
                end
                
                eqTerms(end+1) = struct('coefficient', coeff, 'factors', {varList}, ...
                                         'totalDegree', totalDegree, 'termType', termType);
            end
            
            % Save the equation information into eqStruct
            eqStruct(end+1).derivative = derivVar; 
            eqStruct(end).terms = eqTerms;
        end
    end
    
    if isempty(eqStruct)
        warning('No equations with "=" were parsed. Please check the input file format.');
    end
    
    %%
    nVars = 0;
    for i = 1:length(eqStruct)
        for j = 1:length(eqStruct(i).terms)
            factors = eqStruct(i).terms(j).factors;
            for k = 1:length(factors)
                varStr = factors{k}.var;
                varNum = sscanf(varStr, 'z%d');
                if ~isempty(varNum)
                    nVars = max(nVars, varNum);
                end
            end
        end
    end
    
    % Build a new structure where each equation stores three matrices: constant, linear, and nonlinear.
    % Each row in these matrices is of the form: [coefficient, exp_z1, exp_z2, ..., exp_zN]
    eqStructNumeric = struct('derivative', {}, 'constant', {}, 'linear', {}, 'nonlinear', {}, 'nVars', {});
    for i = 1:length(eqStruct)
        % Initialize matrices
        constTerms = [];
        linearTerms = [];
        nonlinearTerms = [];
        
        for j = 1:length(eqStruct(i).terms)
            term = eqStruct(i).terms(j);
            % Create a row vector of length nVars+1
            numRow = zeros(1, nVars+1);
            numRow(1) = term.coefficient;
            % For each factor in the term, record the corresponding exponent
            for k = 1:length(term.factors)
                factor = term.factors{k};
                varNum = sscanf(factor.var, 'z%d');
                if ~isempty(varNum) && varNum <= nVars
                    numRow(1+varNum) = numRow(1+varNum) + factor.exp;
                end
            end
            
            % Place the row into the corresponding matrix based on the term type
            if strcmp(term.termType, 'constant')
                constTerms = [constTerms; numRow]; 
            elseif strcmp(term.termType, 'linear')
                linearTerms = [linearTerms; numRow]; 
            elseif strcmp(term.termType, 'nonlinear')
                nonlinearTerms = [nonlinearTerms; numRow]; 
            end
        end
        
        eqStructNumeric(i).derivative = eqStruct(i).derivative;
        eqStructNumeric(i).nVars = nVars;
        eqStructNumeric(i).constant = constTerms;
        eqStructNumeric(i).linear = linearTerms;
        eqStructNumeric(i).nonlinear = nonlinearTerms;
    end
    
    %% 
    for i = 1:length(eqStructNumeric)
        % Initialize a cell array for tangent contributions: one cell per state variable
        tangentContrib = cell(1, nVars);
        for j = 1:nVars
             tangentContrib{j} = []; 
        end
        nonlinearTerms = eqStructNumeric(i).nonlinear; % size: [terms x (nVars+1)]
        for r = 1:size(nonlinearTerms, 1)
             termRow = nonlinearTerms(r, :); 
             c = termRow(1);
             exponents = termRow(2:end);
             for j = 1:nVars
                 if exponents(j) > 0
                    newCoeff = c * exponents(j);
                    newExponents = exponents;
                    newExponents(j) = newExponents(j) - 1;
                    newRow = [newCoeff, newExponents];
                    tangentContrib{j} = [tangentContrib{j}; newRow];
                 end
             end
        end
        eqStructNumeric(i).nonlinearTangent = tangentContrib;
    end
    
    %% 
    matFilename = 'Output\classified_eq_numeric.mat';
    save(matFilename, 'eqStructNumeric');
    disp(['Numeric structure saved to file: ', matFilename]);
    
    txtFilename = 'Output\classified_eq_numeric.txt';
    fidOut = fopen(txtFilename, 'w');
    if fidOut == -1
        error('Unable to write to file %s', txtFilename);
    end
    
    fprintf(fidOut, 'Numeric representation of state space equations\n');
    fprintf(fidOut, 'Total number of state variables (nVars): %d\n\n', nVars);
    
    for i = 1:length(eqStructNumeric)
        fprintf(fidOut, 'Equation derivative: %s\n', eqStructNumeric(i).derivative);
        
        % Output constant terms
        fprintf(fidOut, '  Constant terms (each row format: [coefficient, exp_z1, exp_z2, ..., exp_z%d]):\n', nVars);
        if isempty(eqStructNumeric(i).constant)
            fprintf(fidOut, '    None\n');
        else
            for r = 1:size(eqStructNumeric(i).constant,1)
                fprintf(fidOut, '    ');
                fprintf(fidOut, '%g ', eqStructNumeric(i).constant(r, :));
                fprintf(fidOut, '\n');
            end
        end
        
        % Output linear terms
        fprintf(fidOut, '  Linear terms (each row format: [coefficient, exp_z1, exp_z2, ..., exp_z%d]):\n', nVars);
        if isempty(eqStructNumeric(i).linear)
            fprintf(fidOut, '    None\n');
        else
            for r = 1:size(eqStructNumeric(i).linear,1)
                fprintf(fidOut, '    ');
                fprintf(fidOut, '%g ', eqStructNumeric(i).linear(r, :));
                fprintf(fidOut, '\n');
            end
        end
        
        % Output nonlinear terms
        fprintf(fidOut, '  Nonlinear terms (each row format: [coefficient, exp_z1, exp_z2, ..., exp_z%d]):\n', nVars);
        if isempty(eqStructNumeric(i).nonlinear)
            fprintf(fidOut, '    None\n');
        else
            for r = 1:size(eqStructNumeric(i).nonlinear,1)
                fprintf(fidOut, '    ');
                fprintf(fidOut, '%g ', eqStructNumeric(i).nonlinear(r, :));
                fprintf(fidOut, '\n');
            end
        end
        
        % Output tangent contributions from nonlinear terms
        fprintf(fidOut, '  Nonlinear tangent contributions (each cell corresponds to partial derivative with respect to a state variable):\n');
        for j = 1:nVars
            fprintf(fidOut, '    d/dz%d:\n', j);
            if isempty(eqStructNumeric(i).nonlinearTangent{j})
                fprintf(fidOut, '      None\n');
            else
                for r = 1:size(eqStructNumeric(i).nonlinearTangent{j},1)
                    fprintf(fidOut, '      ');
                    fprintf(fidOut, '%g ', eqStructNumeric(i).nonlinearTangent{j}(r, :));
                    fprintf(fidOut, '\n');
                end
            end
        end
        fprintf(fidOut, '\n');
    end
    fclose(fidOut);
    disp(['Numeric representation results written to file: ', txtFilename]);

end