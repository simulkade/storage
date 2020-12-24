function C = read_input(input_file)
    %read_input(path_to_input) reads an input file for the prototype energy storage function
    % and cleans it up by removing comments and extra spaces, and returns a clean table 
    % to be used in matlab subsurface energy storage simulator 
    fid = fopen(input_file, 'r');
    C = textscan(fid, '%s','Delimiter',''); % read the input file into a cell array
    fclose(fid);
    C = C{:}; % get rid of cells in cell C

    % clean the comments:
    ind_comment = cellfun(@(x)(x(1)=='#'), C);
    C(ind_comment) = []; % remove all the comment lines
    % remove comments from other lines
    for i=1:length(C)
        C{i}(strfind(C{i}, '#'):end) = '';
        C{i} = strtrim(C{i});
    end

end % end function