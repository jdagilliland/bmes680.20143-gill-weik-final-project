function metadata = extract_meta(headers, regex, tok_no)
    if ~iscellstr(headers)
        disp('You have not provided a cell array of strings')
        headers_string = cell(size(headers));
        for i = 1:numel(headers)
            headers_string{i} = char(headers{i});
        end
        headers = headers_string;
    end
    [mat, tok] = regexp(headers, regex, 'match', 'tokens');
    metadata = cell(size(tok));
    for i = 1:numel(tok)
        metadata{i} = tok{i}{1}{tok_no};
    end
end
