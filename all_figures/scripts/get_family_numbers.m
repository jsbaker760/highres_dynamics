function FamilyNumbers = get_family_numbers(SIDs)

if isempty(SIDs)
    FamilyNumbers=[];
else

    [i, j] = regexpi(SIDs,'[0-9]+[PAU]');
    if iscell(i)
        FamilyNumbers = arrayfun(@(x) {char(x)}, SIDs);
        FamilyNumbers = arrayfun(@(x) str2double(FamilyNumbers{x}(i{x}:j{x}-1)), 1:numel(i));
    else
        FamilyNumbers = arrayfun(@(x) {char(x)}, SIDs);
        FamilyNumbers = arrayfun(@(x) str2double(FamilyNumbers{x}(i(x):j(x)-1)), 1:numel(i));
    end

end
end