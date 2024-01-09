function annotations_all_table = annotations2supplement(annotations_all)
%% turn into table


annotations_all_table = struct2table(annotations_all);


%% remove superfluous variables before writing to supplement.


annotations_all_table = removevars(annotations_all_table, ["muts_hasmutation","WARNING","locustag_CDHITcluster","locustag1_CDHITcluster","locustag2_CDHITcluster"]);
annotations_all_table = removevars(annotations_all_table, "IsolateNamesWithMutation");
annotations_all_table = removevars(annotations_all_table, "codons");
annotations_all_table = removevars(annotations_all_table, "text");


%% turn into space delimited strings to export as csv while keeping values in the same order


annotations_all_table.protein2(cellfun(@isempty,annotations_all_table.protein2))={''};
annotations_all_table.protein2 = arrayfun(@(x) string(x{:}), annotations_all_table.protein2, 'UniformOutput', false);
annotations_all_table.protein2 =arrayfun(@(x) strjoin(x{:},' '), annotations_all_table.protein2);
annotations_all_table.protein1(cellfun(@isempty,annotations_all_table.protein1))={''};;
annotations_all_table.protein1 = arrayfun(@(x) string(x{:}), annotations_all_table.protein1, 'UniformOutput', false);
annotations_all_table.protein1 =arrayfun(@(x) strjoin(x{:},' '), annotations_all_table.protein1);
annotations_all_table.note = arrayfun(@(x) string(x{:}), annotations_all_table.note, 'UniformOutput', false);
annotations_all_table.note =arrayfun(@(x) strjoin(x{:},' '), annotations_all_table.note);
annotations_all_table.note = arrayfun(@(x) string(x{:}), annotations_all_table.note, 'UniformOutput', false);
annotations_all_table.note =arrayfun(@(x) strjoin(x{:},' '), annotations_all_table.note);


%% Add the number (row)


annotations_all_table.mutation_number = (1:size(annotations_all_table,1))';
annotations_all_table.mutation_number = (1:size(annotations_all_table,1))';