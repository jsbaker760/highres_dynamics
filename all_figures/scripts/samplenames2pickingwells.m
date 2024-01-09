function [pick_plate, pick_well] = samplenames2pickingwells(SampleNames)
%%

x_wells = readtable('x_plate_picking_wells.csv','NumHeaderLines',0);
all_sample_names = readtable('samplenames.csv','NumHeaderLines',0);
x_wells = string(table2array(x_wells(2:end,:)));

%%

[pick_plate, pick_well] = deal(strings(numel(SampleNames),1));

for i = 1:numel(SampleNames)
    SN = SampleNames(i);
    % if this samplename contains an X or Z
    if contains(SN,"Z")||contains(SN,"X")
        splt = strsplit(SN,"_");
        if numel(splt)==2
            splt = char(splt(2));
            well = str2double(string(splt(end-1:end)));
        elseif numel(splt)~=2&~isempty(regexp(SN,'X[0][1-6]','match'))
            well = str2double(string(splt(3)));
            splt = char(splt(2));
        elseif numel(splt)~=2&(~isempty(regexp(SN,'X[0][7-9]','match'))|~isempty(regexp(SN,'X[1][0-9]','match')))
            well = str2double(string(splt(3)));
            splt = char(splt(2));
        elseif ~contains(SN,"Z")
            error('not found')
        end
        if contains(SN,"Z")
            idx = find(string(all_sample_names.SampleName)==SN);
            if isempty(idx)
                error('not found')
            end
            old_name = string(all_sample_names.OLD_SampleName(idx));
            old_name=strsplit(old_name,"_");
            pick_plate(i)=old_name(4);
            c = char(old_name(5));
            if numel(c)==3
                pick_well(i)=str2double(string(c(end)));
            elseif numel(c)==4||numel(c)==5
                pick_well(i)=str2double(string(c(end-1:end)));
            else
                error('not found')
            end
        elseif contains(SN,"X")
            XNumber = str2double(string(strrep(splt,'X','')));
            if XNumber>1700
                well=XNumber-1700;
                XNumber=17;
            end 
            if XNumber<=6
                PlateStartI = (8*(XNumber-1))+1;
                [jj,ii]=ind2sub([12 8],well);
                ii = (PlateStartI-1)+ii;
                pickname = x_wells(ii,jj);
                out = regexp(pickname,'[EAP][FNKC](-r)?','match');
                if isempty(out)
                    out = regexp(pickname,'[P][123]','match');
                end
                if ~isempty(out)
                    s=strsplit(pickname,"-");
                    pick_plate(i)=out;
                    pick_well(i)=s(end);
                    continue
                else
                    error('not found')
                end
            elseif XNumber>6&XNumber<17
                PlateStartI = (8*(XNumber-1))+1;
                [jj,ii]=ind2sub([12 8],well);
                ii = (PlateStartI-1)+ii;
                pickname = x_wells(ii,jj);
                out = regexp(pickname,'S[1-5](-)?[R12345]([1-4])?([ABCD])?','match');
                if ~isempty(out)
                    pick_plate(i)=out;
                    well =regexp(pickname,'[ABCDEFGH][1-9]([0-9])?','match');
                    if ~isempty(well)
                        pick_well(i)=well;
                    else
                        error('not found')
                    end
                else
                    error('not found')
                end
            elseif XNumber>=17&XNumber<20
                pick_well(i)=well;
                pick_plate(i)=string(splt);
            else
                error('not found')
            end
        end
    % for all other samples, the well it was picked in is already in the
    % name
    else
        s = strsplit(SN,"_");
        if numel(s)~=3
            error('not found')
        else
            pick_well(i)=s(3);
            pick_plate(i)=s(2);
        end
    end
end