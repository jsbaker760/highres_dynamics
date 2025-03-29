function cut_clade_assemblies(T)

%% load inputs
SamplesCSV = readtable(T);
Clades = string(unique(SamplesCSV.ReferenceGenome));

%% initialize variables

yaxisfactor =5;
bw = .05;

Ratios = zeros(size(Clades));
AssemblySize = zeros(size(Clades));
CutoffHeights = zeros(size(Clades));
MinimumContigCoverage = zeros(size(Clades));


for s =1:numel(Clades)

    CladeName = Clades(s);
    % gets the header of each entry in the assembly fasta to be parsed
    OldGenome = ['data/Assembly/clades/' char(CladeName) '/contigs.fasta'];
    if isfile(OldGenome)
        % puts the data in the cells
        [Header, Sequence] = fastaread(OldGenome);
        [NodeNumber, Coverage,Lengths]  = parseSpadesHeader(Header);
        if ~isempty(Lengths)&~isempty(Coverage)
            Longest = max(Lengths);
            A = sum(Lengths);
            LongestCov = max(Coverage(Lengths==Longest));
            Ratios(s) = max(Coverage)/LongestCov;
            AssemblySize(s) = sum(Lengths);
            %
            f=figure(s);

            p1=subplot(1,4,1);hold on
            scatter(Lengths, Coverage,'filled')
            xlabel('Lengths');
            ylabel('Coverage');

            p2=subplot(1,4,2);hold on
            scatter(Lengths, Coverage,'filled')
            xlabel('Lengths');
            ylabel('Coverage');
            ylim([0 LongestCov*yaxisfactor])
            xlim([0 Longest])

            p3=subplot(1,4,3);hold on
            xlabel('N contigs');
            ylabel('coverage');
            Bins = 0:bw:round(max(Coverage)+.5);
            BinCenters = (bw/2):bw:round(max(Coverage)+.5);
            H = histcounts(Coverage,Bins);
            GoodBin = H~=0;
            barh(BinCenters(GoodBin),H(GoodBin),'FaceColor',[0 0 1],'EdgeAlpha',0,'BarWidth',1);ylim([0 1])




            IsEmptyBin = ~GoodBin;
            Gaps = [];
            w = 0;
            for g = 1:numel(IsEmptyBin)
                if IsEmptyBin(g)
                    w = w+1;
                elseif ~IsEmptyBin(g) & w > 0
                    Gaps(end+1) = w;
                    w = 0;
                end
            end


            Height2draw = ceil(LongestCov*yaxisfactor);
            nBins2draw = sum((BinCenters<Height2draw)&IsEmptyBin);
            EmptyBinIdx = find(IsEmptyBin);

            X = BinCenters(EmptyBinIdx(1:nBins2draw));
            Y = max(H)*ones(size(X));
            barh(X,Y,'FaceColor',[.8 .8 .8],'EdgeAlpha',0,'BarWidth',1);

            ylim([0 LongestCov*yaxisfactor])
            xlim([0 max(H)])
            hold off;


            p4=subplot(1,4,4);hold on
            X = 1:max(Coverage);
            Y = arrayfun(@(x) sum(Lengths(Coverage>x)), X);
            xlabel('N contigs');
            ylabel('coverage');
            Bins = 0:bw:round(max(Coverage)+.5);
            BinCenters = (bw/2):bw:round(max(Coverage)+.5);
            H = histcounts(Coverage,Bins);
            GoodBin = H~=0;
            plot(Y,X,'Color','black','LineWidth',2);
            plot([2.5e6 2.5e6],[1 X(end)],'LineWidth',2,'Color','red');
            ylim([0 LongestCov*yaxisfactor]);
            xlabel('sum Lengths(cov>=x)')
            %%
            MinimumContigCoverage(s)=min(Coverage);
            fprintf([char(string(min(Coverage))) '\n'])
            close all
            %mkdir(['data/Assembly/clades_filtered_contigs/' char(CladeName) ]);
            %fastawrite(['data/Assembly/clades_filtered_contigs/' char(CladeName) '/genome.fasta'],Header(Coverage>CutoffHeights(s)),Sequence(Coverage>CutoffHeights(s)));

        end
    end
end
end
%%

function [NodeNumber, Coverage,Length]  = parseSpadesHeader(Header)
if ischar(Header)&size(Header,1)==1
    Header = {Header};
else
    Header = string(Header');
end

Header = arrayfun(@(x){strsplit(x{:},'_')}, Header);

NodeNumber = arrayfun(@(x) str2double(x{:}{2}),Header);
Length = arrayfun(@(x) str2double(x{:}{4}),Header);
Coverage = arrayfun(@(x) str2double(x{:}{6}),Header);
end
