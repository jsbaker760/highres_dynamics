function plot_minimum_transmission_n_isolates(TransmissionTable)

NisoSource = arrayfun(@(x) {repmat(TransmissionTable.NisolatesSource{x},numel(TransmissionTable.NisolatesRecipient{x}),1)}, 1:size(TransmissionTable,1));
NisoSource=vertcat(NisoSource{:});
NisoRecip = vertcat(TransmissionTable.NisolatesRecipient{:});

Both = [NisoSource NisoRecip];
MinBoth = min(Both,[],2);
Ntrans =vertcat(TransmissionTable.Ntransmissions{:});
X=MinBoth;Y=Ntrans;

scatter(X,Y,'filled','MarkerFaceColor','white','MarkerEdgeColor','black')
hold on
ft1 = fit(X,Y,'poly1');
[R,pv]=corr(X,Y);
eq = ['y=' char(string(round(ft1.p1,2,'significant'))) 'x+' char(string(round(ft1.p2,2,'significant')))];
plot(ft1)
l=legend;
l.String={eq};
text(80,9,['R^2=' char(string(round(R^2,2,'significant')))])
text(80,7,['P=' char(string(round(pv,2,'significant')))])
xlabel('min numisolates between source and recipient')
ylabel('N genotypes transmitted')
xlim([0 100])
ylim([0 12])
