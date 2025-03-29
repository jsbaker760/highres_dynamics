function ManuscriptColors = load_manuscript_colors(ColorBlind)
% Changing these values will change the colors of the resulting figures
ManuscriptColors = struct;
ManuscriptColors.PhylotypesCacnes = flip(cbrewer2('BrBG',8),1);
ManuscriptColors.PhylotypesSepi = flip(cbrewer2('PiYG',4),1);
ManuscriptColors.SpeciesCacnes = [1 1 1];
ManuscriptColors.SpeciesSepi = [0 0 0];
ManuscriptColors.AgeColorMap = 'Blues';
ManuscriptColors.CutotypeColors=[ColorBlind.LightBlue; ColorBlind.Blue ; ColorBlind.Vermillion];
ManuscriptColors.GreyScale = [1 1 1; .85 .85 .85; .7 .7 .7; .55 .55 .55; .4 .4 .4; .25 .25 .25];